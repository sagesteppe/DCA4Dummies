#' Spatial residual diagnostics and residual-kriging correction for SDM predictions
#'
#' A small pipeline for the question a pooled decision curve cannot answer: are a
#' model's errors spatially structured, and if so can that structure be turned into
#' a bias-correction surface? The stages are (1) \code{point_residuals()} — the
#' signed difference (or Pearson / deviance residual) between an observed occurrence
#' and a predicted suitability; (2) \code{residual_variogram()} — the empirical
#' variogram and an optional fitted model, with a pure-nugget flag; (3)
#' \code{residual_krige()} — ordinary kriging of the residuals onto new locations
#' (the correction surface); (4) \code{bootstrap_correction()} — a point-resampling
#' ensemble of that surface, for uncertainty; and \code{apply_residual_correction()}
#' to add the correction back onto a prediction raster. Depends on sf, gstat, terra.
#'
#' Sign convention throughout: a POSITIVE residual is a site the model UNDER-rates
#' (observed exceeds predicted); a negative residual is one it over-rates.
#'
#' The correction surface is only trustworthy where the model interpolates, so
#' \code{residual_krige()} and \code{bootstrap_correction()} take an \code{aoa} mask
#' (the area of applicability) and predict inside it only; outside the AOA no
#' correction is produced.

#' Restrict prediction locations to the area of applicability
#'
#' @param newdata An \code{sf} point layer or a \code{terra} \code{SpatRaster}.
#' @param aoa The AOA as a \code{SpatRaster} mask (in-AOA cells non-\code{NA} and
#'   non-zero) or an \code{sf}/\code{SpatVector} polygon. \code{NULL} returns
#'   \code{newdata} unchanged.
#' @param crs CRS to align a polygon \code{aoa} against (the points' CRS).
#'
#' @return \code{newdata} masked (raster) or filtered (points) to the AOA.
#' @keywords internal
restrict_to_aoa <- function(newdata, aoa, crs = sf::st_crs(newdata)) {
  if (is.null(aoa)) return(newdata)
  if (inherits(newdata, "SpatRaster")) {
    if (inherits(aoa, "SpatRaster")) {
      a <- if (!terra::compareGeom(newdata, aoa, stopOnError = FALSE))
        terra::resample(aoa, newdata, method = "near") else aoa
      terra::mask(newdata, a)
    } else {
      terra::mask(newdata, if (inherits(aoa, "sf")) terra::vect(aoa) else aoa)
    }
  } else {
    if (inherits(aoa, "SpatRaster")) {
      v <- terra::extract(aoa[[1]], terra::vect(newdata))[[2]]
      newdata[is.finite(v) & v != 0, ]
    } else {
      av <- if (inherits(aoa, "SpatVector")) sf::st_as_sf(aoa) else aoa
      inside <- lengths(sf::st_intersects(newdata, sf::st_transform(av, sf::st_crs(newdata)))) > 0
      newdata[inside, ]
    }
  }
}

#' Signed residuals between observed occurrence and predicted suitability
#'
#' @param points An \code{sf} point layer.
#' @param observed Column name in \code{points}, or a numeric 0/1 vector, of the
#'   observed occurrence.
#' @param predicted The predicted suitability: a \code{terra} \code{SpatRaster}
#'   (sampled at the points), a column name in \code{points}, or a numeric vector.
#' @param type \code{"difference"} (\eqn{y-\hat p}, the default and most
#'   interpretable), \code{"pearson"}, or \code{"deviance"}.
#' @param layer For a raster \code{predicted}, the layer to sample (name or index;
#'   default the first).
#' @param eps Clamp applied to predictions before the Pearson / deviance transforms.
#'
#' @return \code{points} with added columns \code{pred} and \code{resid}.
#' @export
point_residuals <- function(points, observed, predicted,
                            type = c("difference", "pearson", "deviance"),
                            layer = 1, eps = 1e-6) {
  type <- match.arg(type)
  y <- if (is.character(observed)) points[[observed]] else observed
  if (is.null(y)) stop("`observed` not found in `points`.")
  y <- as.numeric(y)

  if (inherits(predicted, "SpatRaster")) {
    ex <- terra::extract(predicted[[layer]], terra::vect(points))
    p <- ex[[ncol(ex)]]
  } else if (is.character(predicted)) {
    p <- as.numeric(points[[predicted]])
  } else {
    p <- as.numeric(predicted)
  }
  if (length(p) != length(y)) stop("`observed` and `predicted` lengths differ.")

  pc <- pmin(pmax(p, eps), 1 - eps)
  resid <- switch(type,
    difference = y - p,
    pearson    = (y - pc) / sqrt(pc * (1 - pc)),
    deviance   = sign(y - pc) * sqrt(-2 * (y * log(pc) + (1 - y) * log(1 - pc)))
  )
  points$pred  <- p
  points$resid <- resid
  points[is.finite(points$resid), ]
}

#' Empirical variogram of point residuals, with an optional fitted model
#'
#' @param points An \code{sf} point layer carrying a residual column.
#' @param resid Name of the residual column (default \code{"resid"}).
#' @param cutoff,width Passed to \code{\link[gstat]{variogram}}; defaults let gstat
#'   choose.
#' @param model Variogram model form for the fit, e.g. \code{"Exp"}, \code{"Sph"}.
#' @param fit Whether to fit \code{model} to the empirical variogram.
#' @param nugget_cut A field is flagged \code{pure_nugget} when the fitted nugget is
#'   at least this fraction of the total sill, or when the fit does not converge.
#'
#' @return A list with \code{empirical} (the gstat variogram data frame),
#'   \code{model} (the fitted model or \code{NULL}), \code{nugget}, \code{psill},
#'   \code{range}, \code{nugget_ratio}, \code{sample_var}, and \code{pure_nugget}.
#' @export
residual_variogram <- function(points, resid = "resid", cutoff = NULL, width = NULL,
                               model = "Exp", fit = TRUE, nugget_cut = 0.9) {
  f <- stats::as.formula(paste(resid, "~ 1"))
  args <- list(object = f, data = points)
  if (!is.null(cutoff)) args$cutoff <- cutoff
  if (!is.null(width))  args$width  <- width
  v <- do.call(gstat::variogram, args)
  sample_var <- stats::var(points[[resid]], na.rm = TRUE)

  m <- NULL; nugget <- NA_real_; psill <- NA_real_; rng <- NA_real_
  if (fit) {
    init <- gstat::vgm(psill = sample_var / 2, model = model,
                       range = max(v$dist) / 3, nugget = sample_var / 2)
    m <- tryCatch(gstat::fit.variogram(v, init), warning = function(w) NULL,
                  error = function(e) NULL)
    if (!is.null(m) && nrow(m) == 2) {
      nugget <- m$psill[1]; psill <- m$psill[2]; rng <- m$range[2]
    }
  }
  nugget_ratio <- if (is.finite(nugget) && is.finite(psill) && (nugget + psill) > 0)
    nugget / (nugget + psill) else NA_real_
  pure_nugget <- is.null(m) || is.na(nugget_ratio) || nugget_ratio >= nugget_cut

  list(empirical = v, model = m, nugget = nugget, psill = psill, range = rng,
       nugget_ratio = nugget_ratio, sample_var = sample_var, pure_nugget = pure_nugget)
}

#' Ordinary kriging of residuals onto new locations (the correction surface)
#'
#' @param points An \code{sf} point layer carrying a residual column.
#' @param newdata Where to predict: an \code{sf} point layer, or a \code{terra}
#'   \code{SpatRaster} template (non-\code{NA} cell centres are used and a raster is
#'   returned). Keep raster templates coarse — kriging every cell of a fine grid is
#'   expensive.
#' @param resid Name of the residual column (default \code{"resid"}).
#' @param model A fitted variogram model (e.g. \code{residual_variogram()$model}).
#' @param aoa Area-of-applicability mask restricting where the correction is
#'   produced: a \code{SpatRaster} mask or an \code{sf}/\code{SpatVector} polygon
#'   (see \code{restrict_to_aoa}). \code{NULL} predicts everywhere in \code{newdata}.
#' @param ... Passed to \code{\link[gstat]{krige}} (e.g. \code{nmax}).
#'
#' @return \code{newdata} with the kriged correction (\code{var1.pred}) and variance
#'   (\code{var1.var}); a \code{SpatRaster} of both when \code{newdata} is a raster.
#'   Locations outside \code{aoa} are dropped (points) or \code{NA} (raster).
#' @export
residual_krige <- function(points, newdata, resid = "resid", model, aoa = NULL, ...) {
  if (is.null(model)) stop("`model` is NULL (pure nugget?) — nothing to interpolate.")
  f <- stats::as.formula(paste(resid, "~ 1"))
  newdata <- restrict_to_aoa(newdata, aoa, sf::st_crs(points))
  ras <- inherits(newdata, "SpatRaster")
  nd <- if (ras) {
    xy <- terra::crds(newdata[[1]], na.rm = TRUE)
    sf::st_as_sf(as.data.frame(xy), coords = c("x", "y"), crs = sf::st_crs(points))
  } else newdata
  k <- gstat::krige(f, points, nd, model = model, debug.level = 0, ...)
  if (!ras) return(k)
  tmpl <- terra::deepcopy(newdata[[1]]); terra::values(tmpl) <- NA
  cells <- terra::cellFromXY(tmpl, terra::crds(newdata[[1]], na.rm = TRUE))
  pred <- var <- tmpl
  pred[cells] <- k$var1.pred; var[cells] <- k$var1.var
  out <- c(pred, var); names(out) <- c("correction", "krige_var"); out
}

#' Bootstrap ensemble of the residual-kriging correction surface
#'
#' Resamples the points with replacement \code{n} times, re-krige each draw onto
#' \code{newdata}, and summarise the ensemble. Duplicate coordinates from resampling
#' are jittered slightly so the kriging system stays non-singular.
#'
#' @param points An \code{sf} point layer carrying a residual column.
#' @param newdata Prediction locations (\code{sf} points or a \code{SpatRaster}).
#' @param resid Name of the residual column (default \code{"resid"}).
#' @param model A fitted variogram model.
#' @param n Number of bootstrap resamples.
#' @param aoa Area-of-applicability mask restricting where the correction is
#'   produced (see \code{residual_krige}). \code{NULL} predicts everywhere.
#' @param jitter Metres of jitter applied to resampled points to break ties.
#' @param seed Optional RNG seed.
#' @param ... Passed to \code{\link[gstat]{krige}}.
#'
#' @return \code{newdata} (or a \code{SpatRaster}) with the ensemble mean
#'   \code{corr_mean} and standard deviation \code{corr_sd} of the correction,
#'   restricted to \code{aoa}.
#' @export
bootstrap_correction <- function(points, newdata, resid = "resid", model, n = 199,
                                 aoa = NULL, jitter = 1, seed = NULL, ...) {
  if (is.null(model)) stop("`model` is NULL (pure nugget?) — nothing to interpolate.")
  if (!is.null(seed)) set.seed(seed)
  f <- stats::as.formula(paste(resid, "~ 1"))
  newdata <- restrict_to_aoa(newdata, aoa, sf::st_crs(points))
  ras <- inherits(newdata, "SpatRaster")
  nd <- if (ras) {
    xy <- terra::crds(newdata[[1]], na.rm = TRUE)
    sf::st_as_sf(as.data.frame(xy), coords = c("x", "y"), crs = sf::st_crs(points))
  } else newdata

  mat <- matrix(NA_real_, nrow = nrow(nd), ncol = n)
  for (b in seq_len(n)) {
    s <- points[sample.int(nrow(points), replace = TRUE), ]
    if (jitter > 0) s <- sf::st_jitter(s, amount = jitter)
    kb <- tryCatch(gstat::krige(f, s, nd, model = model, debug.level = 0, ...),
                   error = function(e) NULL)
    if (!is.null(kb)) mat[, b] <- kb$var1.pred
  }
  corr_mean <- rowMeans(mat, na.rm = TRUE)
  corr_sd   <- apply(mat, 1, stats::sd, na.rm = TRUE)

  if (!ras) {
    nd$corr_mean <- corr_mean; nd$corr_sd <- corr_sd; return(nd)
  }
  tmpl <- terra::deepcopy(newdata[[1]]); terra::values(tmpl) <- NA
  cells <- terra::cellFromXY(tmpl, terra::crds(newdata[[1]], na.rm = TRUE))
  mn <- sdv <- tmpl; mn[cells] <- corr_mean; sdv[cells] <- corr_sd
  out <- c(mn, sdv); names(out) <- c("corr_mean", "corr_sd"); out
}

#' Apply a residual correction to a prediction raster
#'
#' Adds a correction surface to predicted suitability and clamps to a valid range.
#' The correction is typically defined only inside the AOA (see
#' \code{residual_krige}); \code{outside} controls what happens to cells with no
#' correction.
#'
#' @param pred A \code{SpatRaster} of predicted suitability.
#' @param correction A \code{SpatRaster} correction (resampled to \code{pred} if their
#'   grids differ), \code{NA} outside the AOA.
#' @param clamp Range to clamp the corrected values to (default \code{c(0, 1)};
#'   \code{NULL} to skip).
#' @param outside What to do where \code{correction} is \code{NA}: \code{"keep"} the
#'   original prediction uncorrected (default), or \code{"na"} to blank those cells.
#'
#' @return A \code{SpatRaster} of corrected suitability.
#' @export
apply_residual_correction <- function(pred, correction, clamp = c(0, 1),
                                      outside = c("keep", "na")) {
  outside <- match.arg(outside)
  if (!terra::compareGeom(pred, correction, stopOnError = FALSE))
    correction <- terra::resample(correction, pred)
  out <- if (outside == "keep") pred + terra::ifel(is.na(correction), 0, correction)
         else pred + correction
  if (!is.null(clamp)) out <- terra::clamp(out, clamp[1], clamp[2])
  names(out) <- "suitability_corrected"
  out
}
