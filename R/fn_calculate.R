#' Score patches from the downloaded inputs (calculate module)
#'
#' The calculation half of the patch ranking: given the data acquired by
#' \code{\link{rank_download}} plus the on-disk surfaces, it computes every
#' \code{RankIndex} term, samples the raster terms at the patch centroids, and combines
#' them with \code{\link{rank_index}}. No network calls — it reads the files the
#' download module wrote, so it is fast to re-run while tuning weights.
#'
#' Terms are computed only when their inputs are supplied, so the ranking degrades
#' gracefully: omit \code{records} to drop \strong{S}, omit \code{network}/\code{origin}
#' to drop \strong{A}/\strong{C}. \strong{Dg}/\strong{Dr} need only \code{occ} and
#' \code{suitability}.
#'
#' @param patches \pkg{sf} units to score (polygons or points); must carry
#'   \code{area_col} for \strong{C}.
#' @param occ known-population points (\pkg{sf}/\code{SpatVector}) — the sources for
#'   \strong{Dg} and \strong{Dr}.
#' @param suitability habitat-suitability \code{SpatRaster} (or path): the template for
#'   the raster terms and, inverted, the dispersal resistance.
#' @param records GBIF \pkg{sf} with \code{speciesKey} (and \code{basisOfRecord}) for
#'   \strong{S}; \code{NULL} skips it.
#' @param network,origin OSM road/trail \pkg{sf} and the trip origin for \strong{A}/
#'   \strong{C}; \code{NULL} skips them.
#' @param dem optional elevation \code{SpatRaster} for the accessibility off-trail hop.
#' @param resistance optional precomputed friction \code{SpatRaster} for \strong{Dr};
#'   \code{NULL} builds one from inverted suitability (\code{\link{build_resistance}}).
#' @param agg integer factor to coarsen \code{suitability} by before the raster terms
#'   (default \code{10}); the ranking is sampled at centroids, so fine resolution is
#'   wasted here.
#' @param area_col patch-area column (acres) used by \strong{C} (default \code{"acres"}).
#' @param sigma,target_area,k,w_habitat term parameters passed through to
#'   \code{\link{survey_surface}} / \code{\link{complementarity}} /
#'   \code{\link{accessibility}} / \code{\link{build_resistance}}.
#' @param weights per-term multipliers for \code{\link{rank_index}} (default equal).
#' @return \code{patches} with the rescaled term columns, \code{index} and \code{rank}
#'   appended (an \pkg{sf}).
#' @seealso \code{\link{rank_download}}, \code{\link{rank_index}}
#' @export
rank_calculate <- function(patches, occ, suitability, records = NULL,
                           network = NULL, origin = NULL, dem = NULL,
                           resistance = NULL, agg = 10, area_col = "acres",
                           sigma = 100, target_area = 50, k = 3, w_habitat = 5,
                           weights = NULL) {
  patches <- sf::st_as_sf(patches)
  suit <- if (inherits(suitability, "SpatRaster")) suitability
          else terra::rast(suitability)
  tmpl <- if (agg > 1) terra::aggregate(suit, agg, fun = "mean", na.rm = TRUE) else suit
  cen <- terra::centroids(terra::vect(patches))
  samp <- function(r) terra::extract(r, cen)[, 2]

  terms <- list()

  # Dg — geographic distance to nearest known population.
  terms$Dg <- samp(geo_distance(occ, tmpl, invert = TRUE))

  # Dr — resistance distance over inverted-suitability friction.
  fric <- if (is.null(resistance))
            build_resistance(tmpl, habitat = (1 - tmpl), w_habitat = w_habitat)
          else terra::resample(if (inherits(resistance, "SpatRaster")) resistance
                               else terra::rast(resistance), tmpl)
  terms$Dr <- samp(cost_distance(fric, occ, template = tmpl, invert = TRUE))

  # S — survey intensity (if records supplied). Evaluate the gridless field directly
  # at the centroids (no full-raster render); negative intensity so rank_index's
  # rescale makes well-surveyed low and unsurveyed high.
  if (!is.null(records)) {
    records <- sf::st_transform(sf::st_as_sf(records), sf::st_crs(tmpl))
    rc <- sf::st_coordinates(records)[, 1:2, drop = FALSE]
    bas <- if ("basisOfRecord" %in% names(records)) records$basisOfRecord else NULL
    fld <- survey_field(rc, records$speciesKey, sigma = sigma, basis = bas,
                        domain = terra::ext(tmpl))
    terms$S <- -eval_survey(fld, terra::crds(cen))
  }

  # A, C — accessibility and complementarity over the network (if supplied).
  if (!is.null(network) && !is.null(origin)) {
    acc <- accessibility(network, patches, origins = origin, dem = dem, k = k)
    terms$A <- acc$patches$A
    terms$C <- complementarity(acc$net, acc$patch_nodes,
                               areas = patches[[area_col]],
                               target_area = target_area)$C
  }

  scored <- do.call(rank_index, c(terms, list(weights = weights)))
  sf::st_sf(data.frame(sf::st_drop_geometry(patches), scored,
                       check.names = FALSE),
            geometry = sf::st_geometry(patches))
}
