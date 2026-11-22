#' Fetch target-group occurrence records from GBIF for the survey-intensity term
#'
#' Pulls occurrences inside a bounding box from GBIF as \pkg{sf} points, ready to feed
#' \code{\link{survey_surface}} as the target-group records for \strong{S}. The pull
#' uses GBIF's server-side \strong{bulk download} (\code{rgbif::occ_download}) — one
#' prepared job rather than thousands of fragile paged requests, and the route GBIF
#' asks heavy users to take. \code{count_only = TRUE} first (no account needed) to see
#' the filtered total before committing.
#'
#' \strong{Requires a free GBIF account} for the download (not for the count). Register
#' at \url{https://www.gbif.org/}, then set \code{GBIF_USER}, \code{GBIF_PWD} and
#' \code{GBIF_EMAIL} (e.g. in \code{~/.Renviron}) or pass them as arguments. A download
#' is minted under your account and citable by its DOI.
#'
#' \strong{Filtering is the point.} A raw vascular-plant pull is swamped by imprecise
#' iNaturalist observations. \code{year_from} and \code{coord_uncertainty_max} trim to
#' records fit to mark \emph{where a careful collection happened}. The uncertainty
#' filter also drops records with \emph{no} uncertainty recorded, which is most of them
#' — so it defaults to \code{NA} (off); set e.g. \code{100} only if you want to keep
#' solely well-located records.
#'
#' @param bbox the area: an \pkg{sf}/\pkg{terra} object (extent used), or numeric
#'   \code{c(xmin, ymin, xmax, ymax)} in lon/lat.
#' @param taxon_key GBIF backbone key to restrict to (default \code{7707728},
#'   Tracheophyta / vascular plants).
#' @param basis \code{basisOfRecord} values to keep (default preserved specimens +
#'   human observations).
#' @param year_from earliest year to keep (default \code{1990}); \code{NA} disables.
#' @param coord_uncertainty_max max coordinate uncertainty in metres (default
#'   \code{NA} = off, since the field is rarely populated); a number also drops records
#'   lacking the field.
#' @param crs output CRS (default lon/lat).
#' @param count_only if \code{TRUE}, return just the filtered record count — a cheap
#'   request, no account, no download.
#' @param path GeoPackage to write the records to (\code{NULL} ->
#'   \code{gbif_records.gpkg} in the working directory).
#' @param download_dir where to save the GBIF download zip (\code{NULL} -> beside
#'   \code{path}); kept so the costly job need not be re-run.
#' @param user,pwd,email GBIF credentials (default the \code{GBIF_USER}/\code{GBIF_PWD}/
#'   \code{GBIF_EMAIL} environment variables).
#' @param poll seconds between status checks while the download is prepared (default
#'   \code{30}).
#' @return with \code{count_only}, an integer. Otherwise an \pkg{sf} of \code{POINT}s
#'   with \code{speciesKey}, \code{species}, \code{basisOfRecord},
#'   \code{coordinateUncertaintyInMeters} and \code{year}; attributes \code{path} (file
#'   written) and \code{doi} (citable download). Pass its coords + \code{speciesKey} to
#'   \code{\link{survey_surface}}.
#' @seealso \code{\link{survey_surface}}
#' @export
gbif_records <- function(bbox, taxon_key = 7707728,
                         basis = c("PRESERVED_SPECIMEN", "HUMAN_OBSERVATION"),
                         year_from = 1990, coord_uncertainty_max = NA,
                         crs = NULL, count_only = FALSE, path = NULL,
                         download_dir = NULL,
                         user = Sys.getenv("GBIF_USER"),
                         pwd = Sys.getenv("GBIF_PWD"),
                         email = Sys.getenv("GBIF_EMAIL"), poll = 30) {
  stopifnot(requireNamespace("rgbif", quietly = TRUE))
  bb <- .lonlat_bbox(bbox)
  wkt <- rgbif::gbif_bbox2wkt(bb$box[1], bb$box[2], bb$box[3], bb$box[4])
  out_crs <- if (!is.null(crs)) sf::st_crs(crs)
             else if (!is.na(bb$crs)) bb$crs else sf::st_crs(4326)

  # --- count: a cheap pre-check, no account (vector basis trips occ_search, so
  # count each stream and sum) ---
  if (count_only) {
    qargs <- list(taxonKey = taxon_key, geometry = wkt, hasCoordinate = TRUE)
    if (!is.na(year_from))
      qargs$year <- paste0(year_from, ",", as.integer(format(Sys.Date(), "%Y")))
    if (!is.na(coord_uncertainty_max))
      qargs$coordinateUncertaintyInMeters <- paste0("0,", coord_uncertainty_max)
    counts <- vapply(basis, function(b) {
      n <- do.call(rgbif::occ_search,
                   c(qargs, list(basisOfRecord = b, limit = 1)))$meta$count
      if (is.null(n)) 0L else as.integer(n)
    }, integer(1))
    return(sum(counts))
  }

  # Blocking convenience: request, wait, retrieve.
  key <- gbif_request(bbox, taxon_key = taxon_key, basis = basis,
                      year_from = year_from,
                      coord_uncertainty_max = coord_uncertainty_max,
                      user = user, pwd = pwd, email = email)
  gbif_retrieve(key, path = path, download_dir = download_dir, crs = out_crs,
                poll = poll)
}

#' Request a GBIF download (non-blocking) — stage 1 of 2
#'
#' Submits the download predicate to GBIF and returns immediately with the job key,
#' \strong{without} waiting for the server to prepare the archive. Do other work (e.g.
#' the OSM pull) while GBIF builds it, then collect with \code{\link{gbif_retrieve}}.
#' Arguments mirror \code{\link{gbif_records}}.
#'
#' @inheritParams gbif_records
#' @return the \code{occ_download} job key (pass to \code{\link{gbif_retrieve}}).
#' @seealso \code{\link{gbif_retrieve}}, \code{\link{gbif_records}}
#' @export
gbif_request <- function(bbox, taxon_key = 7707728,
                         basis = c("PRESERVED_SPECIMEN", "HUMAN_OBSERVATION"),
                         year_from = 1990, coord_uncertainty_max = NA,
                         user = Sys.getenv("GBIF_USER"),
                         pwd = Sys.getenv("GBIF_PWD"),
                         email = Sys.getenv("GBIF_EMAIL")) {
  stopifnot(requireNamespace("rgbif", quietly = TRUE))
  if (!nzchar(user) || !nzchar(pwd) || !nzchar(email))
    stop("GBIF download needs an account: set GBIF_USER, GBIF_PWD, GBIF_EMAIL ",
         "(register free at gbif.org), or pass user/pwd/email. Use count_only = TRUE ",
         "to preview counts without an account.")
  bb <- .lonlat_bbox(bbox)
  wkt <- rgbif::gbif_bbox2wkt(bb$box[1], bb$box[2], bb$box[3], bb$box[4])

  preds <- list(rgbif::pred("taxonKey", taxon_key),
                rgbif::pred("hasCoordinate", TRUE),
                rgbif::pred_in("basisOfRecord", basis),
                rgbif::pred_within(wkt))
  if (!is.na(year_from))
    preds <- c(preds, list(rgbif::pred_gte("year", year_from)))
  if (!is.na(coord_uncertainty_max))
    preds <- c(preds, list(rgbif::pred_lte("coordinateUncertaintyInMeters",
                                           coord_uncertainty_max)))

  key <- do.call(rgbif::occ_download,
                 c(preds, list(format = "SIMPLE_CSV",
                               user = user, pwd = pwd, email = email)))
  message("GBIF download requested: ", key, " — preparing server-side (non-blocking).")
  key
}

#' Retrieve a requested GBIF download — stage 2 of 2
#'
#' Waits for the \code{\link{gbif_request}} job to finish preparing (usually already
#' done if you did other work in between), downloads the archive, imports it, and
#' writes \pkg{sf} points.
#'
#' @param key the job key from \code{\link{gbif_request}}.
#' @param path,download_dir,crs,poll as in \code{\link{gbif_records}}.
#' @return an \pkg{sf} of points; attributes \code{path} and \code{doi}.
#' @seealso \code{\link{gbif_request}}
#' @export
gbif_retrieve <- function(key, path = NULL, download_dir = NULL, crs = NULL,
                          poll = 30) {
  stopifnot(requireNamespace("rgbif", quietly = TRUE))
  out_crs <- if (is.null(crs)) sf::st_crs(4326) else sf::st_crs(crs)
  rgbif::occ_download_wait(key, status_ping = poll, quiet = FALSE)

  if (is.null(path)) path <- file.path(getwd(), "gbif_records.gpkg")
  if (is.null(download_dir)) download_dir <- dirname(path)
  dir.create(download_dir, showWarnings = FALSE, recursive = TRUE)
  zf <- rgbif::occ_download_get(key, path = download_dir, overwrite = TRUE)
  occ <- rgbif::occ_download_import(zf)

  keep <- !is.na(occ$decimalLongitude) & !is.na(occ$decimalLatitude) &
          !is.na(occ$speciesKey)
  cols <- intersect(c("speciesKey", "species", "basisOfRecord",
                      "coordinateUncertaintyInMeters", "year"), names(occ))
  res <- sf::st_as_sf(as.data.frame(occ[keep, c("decimalLongitude",
                                                "decimalLatitude", cols)]),
                      coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)
  if (out_crs != sf::st_crs(4326)) res <- sf::st_transform(res, out_crs)
  sf::st_write(res, path, quiet = TRUE, delete_dsn = TRUE)
  doi <- tryCatch(as.character(attr(key, "doi")), error = function(e) NA_character_)
  doi <- if (length(doi)) doi[1] else NA_character_   # attr may be NULL / length 0
  message(sprintf("gbif_retrieve: %d records -> %s%s", nrow(res), path,
                  if (!is.na(doi)) paste0("  (DOI: ", doi, ")") else ""))
  structure(res, path = path, doi = doi)
}
