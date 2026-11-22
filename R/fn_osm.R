#' Fetch the OSM access network (roads + trails) for a domain
#'
#' A convenience fetcher for the \strong{A} (accessibility) and \strong{C}
#' (complementarity) terms: pulls the OpenStreetMap highway network inside a bounding
#' box from an Overpass endpoint and returns it as \pkg{sf} lines, ready to pass as
#' \code{sources} to \code{\link{cost_distance}}. OpenTrailMap and similar viewers are
#' just renderers of this same OSM data.
#'
#' \strong{Be a good citizen.} Public Overpass instances are a shared, volunteer-run
#' resource with a tight usage policy. This helper makes only modest promises and
#' tries to behave:
#' \itemize{
#'   \item \strong{Batches} a large domain into \code{tiles}^2 smaller queries rather
#'     than one giant one.
#'   \item \strong{Checkpoints} every tile to \code{cache_dir} the moment it returns,
#'     and on a re-run \strong{skips} any tile already cached — so an interrupted or
#'     retried pull resumes instead of re-hammering the API.
#'   \item \strong{Breathes} between calls (\code{req_throttle}) and \strong{backs off}
#'     exponentially on transient failures / rate limits (429/502/503/504, honouring
#'     any \code{Retry-After}).
#'   \item Sends a plain \strong{User-Agent} and never your identity. If you have a
#'     contact address you are willing to share per OSM etiquette, pass it in
#'     \code{contact} — it is \emph{opt-in}, never read from your session.
#' }
#' For anything heavy, run your \strong{own Overpass instance} (or a provider you have
#' an account with) and point \code{endpoint} at it — do not hammer the public API.
#'
#' Data is \eqn{\copyright} OpenStreetMap contributors, \strong{ODbL}: attribute it,
#' and snapshot the result (the \code{timestamp} attribute records the OSM base time)
#' so analyses stay reproducible as OSM changes.
#'
#' @param bbox the domain: an \pkg{sf}/\pkg{terra} object (its extent is used), or a
#'   numeric \code{c(xmin, ymin, xmax, ymax)} in lon/lat (EPSG:4326).
#' @param crs output CRS (anything \code{sf::st_crs} accepts). Defaults to the input
#'   object's CRS, else lon/lat.
#' @param roads,trails OSM \code{highway} values to treat as each class; the returned
#'   lines carry a \code{class} column ("road"/"trail"). Set either to
#'   \code{character(0)} to skip it.
#' @param path file to write the merged network to (GeoPackage). \code{NULL} (default)
#'   writes \code{osm_access.gpkg} in the working directory.
#' @param cache_dir directory for per-tile checkpoints (one \code{.rds} per tile).
#'   \code{NULL} (default) uses a folder named after \code{path} beside it. Existing
#'   tiles here are reused unless \code{overwrite = TRUE}; delete it to force a fresh
#'   pull.
#' @param overwrite refetch and overwrite cached tiles (default \code{FALSE}).
#' @param endpoint Overpass API URL (default the public instance — override for your
#'   own).
#' @param contact optional contact string appended to the User-Agent (opt-in; e.g. a
#'   project URL or email you choose to share). \code{NULL} sends none.
#' @param tiles split the bbox into a \code{tiles} x \code{tiles} grid of requests
#'   (default \code{1}).
#' @param pause minimum seconds between requests (default \code{1}); also the base of
#'   the retry backoff.
#' @param max_tries attempts per tile before giving up (default \code{4}).
#' @param timeout Overpass server-side \code{[timeout:]} seconds (default \code{50}).
#' @return an \pkg{sf} object of \code{LINESTRING}s with \code{osm_id}, \code{highway}
#'   and \code{class} columns, in \code{crs}; attributes \code{timestamp} (OSM base
#'   time) and \code{path} (where it was written). Pass straight to
#'   \code{\link{cost_distance}} as \code{sources}.
#' @export
osm_access <- function(bbox, crs = NULL,
                       roads = c("motorway", "trunk", "primary", "secondary",
                                 "tertiary", "unclassified", "residential",
                                 "service", "living_street", "road"),
                       trails = c("path", "footway", "track", "bridleway",
                                  "cycleway", "steps"),
                       path = NULL, cache_dir = NULL, overwrite = FALSE,
                       endpoint = "https://overpass-api.de/api/interpreter",
                       contact = NULL, tiles = 1L, pause = 1, max_tries = 4L,
                       timeout = 50) {
  stopifnot(requireNamespace("httr2", quietly = TRUE),
            requireNamespace("jsonlite", quietly = TRUE))
  hw <- c(roads, trails)
  if (!length(hw)) stop("nothing to fetch: both 'roads' and 'trails' are empty")

  bb <- .lonlat_bbox(bbox)                 # list(box = c(xmin,ymin,xmax,ymax), crs)
  out_crs <- if (!is.null(crs)) sf::st_crs(crs)
             else if (!is.na(bb$crs)) bb$crs else sf::st_crs(4326)

  # Output + checkpoint locations.
  if (is.null(path)) path <- file.path(getwd(), "osm_access.gpkg")
  if (is.null(cache_dir))
    cache_dir <- file.path(dirname(path),
                           paste0(tools::file_path_sans_ext(basename(path)), "_tiles"))
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

  ua <- paste0("lsed-osm/0.1 (R; polite)",
               if (!is.null(contact)) paste0(" ", contact))
  rx <- paste(hw, collapse = "|")

  # One throttled, backing-off request for a single lon/lat box -> parsed JSON.
  fetch_box <- function(b) {
    ql <- sprintf(
      "[out:json][timeout:%d];(way[\"highway\"~\"^(%s)$\"](%.6f,%.6f,%.6f,%.6f););out geom;",
      timeout, rx, b[2], b[1], b[4], b[3])        # Overpass order: S,W,N,E
    resp <- httr2::request(endpoint) |>
      httr2::req_user_agent(ua) |>
      httr2::req_body_form(data = ql) |>
      httr2::req_timeout(timeout + 15) |>
      httr2::req_throttle(capacity = 1, fill_time_s = pause) |>
      httr2::req_retry(max_tries = max_tries,
                       backoff = function(i) min(60, pause * 2^i),
                       is_transient = function(r)
                         httr2::resp_status(r) %in% c(429, 502, 503, 504)) |>
      httr2::req_perform()
    jsonlite::fromJSON(httr2::resp_body_string(resp), simplifyVector = FALSE)
  }

  # Parsed Overpass JSON -> sf (0+ rows, lon/lat), carrying the OSM timestamp.
  parse_sf <- function(j) {
    ways <- Filter(function(e) e$type == "way" && !is.null(e$geometry), j$elements)
    ts <- j$osm3s$timestamp_osm_base
    if (!length(ways)) return(structure(.empty_access_sf(), timestamp = ts))
    ids <- vapply(ways, function(w) as.character(w$id), "")
    geoms <- sf::st_sfc(lapply(ways, function(w)
      sf::st_linestring(do.call(rbind, lapply(w$geometry,
                                              function(p) c(p$lon, p$lat))))),
      crs = 4326)
    htype <- vapply(ways, function(w)
      if (!is.null(w$tags$highway)) w$tags$highway else NA_character_, "")
    cls <- ifelse(htype %in% roads, "road", ifelse(htype %in% trails, "trail", NA))
    structure(sf::st_sf(osm_id = ids, highway = htype, class = cls, geometry = geoms),
              timestamp = ts)
  }

  # Batch into tiles; checkpoint each, skip ones already on disk.
  xs <- seq(bb$box[1], bb$box[3], length.out = tiles + 1L)
  ys <- seq(bb$box[2], bb$box[4], length.out = tiles + 1L)
  grid <- expand.grid(i = seq_len(tiles), j = seq_len(tiles))
  tiles_sf <- vector("list", nrow(grid))
  for (g in seq_len(nrow(grid))) {
    i <- grid$i[g]; j <- grid$j[g]
    f <- file.path(cache_dir, sprintf("tile_%03d_%03d.rds", i, j))
    if (file.exists(f) && !overwrite) {
      tiles_sf[[g]] <- readRDS(f)
      next
    }
    sfj <- parse_sf(fetch_box(c(xs[i], ys[j], xs[i + 1L], ys[j + 1L])))
    saveRDS(sfj, f)                                 # persist as it completes
    tiles_sf[[g]] <- sfj
  }

  stamp <- attr(tiles_sf[[1]], "timestamp")
  allw <- do.call(rbind, tiles_sf)
  allw <- allw[!duplicated(allw$osm_id), ]          # tile seams share boundary ways
  if (nrow(allw) && out_crs != sf::st_crs(4326)) allw <- sf::st_transform(allw, out_crs)
  if (!nrow(allw)) warning("no highways returned for this bbox")

  sf::st_write(allw, path, quiet = TRUE, delete_dsn = TRUE)
  message(sprintf("osm_access: %d ways -> %s", nrow(allw), path))
  structure(allw, timestamp = stamp, path = path)
}

# ---------------------------------------------------------------------------
# internal
# ---------------------------------------------------------------------------

# An empty access-network sf with the canonical columns (lon/lat), so empty tiles
# rbind cleanly and cache round-trips preserve schema.
.empty_access_sf <- function() {
  sf::st_sf(osm_id = character(0), highway = character(0), class = character(0),
            geometry = sf::st_sfc(crs = 4326))
}

# Reduce an sf/terra object or a numeric vector to a lon/lat bbox + its source CRS.
.lonlat_bbox <- function(bbox) {
  if (is.numeric(bbox)) {
    stopifnot(length(bbox) == 4)
    return(list(box = as.numeric(bbox), crs = sf::st_crs(NA)))
  }
  src <- tryCatch(sf::st_crs(bbox), error = function(e) sf::st_crs(NA))
  e <- if (inherits(bbox, c("SpatRaster", "SpatVector", "SpatExtent"))) {
    ex <- terra::ext(bbox); c(ex[1], ex[3], ex[2], ex[4])   # xmin,ymin,xmax,ymax
  } else {
    as.numeric(sf::st_bbox(bbox))[c("xmin", "ymin", "xmax", "ymax")]
  }
  box <- sf::st_bbox(c(xmin = e[1], ymin = e[2], xmax = e[3], ymax = e[4]), crs = src)
  if (!is.na(src) && src != sf::st_crs(4326))
    box <- sf::st_bbox(sf::st_transform(sf::st_as_sfc(box), 4326))
  list(box = as.numeric(box)[c("xmin", "ymin", "xmax", "ymax")], crs = src)
}
