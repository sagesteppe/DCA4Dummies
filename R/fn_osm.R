#' Cluster patches into spatial groups
#'
#' Groups patch centroids that lie within \code{eps_km} of one another (a
#' within-distance graph, connected components), so OSM can be fetched only where
#' patches actually are rather than across the empty country between them.
#'
#' @param patches \pkg{sf} points or polygons.
#' @param eps_km join distance in km: patches within this of each other share a cluster
#'   (default \code{5}).
#' @param crs planar CRS (metres) for the distance test; defaults to the patches' CRS.
#' @return the patch centroids as \pkg{sf} points with a \code{cluster} id column.
#' @export
cluster_patches <- function(patches, eps_km = 5, crs = NULL) {
  stopifnot(requireNamespace("igraph", quietly = TRUE))
  if (!is.null(crs)) patches <- sf::st_transform(patches, sf::st_crs(crs))
  cen <- sf::st_as_sf(sf::st_centroid(sf::st_geometry(patches)))
  nb <- sf::st_is_within_distance(cen, dist = eps_km * 1000)
  g <- igraph::graph_from_adj_list(nb, mode = "all")
  cen$cluster <- igraph::components(g)$membership
  cen
}

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
#' @param exclude optional geometry (\pkg{sf}/\pkg{sfc}, any CRS) of ground already
#'   covered elsewhere; any tile falling \emph{entirely} within it is skipped — never
#'   fetched — so overlapping pulls don't re-query the API. Its benefit grows as tiles
#'   shrink toward the excluded regions' size.
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
                       contact = NULL, tiles = 1L, exclude = NULL, pause = 1,
                       max_tries = 4L, timeout = 50) {
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

  # Any ground already covered elsewhere — tiles fully inside it are never fetched.
  excl <- if (!is.null(exclude))
    sf::st_union(sf::st_transform(sf::st_geometry(sf::st_as_sf(exclude)), 4326))

  # Batch into tiles; checkpoint each, skip ones already on disk or excluded.
  xs <- seq(bb$box[1], bb$box[3], length.out = tiles + 1L)
  ys <- seq(bb$box[2], bb$box[4], length.out = tiles + 1L)
  grid <- expand.grid(i = seq_len(tiles), j = seq_len(tiles))
  tiles_sf <- vector("list", nrow(grid))
  for (g in seq_len(nrow(grid))) {
    i <- grid$i[g]; j <- grid$j[g]
    tb <- c(xs[i], ys[j], xs[i + 1L], ys[j + 1L])
    # Skip a tile fully within already-covered ground — don't hit the API for it.
    if (!is.null(excl)) {
      tpoly <- sf::st_as_sfc(sf::st_bbox(c(xmin = tb[1], ymin = tb[2],
                                           xmax = tb[3], ymax = tb[4]), crs = 4326))
      if (as.logical(sf::st_covered_by(tpoly, excl, sparse = FALSE))) {
        tiles_sf[[g]] <- .empty_access_sf(); next
      }
    }
    f <- file.path(cache_dir, sprintf("tile_%03d_%03d.rds", i, j))
    if (file.exists(f) && !overwrite) {
      tiles_sf[[g]] <- readRDS(f)
      next
    }
    sfj <- parse_sf(fetch_box(tb))
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
    as.numeric(sf::st_bbox(bbox)[c("xmin", "ymin", "xmax", "ymax")])  # name, then strip
  }
  box <- sf::st_bbox(c(xmin = e[1], ymin = e[2], xmax = e[3], ymax = e[4]), crs = src)
  if (!is.na(src) && src != sf::st_crs(4326))
    box <- sf::st_bbox(sf::st_transform(sf::st_as_sfc(box), 4326))
  list(box = as.numeric(box[c("xmin", "ymin", "xmax", "ymax")]), crs = src)
}

#' Fetch an OSM network scoped to patch clusters, connected to an origin
#'
#' Pulls dense roads+trails only within boxes around each patch cluster
#' (\code{\link{cluster_patches}}), plus a \strong{roads-only} connective layer over
#' the envelope of the origin and all clusters. The field office then drives in on one
#' connected network to every cluster, without fetching the empty country between them
#' (the hard part of a scattered domain). Reuses \code{\link{osm_access}} for every
#' sub-pull, so caching, throttling and back-off all apply.
#'
#' @param patches \pkg{sf} points/polygons to scope the pull to.
#' @param origin the trip origin (field office): an \pkg{sf} point, or a numeric
#'   \code{c(lon, lat)} in EPSG:4326.
#' @param crs planar CRS (metres) to work in; defaults to the patches' CRS.
#' @param eps_km cluster join distance (km, default \code{5};
#'   see \code{\link{cluster_patches}}).
#' @param buffer_m buffer (m) added around each cluster's bounding box (default
#'   \code{1500}).
#' @param roads,trails OSM \code{highway} classes pulled densely \emph{within} cluster
#'   boxes (\code{trails} only there).
#' @param connect_roads the (major) road classes for the \emph{connective} envelope
#'   layer that links the origin to the clusters — default the through-routes
#'   \code{motorway/trunk/primary/secondary/tertiary}. Residential/service/track are
#'   left out here: they are local bloat and already come from the cluster pulls.
#' @param path GeoPackage for the merged network (\code{NULL} -> cwd).
#' @param cache_dir folder for the per-pull checkpoints (\code{NULL} -> beside
#'   \code{path}); existing pulls are reused, so the whole fetch is resumable.
#' @param tiles_connect tiles for the connective roads pull (\code{NULL} -> ~1 per
#'   30 km of envelope width).
#' @param ... passed to \code{\link{osm_access}} (e.g. \code{endpoint}, \code{contact},
#'   \code{pause}).
#' @return a list: \code{network} (merged \pkg{sf} lines, deduped), \code{origin}
#'   (\pkg{sf} point in \code{crs}) and \code{clusters} (centroids with \code{cluster}
#'   ids); attribute \code{path} is where the network was written. Feed
#'   \code{network} + \code{origin} to \code{\link{accessibility}}.
#' @export
osm_network_for <- function(patches, origin, crs = NULL, eps_km = 5, buffer_m = 1500,
                            roads = c("motorway", "trunk", "primary", "secondary",
                                      "tertiary", "unclassified", "residential",
                                      "service", "living_street", "road", "track"),
                            trails = c("path", "footway", "bridleway", "cycleway",
                                       "steps"),
                            connect_roads = c("motorway", "trunk", "primary",
                                              "secondary", "tertiary"),
                            path = NULL, cache_dir = NULL, tiles_connect = NULL, ...) {
  pcrs <- if (!is.null(crs)) sf::st_crs(crs) else sf::st_crs(patches)
  patches <- sf::st_transform(patches, pcrs)
  org <- if (inherits(origin, c("sf", "sfc", "SpatVector")))
           sf::st_transform(sf::st_as_sf(origin), pcrs)
         else sf::st_transform(
           sf::st_as_sf(sf::st_sfc(sf::st_point(c(origin[1], origin[2])), crs = 4326)),
           pcrs)

  cc <- cluster_patches(patches, eps_km = eps_km, crs = pcrs)
  boxes <- lapply(split(seq_len(nrow(cc)), cc$cluster), function(ix)
    sf::st_buffer(sf::st_as_sfc(sf::st_bbox(cc[ix, ])), buffer_m))

  if (is.null(path)) path <- file.path(getwd(), "osm_network_for.gpkg")
  if (is.null(cache_dir))
    cache_dir <- file.path(dirname(path),
                           paste0(tools::file_path_sans_ext(basename(path)), "_cache"))
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

  # local dense pulls (roads + trails) per cluster box
  locals <- lapply(seq_along(boxes), function(i)
    osm_access(boxes[[i]], crs = pcrs, roads = roads, trails = trails,
               path = file.path(cache_dir, sprintf("cluster_%03d.gpkg", i)),
               cache_dir = file.path(cache_dir, sprintf("cluster_%03d_tiles", i)), ...))

  # connective roads-only pull over the {origin + clusters} envelope
  env <- sf::st_as_sfc(sf::st_bbox(do.call(c, c(list(sf::st_geometry(org)),
                                               lapply(boxes, sf::st_geometry)))))
  if (is.null(tiles_connect)) {
    w <- as.numeric(diff(sf::st_bbox(env)[c("xmin", "xmax")]))
    tiles_connect <- max(1L, ceiling(w / 30000))
  }
  conn <- osm_access(env, crs = pcrs, roads = connect_roads, trails = character(0),
                     tiles = tiles_connect,
                     exclude = do.call(c, lapply(boxes, sf::st_geometry)),
                     path = file.path(cache_dir, "connect.gpkg"),
                     cache_dir = file.path(cache_dir, "connect_tiles"), ...)

  net <- do.call(rbind, c(locals, list(conn)))
  net <- net[!duplicated(net$osm_id), , drop = FALSE]
  sf::st_write(net, path, quiet = TRUE, delete_dsn = TRUE)
  message(sprintf("osm_network_for: %d ways, %d clusters -> %s",
                  nrow(net), length(boxes), path))
  structure(list(network = net, origin = org, clusters = cc), path = path)
}
