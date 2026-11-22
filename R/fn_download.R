#' Acquire every external input for the patch ranking (download module)
#'
#' The data-acquisition half of the patch ranking, kept separate from the calculation
#' so the slow, network-bound, API-dependent work runs once and its outputs can be
#' re-used. It overlaps the two big pulls: GBIF's bulk download is \strong{requested}
#' first (and prepared server-side while we wait), the OSM road/trail network is
#' fetched in that window, and only then do we come \strong{back} to collect the
#' finished GBIF archive — so the GBIF prep time is spent fetching OSM, not idling.
#'
#' \preformatted{
#'   gbif_request()  --->  (GBIF prepares server-side)  ---.
#'                                                          |
#'   osm_network_for()  (runs during the wait)             |
#'                                                          v
#'                                                   gbif_retrieve()
#' }
#'
#' Needs a GBIF account (see \code{\link{gbif_request}}). Everything is written to disk
#' (GeoPackages), so a later \code{rank} calculation reads files, not live APIs.
#'
#' @param patches \pkg{sf} patches (points/polygons) that scope both pulls.
#' @param origin trip origin for accessibility — \pkg{sf} point or \code{c(lon, lat)}
#'   (see \code{\link{osm_network_for}}).
#' @param crs planar working CRS (metres); defaults to the patches' CRS.
#' @param taxon_key,basis,year_from,coord_uncertainty_max GBIF filters
#'   (\code{\link{gbif_request}}).
#' @param gbif_path,osm_path output GeoPackages for the records and the network.
#' @param eps_km,buffer_m cluster-scoping for the OSM pull
#'   (\code{\link{osm_network_for}}).
#' @param ... further arguments to \code{\link{osm_network_for}} (e.g.
#'   \code{connect_roads}, \code{endpoint}, \code{contact}, \code{pause}).
#' @return a list: \code{records} (GBIF \pkg{sf}), \code{network} (OSM \pkg{sf}),
#'   \code{origin} (\pkg{sf} point), \code{clusters}, and \code{doi} — everything the
#'   calculation step needs.
#' @seealso \code{\link{gbif_request}}, \code{\link{osm_network_for}},
#'   \code{\link{rank_index}}
#' @export
rank_download <- function(patches, origin, crs = NULL,
                          taxon_key = 7707728,
                          basis = c("PRESERVED_SPECIMEN", "HUMAN_OBSERVATION"),
                          year_from = 1990, coord_uncertainty_max = NA,
                          gbif_path = "data/gbif_records.gpkg",
                          osm_path = "data/osm_network.gpkg",
                          eps_km = 5, buffer_m = 1500, ...) {
  pcrs <- if (!is.null(crs)) sf::st_crs(crs) else sf::st_crs(patches)

  # 1. Kick off the GBIF job — returns at once, GBIF prepares it server-side.
  key <- gbif_request(patches, taxon_key = taxon_key, basis = basis,
                      year_from = year_from,
                      coord_uncertainty_max = coord_uncertainty_max)

  # 2. Fetch OSM while GBIF prepares (this is the overlap that saves the wait).
  message("Fetching OSM network while GBIF prepares the download...")
  net <- osm_network_for(patches, origin, crs = pcrs, eps_km = eps_km,
                         buffer_m = buffer_m, path = osm_path, ...)

  # 3. Come back for the finished GBIF archive (usually ready by now).
  message("Collecting the GBIF download...")
  recs <- gbif_retrieve(key, path = gbif_path, crs = pcrs)

  list(records = recs, network = net$network, origin = net$origin,
       clusters = net$clusters, doi = attr(recs, "doi"))
}
