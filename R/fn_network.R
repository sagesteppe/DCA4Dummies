#' Accessibility (A) and complementarity (C) on a road/trail network graph
#'
#' A vector-first engine for the operational terms. Travel runs \strong{on the
#' linestrings} — roads and trails as a weighted \pkg{igraph}, each edge timed by its
#' class speed — and only the \strong{last hop} from the network to a patch centroid
#' drops to an off-trail Tobler cost. This sidesteps the resolution trap of a friction
#' raster (at 800 m a road is one fat cell) and, because the patches become graph
#' nodes, the same network answers \strong{C} (\code{\link{complementarity}}) for free.
#'
#' \strong{Model.} Roads are the drive-to origins (reachable at zero walking cost);
#' from there a crew walks trails (\code{trail_speed}) and finally cuts cross-country
#' to the centroid at a slope-dependent Tobler speed. Accessibility of a patch is the
#' least-cost path origin → network → hop, in hours.
#'
#' \strong{Noding.} OSM ways share vertex coordinates at junctions, so the graph is
#' built by snapping segment endpoints to a \code{snap}-metre grid — correct for OSM
#' topology. For arbitrary linework that crosses without shared vertices, pre-node it
#' (e.g. \code{sf::st_node}) before calling.
#'
#' @name network
NULL

# ---------------------------------------------------------------------------
# Network construction
# ---------------------------------------------------------------------------

#' Build a travel-time graph from road/trail linestrings
#'
#' @param lines \pkg{sf} LINESTRINGs in a projected CRS (metres), with a \code{class}
#'   column ("road"/"trail"); \code{\link{osm_access}} output fits directly.
#' @param road_speed,trail_speed class speeds in km/h (defaults 30, 4). Roads mainly
#'   serve as origins, so their speed rarely matters.
#' @param snap tight grid (metres) that segment endpoints snap to when forming nodes
#'   (default 1) — within-class noding, so parallel switchbacks/lanes are not fused.
#' @param connect liberal tolerance (metres, default 50) for bridging a \strong{stub
#'   end} of one class to the nearest node of the \emph{other} class — i.e. closing
#'   trailhead gaps where a trail dead-ends near a road (or a road near a trail)
#'   without a shared OSM node. Set \code{0} to disable. Only degree-1 ends are
#'   bridged, so it never fuses same-class features.
#' @return a list: \code{graph} (an undirected \pkg{igraph} with edge \code{weight} in
#'   hours and \code{class} — road/trail/connect, vertex coords \code{x}/\code{y}) and
#'   \code{nodes} (data.frame of node \code{name}, \code{x}, \code{y}).
#' @keywords internal
#' @export
build_network <- function(lines, road_speed = 30, trail_speed = 4, snap = 1,
                          connect = 50) {
  stopifnot(requireNamespace("igraph", quietly = TRUE))
  cls <- if (!is.null(lines$class)) lines$class else "trail"
  spd <- ifelse(cls == "road", road_speed, trail_speed)   # km/h, per line

  geo <- sf::st_geometry(lines)
  key <- function(m) sprintf("%.0f_%.0f", round(m[, 1] / snap), round(m[, 2] / snap))

  parts <- lapply(seq_along(geo), function(i) {
    m <- sf::st_coordinates(geo[[i]])[, 1:2, drop = FALSE]
    n <- nrow(m)
    if (n < 2) return(NULL)
    a <- seq_len(n - 1); b <- a + 1L
    len <- sqrt((m[b, 1] - m[a, 1])^2 + (m[b, 2] - m[a, 2])^2)   # metres
    data.frame(from = key(m[a, , drop = FALSE]), to = key(m[b, , drop = FALSE]),
               weight = (len / 1000) / spd[i], class = cls[i],   # hours
               fx = m[a, 1], fy = m[a, 2], tx = m[b, 1], ty = m[b, 2],
               stringsAsFactors = FALSE)
  })
  ed <- do.call(rbind, parts)
  ed <- ed[ed$from != ed$to & ed$weight > 0, , drop = FALSE]     # drop zero-length

  nodes <- rbind(data.frame(name = ed$from, x = ed$fx, y = ed$fy),
                 data.frame(name = ed$to,   x = ed$tx, y = ed$ty))
  nodes <- nodes[!duplicated(nodes$name), , drop = FALSE]
  rownames(nodes) <- nodes$name

  ge <- ed[, c("from", "to", "weight", "class")]

  # --- class-aware bridging of stub ends (trailhead gaps) ---
  if (connect > 0) {
    deg <- table(c(ed$from, ed$to))
    stubs <- names(deg)[deg == 1]
    road_nodes  <- unique(c(ed$from[ed$class == "road"],  ed$to[ed$class == "road"]))
    trail_nodes <- unique(c(ed$from[ed$class == "trail"], ed$to[ed$class == "trail"]))

    bridge <- function(src_cand, dst_names) {
      src <- intersect(stubs, src_cand)
      if (!length(src) || !length(dst_names)) return(NULL)
      s_sf <- sf::st_as_sf(nodes[src, ], coords = c("x", "y"), crs = sf::st_crs(lines))
      d_sf <- sf::st_as_sf(nodes[dst_names, ], coords = c("x", "y"), crs = sf::st_crs(lines))
      idx <- sf::st_nearest_feature(s_sf, d_sf)
      gap <- as.numeric(sf::st_distance(s_sf, d_sf[idx, ], by_element = TRUE))
      j <- dst_names[idx]
      keep <- gap <= connect & src != j
      if (!any(keep)) return(NULL)
      data.frame(from = src[keep], to = j[keep],
                 weight = (gap[keep] / 1000) / trail_speed, class = "connect",
                 stringsAsFactors = FALSE)
    }
    br <- rbind(bridge(trail_nodes, road_nodes), bridge(road_nodes, trail_nodes))
    if (!is.null(br) && nrow(br)) {
      dk <- paste(pmin(br$from, br$to), pmax(br$from, br$to), sep = "|")
      ge <- rbind(ge, br[!duplicated(dk), , drop = FALSE])
    }
  }

  g <- igraph::graph_from_data_frame(ge, directed = FALSE, vertices = nodes)
  list(graph = g, nodes = nodes)
}

# ---------------------------------------------------------------------------
# Accessibility (A)
# ---------------------------------------------------------------------------

#' Accessibility of patches over the network (A)
#'
#' Finds the travel time from the \strong{origin} (field office / trailhead) to each
#' patch: drive the roads, walk the trails, then an off-trail Tobler hop to the
#' centroid. Each patch is joined to its \code{k} nearest \strong{origin-reachable}
#' nodes and takes the cheapest total route, so a patch beside an isolated trail
#' fragment bushwhacks to the usable network (a finite cost) instead of coming back
#' \code{NA}, and can enter via a connected trail it passes rather than a dead-end
#' stub. Returns the patches with a rescaled score \emph{and} the augmented graph (with
#' patches as nodes), so \code{\link{complementarity}} runs on the same structure.
#'
#' @param lines road/trail \pkg{sf} lines (see \code{\link{build_network}}).
#' @param patches \pkg{sf} points or polygons; polygons are reduced to centroids.
#' @param dem optional elevation \code{SpatRaster} (same CRS) for the Tobler hop slope;
#'   without it the hop uses \code{offtrail_speed} on the flat.
#' @param origins trip start as \pkg{sf} point(s) — typically the \strong{field
#'   office}, snapped free to the nearest road, after which road driving accrues time
#'   (often the dominant cost). \code{NULL} (default) treats every road node as a
#'   zero-cost origin ("distance from the road network", ignoring the drive).
#' @param k nearest origin-reachable nodes each patch may enter through (default
#'   \code{3}); the route uses whichever minimises drive/walk + hop.
#' @param road_speed,trail_speed,snap,connect passed to \code{\link{build_network}}
#'   (\code{connect} bridges trailhead gaps between roads and trails; default 50 m).
#' @param offtrail_speed unpenalised flat-ground off-trail speed in km/h when
#'   \code{dem} is absent (default \code{5}, ~Tobler on the level).
#' @param offtrail_penalty multiplier on the off-trail hop speed (default \code{0.5}):
#'   cross-country through brush/talus is slower than the slope-only Tobler ideal, so
#'   this keeps trails faster than bushwhacking. \code{1} removes the penalty.
#' @param invert if \code{TRUE} (default), the returned \code{A} is
#'   \eqn{1 - \tilde t} (accessible = high), the index form; raw hours are kept in
#'   \code{access_hr}.
#' @return a list: \code{patches} (input patches + \code{access_hr} and \code{A}),
#'   \code{net} (the augmented \code{build_network} list), and \code{patch_nodes} (the
#'   patch vertex names) for reuse by \code{\link{complementarity}}.
#' @export
accessibility <- function(lines, patches, dem = NULL, origins = NULL, k = 3,
                          road_speed = 30, trail_speed = 4, snap = 1, connect = 50,
                          offtrail_speed = 5, offtrail_penalty = 0.5, invert = TRUE) {
  stopifnot(requireNamespace("igraph", quietly = TRUE),
            requireNamespace("nabor", quietly = TRUE))
  net <- build_network(lines, road_speed, trail_speed, snap, connect)
  g <- net$graph
  nodes <- net$nodes
  node_pts <- sf::st_as_sf(nodes, coords = c("x", "y"), crs = sf::st_crs(lines))

  # origins: field-office points snapped to nearest ROAD node (free), else all roads.
  road_names <- unique(as.vector(igraph::ends(
    g, igraph::E(g)[igraph::E(g)$class == "road"], names = TRUE)))
  if (!length(road_names)) road_names <- nodes$name          # no roads -> any node
  if (is.null(origins)) {
    onames <- road_names
  } else {
    origins <- sf::st_transform(origins, sf::st_crs(lines))
    road_pts <- sf::st_as_sf(nodes[road_names, ], coords = c("x", "y"),
                             crs = sf::st_crs(lines))
    onames <- unique(road_names[sf::st_nearest_feature(origins, road_pts)])
  }

  # Cost origin -> every network node once; the reachable set is the usable network.
  d0 <- apply(igraph::distances(g, v = onames, weights = igraph::E(g)$weight), 2, min)
  reach <- names(d0)[is.finite(d0)]
  if (!length(reach)) stop("no network nodes are reachable from the origin(s)")
  reach_xy <- as.matrix(nodes[reach, c("x", "y")])

  # patch centroids
  patches <- sf::st_transform(patches, sf::st_crs(lines))
  cen <- if (all(sf::st_geometry_type(patches) %in% c("POINT", "MULTIPOINT")))
           patches else sf::st_centroid(sf::st_geometry(patches))
  pc <- sf::st_coordinates(cen)[, 1:2, drop = FALSE]

  # k nearest reachable entry nodes per patch; total cost = drive/walk-in + Tobler hop.
  kk <- min(k, length(reach))
  knn <- nabor::knn(reach_xy, pc, k = kk)
  cost <- matrix(Inf, nrow(pc), kk)
  hops <- matrix(0,   nrow(pc), kk)
  tgt  <- matrix("",  nrow(pc), kk)
  for (c in seq_len(kk)) {
    idx <- knn$nn.idx[, c]; run <- knn$nn.dists[, c]
    nm  <- reach[idx]
    spd <- .hop_speed(pc, reach_xy[idx, , drop = FALSE], run, dem,
                      offtrail_speed, offtrail_penalty)
    hops[, c] <- (run / 1000) / spd
    tgt[, c]  <- nm
    cost[, c] <- d0[nm] + hops[, c]
  }
  acc <- apply(cost, 1, min)                 # hours via the cheapest entry

  # Add patches as nodes with their k connector edges (for complementarity reuse).
  pnames <- paste0("patch_", seq_len(nrow(pc)))
  g <- igraph::add_vertices(g, nrow(pc), name = pnames, x = pc[, 1], y = pc[, 2])
  from <- rep(pnames, times = kk)            # column-major, aligns with as.vector(tgt)
  g <- igraph::add_edges(g, as.vector(rbind(from, as.vector(tgt))),
                         weight = as.vector(hops), class = "offtrail")

  patches$access_hr <- acc
  patches$A <- .rescale01_vec(acc, invert)
  list(patches = patches,
       net = list(graph = g,
                  nodes = rbind(nodes, data.frame(name = pnames,
                                x = pc[, 1], y = pc[, 2]))),
       patch_nodes = pnames)
}

# ---------------------------------------------------------------------------
# internal
# ---------------------------------------------------------------------------

# Off-trail Tobler speed (km/h) for each hop, times the off-trail penalty. With a DEM,
# slope = |dz| / run from the centroid to its network node; without, a flat
# offtrail_speed. The penalty makes bushwhacking slower than the slope-only ideal.
.hop_speed <- function(pc, node_xy, run, dem, offtrail_speed, penalty = 0.5) {
  if (is.null(dem) || !inherits(dem, "SpatRaster"))
    return(rep(offtrail_speed * penalty, nrow(pc)))
  z1 <- terra::extract(dem, pc)[, 1]
  z2 <- terra::extract(dem, as.matrix(node_xy))[, 1]
  S <- ifelse(run > 0, abs(z2 - z1) / run, 0)
  spd <- 6 * exp(-3.5 * abs(S + 0.05))           # Tobler km/h
  spd[!is.finite(spd) | spd <= 0] <- offtrail_speed
  spd * penalty
}

# Min-max rescale a numeric vector to [0, 1], NA-safe, optionally complemented.
.rescale01_vec <- function(x, invert = FALSE) {
  rng <- range(x, na.rm = TRUE)
  out <- if (is.finite(rng[2]) && rng[2] > rng[1]) (x - rng[1]) / (rng[2] - rng[1])
         else x * 0
  if (invert) out <- 1 - out
  out
}
