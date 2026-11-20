#' Patch and subpatch delineation for survey-unit mapping
#'
#' The delineation pipeline behind the manuscript's survey units, factored out of
#' \code{Analyses.qmd} so the stages are reusable and testable. The flow is
#' \strong{watershed -> patches -> snap/merge -> contour (band) -> k-means cap}:
#'
#' \enumerate{
#'   \item \code{\link{delineate_patches}} binarises the probability surface at a
#'     suitability threshold and runs a marker-controlled watershed (via
#'     WhiteboxTools) to break the suitable network at low necks, with a rook's-case
#'     component fallback for cells that drain off-mask, then drops sub-\code{min_area}
#'     slivers. This is the shared front end for both subdivision levers.
#'   \item \code{\link{band_subpatches}} subdivides each oversized patch: it
#'     \emph{contours} the raw surface into fixed suitability bands, finds same-band
#'     connected pieces, \emph{snaps} sub-threshold pieces into their highest-band
#'     contiguous neighbour (\code{\link{snap_pieces}}), then geometrically caps any
#'     surviving unit still above the maximum size (\code{\link{kcap}}).
#' }
#'
#' Base R + \pkg{terra} + \pkg{whitebox} only -- no \pkg{igraph}; the region-adjacency
#' and connected-component work is done with the small union-find helpers below
#' (\code{\link{rook_edges}}, \code{\link{uf_comp}}, \code{\link{snap_pieces}}).
#' These surfaces are \strong{raw / not Platt-corrected}; the taxon is a ridge/talus
#' endemic, so linear high-suitability filaments are real habitat and are never
#' smoothed.
#'
#' @name patches
NULL

# ---------------------------------------------------------------------------
# Stage 1: watershed -> patches
# ---------------------------------------------------------------------------

#' Delineate patches by marker-controlled watershed
#'
#' Binarise \code{pr} at \code{sens_spec}, then split the suitable network with a
#' marker-controlled watershed on the \emph{inverted} surface (high-probability cores
#' become pits; necks shallower than \code{delta} are merged by depression filling,
#' surviving pits seed the watershed). Suitable cells that drain off-mask (small
#' patches with no interior pit) keep an identity through a rook's-case (4-connected)
#' component fallback, so none are dropped. Patches at or below \code{min_area} m^2
#' are removed as un-surveyable slivers. Watershed intermediates are written to
#' \code{rdir} so a partial re-run can reuse finished steps.
#'
#' @param pr A \pkg{terra} \code{SpatRaster} of predicted probability (one layer).
#' @param sens_spec Suitable-habitat threshold in probability units; cells
#'   \code{>= sens_spec} are suitable.
#' @param delta Prominence (minimum neck depth, probability units) below which two
#'   cores are merged rather than split. Passed to \code{wbt_fill_depressions(max_depth=)}.
#' @param min_area Patches with total area at or below this (m^2) are dropped.
#' @param rdir Directory for the watershed intermediate rasters (dem/fill/sink/pntr/seg).
#'
#' @return A \code{SpatRaster} of integer patch ids (\code{NA} off-patch), aligned to \code{pr}.
#' @seealso \code{\link{band_subpatches}}
#' @export
delineate_patches <- function(pr, sens_spec, delta, min_area, rdir) {
  mask <- pr >= sens_spec

  # Hydrological inversion trick: cores -> pits, fill necks shallower than delta,
  # label surviving pits as seeds, route every suitable cell to the core it drains to.
  dem <- terra::ifel(mask, 1 - pr, NA)
  terra::writeRaster(dem, file.path(rdir, "patch-dem.tif"), overwrite = TRUE)
  whitebox::wbt_fill_depressions(
    file.path(rdir, "patch-dem.tif"),
    file.path(rdir, "patch-fill.tif"),
    max_depth = delta
  )
  whitebox::wbt_sink(
    file.path(rdir, "patch-fill.tif"),
    file.path(rdir, "patch-sink.tif"),
    zero_background = TRUE
  )
  whitebox::wbt_d8_pointer(
    file.path(rdir, "patch-fill.tif"),
    file.path(rdir, "patch-pntr.tif")
  )
  whitebox::wbt_watershed(
    file.path(rdir, "patch-pntr.tif"),
    file.path(rdir, "patch-sink.tif"),
    file.path(rdir, "patch-seg.tif")
  )

  seg <- terra::rast(file.path(rdir, "patch-seg.tif"))
  seg[seg <= 0] <- NA

  # Fallback for suitable cells that drain off-mask: rook's-case components, id-shifted
  # above the watershed labels so the two id spaces never collide.
  unlab <- mask & is.na(seg)
  cc <- terra::patches(
    terra::ifel(unlab, 1, NA),
    directions = 4,
    zeroAsNA = TRUE
  ) +
    terra::global(seg, "max", na.rm = TRUE)[[1]]
  patch_id <- terra::mask(terra::cover(seg, cc), mask, maskvalues = c(NA, 0))

  # Drop sub-threshold slivers (too small to be a meaningful survey unit).
  cell_area <- prod(terra::res(patch_id))
  fr <- terra::freq(patch_id)
  keep_ids <- fr$value[fr$count * cell_area > min_area]
  terra::ifel(patch_id %in% keep_ids, patch_id, NA)
}

# ---------------------------------------------------------------------------
# Low-level helpers: adjacency, components, snap/merge, geometric cap
# ---------------------------------------------------------------------------

#' Rook's-case adjacency among a set of raster cells
#'
#' East/south neighbour edges among a set of global cell indices, returned as local
#' index pairs (positions within \code{gc}). East+south only, since each undirected
#' rook edge is counted once from its lower-index endpoint.
#'
#' @param gc Integer vector of global (1-based) cell indices.
#' @param nc,nr Number of columns and rows in the parent raster.
#' @return A list with integer vectors \code{i} and \code{j}: local index pairs of adjacent cells.
#' @keywords internal
#' @export
rook_edges <- function(gc, nc, nr) {
  r <- (gc - 1L) %/% nc
  cc <- (gc - 1L) %% nc
  eE <- gc[cc < nc - 1L]
  jE <- match(eE + 1L, gc)
  okE <- !is.na(jE)
  eS <- gc[r < nr - 1L]
  jS <- match(eS + nc, gc)
  okS <- !is.na(jS)
  list(i = c(match(eE[okE], gc), match(eS[okS], gc)), j = c(jE[okE], jS[okS]))
}

#' Connected components by union-find
#'
#' Label \code{n} nodes 1..m by connected component given an undirected edge subset,
#' with path-compressed union-find. Labels are assigned in order of first appearance.
#'
#' @param n Number of nodes.
#' @param ei,ej Integer vectors of equal length: the endpoints (node indices) of each edge.
#' @return Integer vector of length \code{n} with component labels \code{1..m}.
#' @keywords internal
#' @export
uf_comp <- function(n, ei, ej) {
  parent <- seq_len(n)
  find <- function(x) {
    while (parent[x] != x) {
      parent[x] <<- parent[parent[x]]
      x <- parent[x]
    }
    x
  }
  if (length(ei)) {
    for (k in seq_along(ei)) {
      a <- find(ei[k])
      b <- find(ej[k])
      if (a != b) parent[a] <- b
    }
  }
  roots <- vapply(seq_len(n), find, integer(1))
  match(roots, unique(roots))
}

#' Snap sub-threshold band-pieces into their best contiguous neighbour
#'
#' Sequential snap/merge over a region-adjacency graph: repeatedly absorb each piece
#' whose (accumulated) area is below its band's floor into its contiguous neighbour,
#' preferring the highest band and then the largest area. The floor
#' \code{minac[band]} \strong{tightens as suitability rises}, so small high-value
#' ridge units survive while broad low aprons are cleaned of slivers. Passes repeat
#' until a pass makes no merges (or \code{npass} is reached), letting chains collapse.
#'
#' @param piece Integer vector (length = n cells) of band-piece labels \code{1..np}
#'   (e.g. from \code{\link{uf_comp}} on same-band edges).
#' @param ei,ej Integer vectors: \emph{all} rook edges as local cell-index pairs
#'   (from \code{\link{rook_edges}}); cross-piece edges define adjacency.
#' @param area Numeric vector (length \code{np}) of per-piece area (acres).
#' @param band Integer vector (length \code{np}) of per-piece band (1 low .. high).
#' @param minac Numeric vector of per-band minimum size (acres); a piece whose area
#'   is below \code{minac[band]} is a merge candidate.
#' @param npass Maximum number of snap passes.
#' @return Integer vector of length \code{np}: the resolved root piece for each piece.
#'   Map per-cell with \code{root[piece]}.
#' @keywords internal
#' @export
snap_pieces <- function(piece, ei, ej, area, band, minac, npass = 3L) {
  np <- length(area)

  # Region-adjacency graph: undirected neighbour list among pieces (cross-piece edges).
  dif <- piece[ei] != piece[ej]
  rag <- if (any(dif)) {
    unique(t(apply(cbind(piece[ei[dif]], piece[ej[dif]]), 1, sort)))
  } else {
    matrix(0L, 0, 2)
  }
  adj <- vector("list", np)
  if (nrow(rag)) {
    for (e in seq_len(nrow(rag))) {
      a <- rag[e, 1]
      b <- rag[e, 2]
      adj[[a]] <- c(adj[[a]], b)
      adj[[b]] <- c(adj[[b]], a)
    }
  }

  parent <- seq_len(np)
  arr <- area
  brr <- band
  find <- function(x) {
    while (parent[x] != x) {
      parent[x] <<- parent[parent[x]]
      x <- parent[x]
    }
    x
  }
  for (pass in seq_len(npass)) {
    did <- 0L
    for (p in order(arr)) {
      rp <- find(p)
      if (arr[rp] >= minac[brr[rp]]) {
        next
      }
      nb <- adj[[p]]
      if (is.null(nb)) {
        next
      }
      rr <- unique(vapply(nb, find, integer(1)))
      rr <- rr[rr != rp]
      if (!length(rr)) {
        next
      }
      tgt <- rr[order(-brr[rr], -arr[rr])][1]
      parent[rp] <- tgt
      arr[tgt] <- arr[tgt] + arr[rp]
      did <- did + 1L
    }
    if (did == 0L) break
  }
  vapply(seq_len(np), find, integer(1))
}

#' Geometric k-means cap on a cell set
#'
#' Split a set of cells into \code{ceil(n / capn)} compact pieces, each at or below
#' \code{capn} cells, by k-means on cell (col, row) positions; \code{k} is bumped up
#' until no cluster exceeds the cap (k-means targets equal sizes but does not bound
#' them). Sets under the cap pass through as a single piece.
#'
#' @param gc Integer vector of global (1-based) cell indices.
#' @param capn Maximum pieces size in cells.
#' @param nc Number of columns in the parent raster (to recover col/row from \code{gc}).
#' @return Integer cluster label per cell (all \code{1L} when \code{n <= capn}).
#' @keywords internal
#' @export
kcap <- function(gc, capn, nc) {
  n <- length(gc)
  if (n <= capn) {
    return(rep(1L, n))
  }
  xy <- cbind((gc - 1L) %% nc, (gc - 1L) %/% nc)
  k <- as.integer(ceiling(n / capn))
  repeat {
    km <- stats::kmeans(xy, centers = k, iter.max = 100)
    if (max(tabulate(km$cluster)) <= capn) {
      break
    }
    k <- k + 1L
  }
  km$cluster
}

# ---------------------------------------------------------------------------
# Stage 2: contour (band) -> snap -> k-means cap
# ---------------------------------------------------------------------------

#' Subdivide oversized patches by suitability banding
#'
#' For each patch larger than \code{trig_ac} acres: cut cells into fixed suitability
#' \emph{bands} on the raw surface, find same-band connected pieces, snap
#' sub-threshold pieces into their highest-band contiguous neighbour
#' (\code{\link{snap_pieces}}), then geometrically cap any surviving unit above
#' \code{cap_cells} (\code{\link{kcap}}). Patches at or below \code{trig_ac} pass
#' through as a single unit. Units are numbered sequentially across all patches.
#'
#' @param cells Integer vector of the global cell indices of all suitable cells
#'   (i.e. \code{which(!is.na(values(patch_id)))}).
#' @param memb Integer vector (same length as \code{cells}) of the parent patch id
#'   per suitable cell.
#' @param prv Numeric vector of probability values for the \emph{whole} raster
#'   (\code{values(pr)[, 1]}); indexed by global cell id.
#' @param nc,nr Columns and rows of the parent raster.
#' @param cell_area Area of one cell (m^2).
#' @param band_brks Numeric band edges on the raw surface (length = n bands + 1),
#'   passed to \code{cut(right = FALSE)}.
#' @param band_minac Numeric per-band minimum size (acres) for \code{\link{snap_pieces}}.
#' @param trig_ac Only patches larger than this (acres) are banded.
#' @param cap_cells Maximum final unit size in cells (geometric cap).
#' @param seed Random seed for the k-means cap (reproducible splits).
#' @return Integer vector (same length as \code{cells}) of subpatch unit ids.
#' @seealso \code{\link{delineate_patches}}, \code{\link{snap_pieces}}, \code{\link{kcap}}
#' @export
band_subpatches <- function(
  cells,
  memb,
  prv,
  nc,
  nr,
  cell_area,
  band_brks,
  band_minac,
  trig_ac,
  cap_cells,
  seed = 1
) {
  set.seed(seed)
  ac_per_m2 <- 1 / 4046.856
  allids <- as.integer(names(table(memb)))
  unit <- integer(length(cells))
  nu <- 0L

  for (P in allids) {
    sel <- which(memb == P)
    gc <- cells[sel]

    # Small patch -> one survey unit, untouched.
    if (length(gc) * cell_area * ac_per_m2 <= trig_ac) {
      nu <- nu + 1L
      unit[sel] <- nu
      next
    }

    # Contour into bands, then same-band connected pieces.
    vb <- as.integer(cut(prv[gc], band_brks, right = FALSE))
    E <- rook_edges(gc, nc, nr)
    ei <- E$i
    ej <- E$j
    same <- vb[ei] == vb[ej]
    piece <- uf_comp(length(gc), ei[same], ej[same])
    np <- max(piece)
    pband <- as.integer(tapply(vb, piece, `[`, 1)[as.character(seq_len(np))])
    parea <- as.numeric(tabulate(piece)) * cell_area * ac_per_m2

    # Snap sub-threshold pieces into their best neighbour.
    root <- snap_pieces(piece, ei, ej, parea, pband, band_minac)
    uroot <- root[piece]

    # Geometric cap on any surviving unit still too big.
    for (u in unique(uroot)) {
      iu <- which(uroot == u)
      if (length(iu) <= cap_cells) {
        nu <- nu + 1L
        unit[sel[iu]] <- nu
      } else {
        cl <- kcap(gc[iu], cap_cells, nc)
        for (k in unique(cl)) {
          nu <- nu + 1L
          unit[sel[iu[cl == k]]] <- nu
        }
      }
    }
  }
  unit
}
