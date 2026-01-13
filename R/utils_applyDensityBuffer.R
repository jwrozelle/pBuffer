# ======================================================================
# apply_densityBuffer.R
#
# Purpose:
#   Internal helper that applies a pb_densityBuffer to a single displaced point,
#   producing an aligned weight surface (RasterLayer) optionally trimmed to an
#   administrative boundary and renormalized to sum to 1.
#
# Why this exists:
#   Every join type (raster/polygon/point-count/nearest-ID probabilities) needs
#   exactly the same "apply kernel at point" logic. Centralizing that logic avoids
#   divergence and makes mixture support trivial across join types.
#
# Key responsibilities:
#   - choose correct kernel(s): single radius or mixture
#   - handle rural_only mixture by checking urban/rural label
#   - align origin-centered kernel to the point location
#   - optionally mask to admin polygon containing the point
#   - renormalize weights after trimming
# ======================================================================

#' Apply a density buffer at a point to produce an aligned weight raster
#'
#' @param point_sf Single-row sf POINT.
#' @param densityBuffer pb_densityBuffer object.
#' @param templateRaster RasterLayer used for CRS and alignment context (the join raster).
#' @param radius_m Optional override radius (single-kernel mode).
#' @param adminBound Optional admin polygons; if provided weights trimmed to polygon.
#' @param adminID Admin ID column in adminBound.
#'
#' @return A RasterLayer with weights (summing to 1 over non-NA cells), aligned near the point.
#' @keywords internal
pb_apply_densityBuffer <- function(point_sf,
                                   densityBuffer,
                                   templateRaster,
                                   radius_m = NULL,
                                   adminBound = NULL,
                                   adminID = "ID_2") {
  
  if (!inherits(point_sf, "sf") || nrow(point_sf) != 1) {
    stop("point_sf must be a single-row sf POINT.", call. = FALSE)
  }
  if (!inherits(densityBuffer, "pb_densityBuffer")) {
    stop("densityBuffer must be a pb_densityBuffer object.", call. = FALSE)
  }
  if (!inherits(templateRaster, "RasterLayer")) {
    stop("templateRaster must be a raster::RasterLayer.", call. = FALSE)
  }
  
  xy <- sf::st_coordinates(point_sf)
  if (nrow(xy) != 1 || anyNA(xy[1, c("X","Y")])) return(NULL)
  
  # Helper: determine whether point is urban (affects rural_only mixture)
  urban_col <- densityBuffer$urban_col %||% "URBAN_RURA"
  urban_val <- densityBuffer$urban_value %||% "U"
  is_urban <- FALSE
  if (urban_col %in% names(point_sf)) {
    v <- as.character(sf::st_drop_geometry(point_sf)[[urban_col]][1])
    is_urban <- isTRUE(!is.na(v) && v == urban_val)
  }
  
  # Determine kernel components:
  # - default: single radius kernel
  # - if mixture enabled and applicable: mixture radii/weights
  comps <- NULL
  
  mix <- densityBuffer$mixture
  if (!is.null(mix) && isTRUE(mix$enabled)) {
    
    applies_to <- mix$applies_to %||% "rural_only"
    use_mix <- if (applies_to == "global") TRUE else !is_urban
    
    if (use_mix) {
      comps <- Map(function(w, r) list(weight = w, radius_m = r),
                   mix$weights, mix$radii_m)
    }
  }
  
  if (is.null(comps)) {
    # Single-kernel mode
    if (is.null(radius_m)) radius_m <- max(densityBuffer$radii_m)
    comps <- list(list(weight = 1, radius_m = radius_m))
  }
  
  # Compute a weighted sum of aligned kernel rasters.
  # This avoids duplicating mixture logic in every join function.
  w_raster_sum <- NULL
  
  for (comp in comps) {
    
    r_key <- as.character(comp$radius_m)
    if (!r_key %in% names(densityBuffer$kernels)) {
      stop("densityBuffer missing kernel for radius ", comp$radius_m, call. = FALSE)
    }
    
    Kmat <- densityBuffer$kernels[[r_key]]
    
    # Build origin-centered kernel raster and place it at the point location.
    pb_r <- raster::raster(Kmat)
    raster::crs(pb_r) <- raster::crs(templateRaster)
    
    cell_m <- densityBuffer$cell_m
    
    # Center kernel extent on point coordinate
    pb_xmin <- xy[1, "X"] - cell_m * (ncol(Kmat) / 2)
    pb_xmax <- xy[1, "X"] + cell_m * (ncol(Kmat) / 2)
    pb_ymin <- xy[1, "Y"] - cell_m * (nrow(Kmat) / 2)
    pb_ymax <- xy[1, "Y"] + cell_m * (nrow(Kmat) / 2)
    raster::extent(pb_r) <- c(pb_xmin, pb_xmax, pb_ymin, pb_ymax)
    
    # Weight by mixture component weight
    pb_r <- pb_r * comp$weight
    
    w_raster_sum <- if (is.null(w_raster_sum)) pb_r else (w_raster_sum + pb_r)
  }
  
  # Optional admin trimming: mask weights to polygon containing the point
  if (!is.null(adminBound)) {
    if (!inherits(adminBound, "sf")) stop("adminBound must be sf.", call. = FALSE)
    if (!adminID %in% names(adminBound)) stop("adminBound missing adminID.", call. = FALSE)
    if (sf::st_crs(adminBound) != sf::st_crs(point_sf)) {
      stop("adminBound and point_sf must share CRS.", call. = FALSE)
    }
    
    j <- suppressWarnings(sf::st_join(point_sf, adminBound[, adminID, drop = FALSE], sf::st_within, left = TRUE))
    a_val <- sf::st_drop_geometry(j)[[adminID]][1]
    
    if (!is.na(a_val)) {
      poly <- adminBound[adminBound[[adminID]] == a_val, ]
      w_raster_sum <- suppressWarnings(raster::mask(w_raster_sum, methods::as(poly, "Spatial")))
    }
  }
  
  # Renormalize weights to sum to 1 over non-NA cells
  w <- raster::values(w_raster_sum)
  ok <- !is.na(w)
  s <- sum(w[ok])
  if (!is.finite(s) || s <= 0) return(NULL)
  
  raster::values(w_raster_sum)[ok] <- w[ok] / s
  w_raster_sum
}
