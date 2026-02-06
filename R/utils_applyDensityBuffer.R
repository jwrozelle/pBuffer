# ======================================================================
# utils_applyDensityBuffer.R
#
# Purpose:
#   Internal helper that applies a pb_densityBuffer to displaced point(s),
#   producing aligned weight surface(s) (RasterLayer) optionally trimmed to an
#   administrative boundary and renormalized to sum to 1.
#
# Critical behavior (used across join_*):
#   - If point_sf has 1 row: returns a single RasterLayer (or NULL)
#   - If point_sf has >1 row: returns a list of RasterLayer/NULL, length nrow(point_sf)
#
# Why:
#   pb_pointProbJoin() and pb_polyProbsJoin() need per-point kernels for many
#   points; pb_rasterValJoin() calls this for one point at a time.
# ======================================================================

# Local fallback: utils_validation.R also defines this, but keeping a guard here
# prevents runtime surprises if file load order changes.
if (!exists("%||%", mode = "function")) {
  `%||%` <- function(x, y) if (is.null(x)) y else x
}

#' Apply a density buffer at point(s) to produce aligned weight raster(s)
#'
#' @param point_sf sf POINT (one or more rows).
#' @param densityBuffer pb_densityBuffer object.
#' @param templateRaster Optional raster::RasterLayer used for CRS context.
#'   If NULL, a minimal in-memory RasterLayer is created using point_sf's CRS.
#' @param radius_m Optional override radius (single-kernel mode).
#' @param adminBound Optional admin polygons; if provided weights trimmed to polygon.
#' @param adminID Admin ID column in adminBound.
#'
#' @return RasterLayer (if nrow(point_sf)==1) or list of RasterLayer/NULL.
#' @keywords internal
pb_apply_densityBuffer <- function(point_sf,
                                   densityBuffer,
                                   templateRaster = NULL,
                                   radius_m = NULL,
                                   adminBound = NULL,
                                   adminID = "ID_2") {
  
  if (!inherits(point_sf, "sf")) {
    stop("point_sf must be an sf object.", call. = FALSE)
  }
  if (nrow(point_sf) < 1) {
    stop("point_sf must have at least 1 row.", call. = FALSE)
  }
  if (!inherits(densityBuffer, "pb_densityBuffer")) {
    stop("densityBuffer must be a pb_densityBuffer object.", call. = FALSE)
  }
  
  # We only need a CRS context to stamp onto the kernel raster.
  # For raster joins, callers pass the input raster (preferred).
  if (is.null(templateRaster)) {
    crs_wkt <- tryCatch(sf::st_crs(point_sf)$wkt, error = function(e) NA_character_)
    if (is.na(crs_wkt) || !nzchar(crs_wkt)) {
      stop("templateRaster is NULL and point_sf has no CRS; cannot set CRS for kernel rasters.", call. = FALSE)
    }
    templateRaster <- raster::raster()
    raster::crs(templateRaster) <- crs_wkt
  } else {
    if (!inherits(templateRaster, "RasterLayer")) {
      stop("templateRaster must be a raster::RasterLayer (or NULL).", call. = FALSE)
    }
  }
  
  # Internal single-row worker
  one_point <- function(pt) {
    
    xy <- sf::st_coordinates(pt)
    if (nrow(xy) != 1 || anyNA(xy[1, c("X", "Y")])) return(NULL)
    
    # Determine whether point is urban (affects rural_only mixture)
    urban_col <- densityBuffer$urban_col %||% "URBAN_RURA"
    urban_val <- densityBuffer$urban_value %||% "U"
    is_urban <- FALSE
    if (urban_col %in% names(pt)) {
      v <- as.character(sf::st_drop_geometry(pt)[[urban_col]][1])
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
        comps <- Map(
          function(w, r) list(weight = w, radius_m = r),
          mix$weights, mix$radii_m
        )
      }
    }
    
    if (is.null(comps)) {
      # Single-kernel mode
      if (is.null(radius_m)) radius_m <- max(densityBuffer$radii_m)
      comps <- list(list(weight = 1, radius_m = radius_m))
    }
    
    # Weighted sum of aligned kernel rasters
    w_raster_sum <- NULL
    cell_m <- densityBuffer$cell_m
    
    for (comp in comps) {
      
      r_key <- as.character(comp$radius_m)
      if (!r_key %in% names(densityBuffer$kernels)) {
        stop("densityBuffer missing kernel for radius ", comp$radius_m, call. = FALSE)
      }
      
      Kmat <- densityBuffer$kernels[[r_key]]
      
      pb_r <- raster::raster(Kmat)
      raster::crs(pb_r) <- raster::crs(templateRaster)
      
      # Center kernel extent on point coordinate
      pb_xmin <- xy[1, "X"] - cell_m * (ncol(Kmat) / 2)
      pb_xmax <- xy[1, "X"] + cell_m * (ncol(Kmat) / 2)
      pb_ymin <- xy[1, "Y"] - cell_m * (nrow(Kmat) / 2)
      pb_ymax <- xy[1, "Y"] + cell_m * (nrow(Kmat) / 2)
      raster::extent(pb_r) <- c(pb_xmin, pb_xmax, pb_ymin, pb_ymax)
      
      pb_r <- pb_r * comp$weight
      w_raster_sum <- if (is.null(w_raster_sum)) pb_r else (w_raster_sum + pb_r)
    }
    
    # Optional admin trimming: mask weights to polygon containing the point
    if (!is.null(adminBound)) {
      if (!inherits(adminBound, "sf")) stop("adminBound must be sf.", call. = FALSE)
      if (!adminID %in% names(adminBound)) stop("adminBound missing adminID.", call. = FALSE)
      if (sf::st_crs(adminBound) != sf::st_crs(pt)) {
        stop("adminBound and point_sf must share CRS.", call. = FALSE)
      }
      
      j <- suppressWarnings(sf::st_join(pt, adminBound[, adminID, drop = FALSE], sf::st_within, left = TRUE))
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
  
  # Vectorized behavior: list for multi-row, RasterLayer for single-row
  if (nrow(point_sf) == 1) {
    return(one_point(point_sf))
  }
  
  out <- vector("list", nrow(point_sf))
  for (i in seq_len(nrow(point_sf))) {
    out[[i]] <- one_point(point_sf[i, , drop = FALSE])
  }
  out
}
