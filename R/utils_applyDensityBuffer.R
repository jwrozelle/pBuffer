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
                                   templateRaster,
                                   radius_m = NULL,
                                   adminBound = NULL,
                                   adminID = "ID_2") {
  
  if (!inherits(point_sf, "sf") || nrow(point_sf) != 1L) {
    stop("point_sf must be a single-row sf POINT.", call. = FALSE)
  }
  if (!inherits(densityBuffer, "pb_densityBuffer")) {
    stop("densityBuffer must be a pb_densityBuffer object.", call. = FALSE)
  }
  if (!inherits(templateRaster, "RasterLayer")) {
    stop("templateRaster must be a raster::RasterLayer.", call. = FALSE)
  }
  
  xy <- sf::st_coordinates(point_sf)
  if (nrow(xy) != 1L || anyNA(xy[1, c("X", "Y")])) return(NULL)
  
  urban_col   <- if (is.null(densityBuffer$urban_col)) "URBAN_RURA" else densityBuffer$urban_col
  urban_value <- if (is.null(densityBuffer$urban_value)) "U" else densityBuffer$urban_value
  
  is_urban <- FALSE
  if (urban_col %in% names(point_sf)) {
    v <- as.character(sf::st_drop_geometry(point_sf)[[urban_col]][1])
    is_urban <- isTRUE(!is.na(v) && v == urban_value)
  }
  
  comps <- NULL
  mix <- densityBuffer$mixture
  
  if (!is.null(mix) && isTRUE(mix$enabled)) {
    applies_to <- if (is.null(mix$applies_to)) "rural_only" else mix$applies_to
    use_mix <- if (identical(applies_to, "global")) TRUE else !is_urban
    
    if (use_mix) {
      comps <- Map(
        function(w, r) list(weight = w, radius_m = r),
        mix$weights,
        mix$radii_m
      )
    }
  }
  
  if (is.null(comps)) {
    if (is.null(radius_m)) radius_m <- max(densityBuffer$radii_m)
    comps <- list(list(weight = 1, radius_m = radius_m))
  }
  
  build_component <- function(comp) {
    r_key <- as.character(comp$radius_m)
    if (!r_key %in% names(densityBuffer$kernels)) {
      stop("densityBuffer missing kernel for radius ", comp$radius_m, call. = FALSE)
    }
    
    Kmat <- densityBuffer$kernels[[r_key]]
    if (!is.matrix(Kmat)) Kmat <- as.matrix(Kmat)
    
    pb_r <- raster::raster(Kmat)
    raster::crs(pb_r) <- raster::crs(templateRaster)
    
    cell_m <- densityBuffer$cell_m
    raster::extent(pb_r) <- raster::extent(
      xy[1, "X"] - cell_m * (ncol(Kmat) / 2),
      xy[1, "X"] + cell_m * (ncol(Kmat) / 2),
      xy[1, "Y"] - cell_m * (nrow(Kmat) / 2),
      xy[1, "Y"] + cell_m * (nrow(Kmat) / 2)
    )
    
    pb_r * comp$weight
  }
  
  parts <- lapply(comps, build_component)
  
  if (length(parts) == 1L) {
    w_raster_sum <- parts[[1]]
  } else {
    exts <- lapply(parts, raster::extent)
    
    union_ext <- raster::extent(
      min(vapply(exts, function(e) e@xmin, numeric(1))),
      max(vapply(exts, function(e) e@xmax, numeric(1))),
      min(vapply(exts, function(e) e@ymin, numeric(1))),
      max(vapply(exts, function(e) e@ymax, numeric(1)))
    )
    
    parts <- lapply(parts, raster::extend, y = union_ext, value = 0)
    w_raster_sum <- Reduce(`+`, parts)
  }
  
  if (!is.null(adminBound)) {
    if (!inherits(adminBound, "sf")) stop("adminBound must be sf.", call. = FALSE)
    if (!adminID %in% names(adminBound)) stop("adminBound missing adminID.", call. = FALSE)
    if (sf::st_crs(adminBound) != sf::st_crs(point_sf)) {
      stop("adminBound and point_sf must share CRS.", call. = FALSE)
    }
    
    poly <- NULL
    
    if (nrow(adminBound) == 1L) {
      poly <- adminBound
    } else {
      point_for_join <- point_sf[, setdiff(names(point_sf), adminID), drop = FALSE]
      
      j <- suppressWarnings(
        sf::st_join(
          point_for_join,
          adminBound[, adminID, drop = FALSE],
          join = sf::st_within,
          left = TRUE
        )
      )
      
      vals <- sf::st_drop_geometry(j)[[adminID]]
      a_val <- if (length(vals)) vals[1] else NA
      
      if (!is.na(a_val)) {
        poly <- adminBound[adminBound[[adminID]] == a_val, , drop = FALSE]
      }
    }
    
    if (!is.null(poly) && nrow(poly) > 0) {
      w_raster_sum <- suppressWarnings(
        raster::mask(w_raster_sum, methods::as(poly, "Spatial"))
      )
    }
  }
  
  w <- raster::values(w_raster_sum)
  ok <- !is.na(w)
  s <- sum(w[ok])
  
  if (!is.finite(s) || s <= 0) return(NULL)
  
  raster::values(w_raster_sum)[ok] <- w[ok] / s
  w_raster_sum
}
