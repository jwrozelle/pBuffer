#' Generate a Monte Carlo displacement cloud under DHS-style rules
#'
#' Creates a cloud of candidate locations for a single DHS cluster by drawing
#' many displacement realizations consistent with DHS public-use rules:
#' urban: max 2 km; rural: max 5 km with a small fraction (rounded down) allowed up to 10 km.
#'
#' The output is an sf POINT object that can be fed into pBuffer join functions
#' directly, or binned into a grid via `pb_cloud_to_grid()`.
#'
#' @param point_sf A single-row `sf` POINT object (one cluster). Must include `urban_col`.
#' @param n_draws Number of Monte Carlo draws (default 200000).
#' @param urban_col Name of urban/rural column (default `"URBAN_RURA"`, expecting `"U"`/`"R"`).
#' @param p_rural_hi Fraction of rural draws with higher max displacement (default 0.01).
#'   DHS protocol specifies 1% of rural clusters, rounded down; implemented as `floor(n_draws*p_rural_hi)`.
#'   Set to 0 to ignore the rare 10 km tail.
#' @param rural_max Rural max displacement in meters (default 5000).
#' @param rural_max_hi Rare rural max displacement in meters (default 10000).
#' @param urban_max Urban max displacement in meters (default 2000).
#' @param uniform_area If FALSE (default), distances are uniform in radius. If TRUE, uniform over area.
#' @param weightsCol Name of the weight column in the output (default `"layer"`).
#' @param keep_offsets If TRUE, includes `dx` and `dy` columns in output (default TRUE).
#'
#' @return An `sf` POINT object with `n_draws` rows. Contains a weight column (`weightsCol`)
#'   summing to 1 and geometry representing the displaced locations.
#'
#' @export
pb_cloud_dhs_mc <- function(point_sf,
                            n_draws = 200000,
                            urban_col = "URBAN_RURA",
                            p_rural_hi = 0.01,
                            rural_max = 5000,
                            rural_max_hi = 10000,
                            urban_max = 2000,
                            uniform_area = FALSE,
                            weightsCol = "layer",
                            keep_offsets = TRUE) {
  
  if (!inherits(point_sf, "sf")) stop("point_sf must be an sf object.")
  if (nrow(point_sf) != 1) stop("point_sf must have exactly 1 row.")
  if (!inherits(sf::st_geometry(point_sf), "sfc_POINT")) stop("point_sf must have POINT geometry.")
  if (!urban_col %in% names(point_sf)) stop("Missing urban/rural column: ", urban_col)
  
  crs_in <- sf::st_crs(point_sf)
  if (is.null(crs_in)) stop("point_sf has no CRS; set it before using this function.")
  xy <- sf::st_coordinates(point_sf)[1, ]
  if (anyNA(xy)) stop("Could not read coordinates from point_sf.")
  
  ur <- point_sf[[urban_col]]
  if (!ur %in% c("U", "R")) stop("urban_col must contain 'U' or 'R' for the provided point.")
  
  n <- as.integer(n_draws)
  if (n < 1) stop("n_draws must be >= 1.")
  
  # Draw angles
  theta <- stats::runif(n, 0, 2 * pi)
  
  # Determine max radius per draw
  if (ur == "U") {
    max_r <- rep(urban_max, n)
  } else {
    max_r <- rep(rural_max, n)
    
    # DHS: 1% of rural clusters, rounded down (generalized via p_rural_hi)
    if (is.finite(p_rural_hi) && p_rural_hi > 0 && rural_max_hi > rural_max) {
      k <- as.integer(floor(n * p_rural_hi))
      if (k >= 1L) {
        idx <- sample.int(n, size = k, replace = FALSE)
        max_r[idx] <- rural_max_hi
      }
    }
  }
  
  u <- stats::runif(n, 0, 1)
  r <- if (uniform_area) sqrt(u) * max_r else u * max_r
  
  dx <- cos(theta) * r
  dy <- sin(theta) * r
  
  # Build sf points (absolute coords)
  df <- data.frame(
    draw_id = seq_len(n),
    dx = dx,
    dy = dy,
    x = xy[1] + dx,
    y = xy[2] + dy
  )
  
  # Equal weights
  df[[weightsCol]] <- rep(1 / n, n)
  
  out <- sf::st_as_sf(df, coords = c("x", "y"), crs = crs_in)
  
  if (!keep_offsets) {
    out$dx <- NULL
    out$dy <- NULL
  }
  
  out
}


#' Bin a displacement cloud into a probability grid (raster)
#'
#' Converts a set of candidate point locations (an sf POINT cloud) into a discrete
#' probability surface on a regular grid. By default this uses simple cell binning
#' (counts) and then normalizes to sum to 1. If a weights column is provided, cell
#' values are the sum of weights.
#'
#' Optionally trims/masks the grid to an administrative boundary and renormalizes
#' via `trimProbBuff()`.
#'
#' @param cloud_sf An `sf` POINT object (e.g., output of `pb_cloud_dhs_mc()`).
#' @param cell_m Grid resolution in meters (default 50).
#' @param weightsCol Name of the weight column in `cloud_sf` (default `"layer"`). If missing,
#'   points are treated as equal weight.
#' @param extent_m Optional numeric scalar giving half-width of the square extent centered on
#'   `center_xy`. If NULL, extent is derived from the cloud bounding box.
#' @param center_xy Optional numeric length-2 vector c(x,y) for grid centering. If NULL,
#'   uses the centroid of the cloud’s bbox.
#' @param adminBound Optional `sf` polygon used to trim/mask the grid (default NULL).
#' @param drop_zeros If TRUE (default), zeros remain zeros in raster; conversion to points
#'   can drop them later.
#'
#' @return A `raster::RasterLayer` whose cell values sum to 1.
#'
#' @export
pb_cloud_to_grid <- function(cloud_sf,
                             cell_m = 50,
                             weightsCol = "layer",
                             extent_m = NULL,
                             center_xy = NULL,
                             adminBound = NULL,
                             drop_zeros = TRUE) {
  
  if (!inherits(cloud_sf, "sf")) stop("cloud_sf must be an sf object.")
  if (!inherits(sf::st_geometry(cloud_sf), "sfc_POINT")) stop("cloud_sf must have POINT geometry.")
  if (nrow(cloud_sf) < 1) stop("cloud_sf has zero rows.")
  
  crs_in <- sf::st_crs(cloud_sf)
  if (is.null(crs_in)) stop("cloud_sf has no CRS; set it before using this function.")
  
  # Determine grid center
  bb <- sf::st_bbox(cloud_sf)
  if (is.null(center_xy)) {
    center_xy <- c((bb["xmin"] + bb["xmax"]) / 2, (bb["ymin"] + bb["ymax"]) / 2)
  }
  
  # Determine extent
  if (is.null(extent_m)) {
    extent_m <- max(bb["xmax"] - bb["xmin"], bb["ymax"] - bb["ymin"]) / 2
  }
  
  r0 <- raster::raster(
    xmn = center_xy[1] - extent_m, xmx = center_xy[1] + extent_m,
    ymn = center_xy[2] - extent_m, ymx = center_xy[2] + extent_m,
    res = cell_m,
    crs = crs_in
  )
  
  # Rasterize: sum weights if present; else count points
  if (weightsCol %in% names(cloud_sf)) {
    sp_pts <- as(cloud_sf, "Spatial")
    rr <- raster::rasterize(sp_pts, r0, field = weightsCol, fun = "sum", background = 0)
  } else {
    sp_pts <- as(cloud_sf, "Spatial")
    rr <- raster::rasterize(sp_pts, r0, field = 1, fun = "count", background = 0)
  }
  
  # Optional trim + renormalize
  if (!is.null(adminBound)) {
    if (!inherits(adminBound, "sf")) stop("adminBound must be an sf object if provided.")
    rr <- trimProbBuff(probRaster = rr, adminBound = adminBound)
  } else {
    s <- sum(rr[], na.rm = TRUE)
    if (is.na(s) || s == 0) stop("Grid sums to 0/NA; cannot normalize.")
    rr[] <- rr[] / s
  }
  
  rr
}


#' Convert a probability grid (raster) to weighted points
#'
#' Converts a probability raster (e.g., output of `pb_cloud_to_grid()`) into an `sf`
#' POINT layer located at cell centers, with a weight column containing the raster
#' cell values.
#'
#' @param prob_raster A `raster::RasterLayer` whose values represent probabilities.
#' @param weightsCol Name of output weight column (default `"layer"`).
#' @param drop_zeros If TRUE (default), drops cells with zero/NA probability.
#'
#' @return An `sf` POINT object with weight column `weightsCol` summing to 1.
#' @export
pb_grid_to_points <- function(prob_raster,
                              weightsCol = "layer",
                              drop_zeros = TRUE) {
  
  if (!inherits(prob_raster, "RasterLayer")) stop("prob_raster must be a raster::RasterLayer.")
  
  df <- raster::rasterToPoints(prob_raster) |> as.data.frame()
  val_col <- setdiff(names(df), c("x", "y"))
  if (length(val_col) != 1) stop("Unexpected rasterToPoints output (expected one value column).")
  
  if (drop_zeros) {
    df <- df[!is.na(df[[val_col]]) & df[[val_col]] > 0, , drop = FALSE]
  }
  
  pts <- sf::st_as_sf(df, coords = c("x", "y"), crs = sf::st_crs(prob_raster))
  names(pts)[names(pts) == val_col] <- weightsCol
  
  # Renormalize (important if dropped zeros)
  s <- sum(pts[[weightsCol]], na.rm = TRUE)
  if (is.na(s) || s == 0) stop("Point weights sum to 0/NA after conversion.")
  pts[[weightsCol]] <- pts[[weightsCol]] / s
  
  pts
}