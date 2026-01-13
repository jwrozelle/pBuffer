# ======================================================================
# prob_points_mc.R
#
# Purpose:
#   Monte Carlo utilities for DHS-style positional uncertainty:
#     - pb_cloud_dhs_mc(): simulate a cloud of candidate locations around a point
#     - pb_cloud_to_grid(): bin cloud to a probability raster (sums to 1)
#     - pb_grid_to_points(): convert probability raster to weighted points
#
# Notes:
#   These are useful when you want an MC-derived approximation to the density buffer,
#   or when you want to empirically approximate complex trimming constraints.
# ======================================================================

#' Generate a Monte Carlo displacement cloud under DHS-style rules
#'
#' @description
#' Creates a cloud of candidate locations for a single cluster by drawing many
#' displacement realizations:
#' \itemize{
#'   \item urban: 0–2km
#'   \item rural: 99% 0–5km, 1% 0–10km (configurable)
#' }
#'
#' This cloud can be binned to a grid using \code{pb_cloud_to_grid()}.
#'
#' @param point_sf Single-row \code{sf} POINT object.
#' @param n_draws Number of Monte Carlo draws.
#' @param urban_col Name of urban/rural column (expects "U"/"R").
#' @param p_rural_hi Tail probability for rural 10km draws.
#' @param rural_max Rural max displacement (m).
#' @param rural_max_hi Rural tail max displacement (m).
#' @param urban_max Urban max displacement (m).
#' @param uniform_area If TRUE, draws uniform in area: r = sqrt(U)*R. If FALSE, uniform in radius: r = U*R.
#' @param weightsCol Output weight column name (default "layer").
#' @param keep_offsets If TRUE, keep dx/dy columns for debugging.
#'
#' @return \code{sf} POINT with \code{n_draws} rows and weights summing to 1.
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
  
  if (!inherits(point_sf, "sf")) stop("point_sf must be an sf object.", call. = FALSE)
  if (nrow(point_sf) != 1) stop("point_sf must have exactly 1 row.", call. = FALSE)
  if (!inherits(sf::st_geometry(point_sf), "sfc_POINT")) stop("point_sf must have POINT geometry.", call. = FALSE)
  if (!urban_col %in% names(point_sf)) stop("Missing urban/rural column: ", urban_col, call. = FALSE)
  
  crs_in <- sf::st_crs(point_sf)
  if (is.null(crs_in)) stop("point_sf has no CRS; set it before use.", call. = FALSE)
  
  xy <- sf::st_coordinates(point_sf)[1, ]
  if (anyNA(xy)) stop("Could not read coordinates from point_sf.", call. = FALSE)
  
  ur <- as.character(sf::st_drop_geometry(point_sf)[[urban_col]][1])
  if (!ur %in% c("U", "R")) stop("urban_col must be 'U' or 'R' for the provided point.", call. = FALSE)
  
  # Determine max radius per draw under DHS mixture
  if (ur == "U") {
    max_m <- rep(urban_max, n_draws)
  } else {
    # Implement tail as a fraction of draws
    n_hi <- floor(n_draws * p_rural_hi)
    max_m <- c(rep(rural_max_hi, n_hi), rep(rural_max, n_draws - n_hi))
    max_m <- sample(max_m)  # shuffle
  }
  
  theta <- stats::runif(n_draws, 0, 2 * pi)
  u     <- stats::runif(n_draws, 0, 1)
  r     <- if (isTRUE(uniform_area)) sqrt(u) * max_m else u * max_m
  
  dx <- cos(theta) * r
  dy <- sin(theta) * r
  
  out <- data.frame(
    x = xy["X"] + dx,
    y = xy["Y"] + dy,
    stringsAsFactors = FALSE
  )
  
  # Equal-weight cloud; later binning turns this into probability mass
  out[[weightsCol]] <- rep(1 / n_draws, n_draws)
  
  if (isTRUE(keep_offsets)) {
    out$dx <- dx
    out$dy <- dy
    out$max_m <- max_m
  }
  
  sf::st_as_sf(out, coords = c("x", "y"), crs = crs_in)
}

#' Bin a Monte Carlo cloud to a probability raster
#'
#' @description
#' Converts a cloud of candidate points (with optional weights) into a RasterLayer
#' where cell values represent probability mass and sum to 1 (after optional trimming).
#'
#' @param cloud_sf \code{sf} POINT cloud.
#' @param cell_m Raster resolution in meters.
#' @param weightsCol Weight column name in \code{cloud_sf}. If missing, points treated as equal weight.
#' @param adminBound Optional \code{sf} polygon to mask/trim the probability raster.
#' @param adminID Admin ID column in adminBound.
#'
#' @return \code{raster::RasterLayer} whose cell values sum to 1.
#' @export
pb_cloud_to_grid <- function(cloud_sf,
                             cell_m = 50,
                             weightsCol = "layer",
                             adminBound = NULL,
                             adminID = "ID_2") {
  
  if (!inherits(cloud_sf, "sf")) stop("cloud_sf must be an sf object.", call. = FALSE)
  if (!inherits(sf::st_geometry(cloud_sf), "sfc_POINT")) stop("cloud_sf must have POINT geometry.", call. = FALSE)
  
  crs_in <- sf::st_crs(cloud_sf)
  if (is.null(crs_in)) stop("cloud_sf has no CRS.", call. = FALSE)
  
  w <- NULL
  if (weightsCol %in% names(cloud_sf)) {
    w <- sf::st_drop_geometry(cloud_sf)[[weightsCol]]
  } else {
    w <- rep(1, nrow(cloud_sf))
  }
  w[is.na(w)] <- 0
  
  # Build template raster over cloud extent
  bb <- sf::st_bbox(cloud_sf)
  r <- raster::raster(
    raster::extent(bb["xmin"], bb["xmax"], bb["ymin"], bb["ymax"]),
    res = cell_m,
    crs = crs_in$wkt
  )
  
  # Rasterize weighted points to probability mass per cell
  sp_pts <- methods::as(cloud_sf, "Spatial")
  r <- raster::rasterize(sp_pts, r, field = w, fun = sum, background = 0)
  
  # Optional admin trimming
  if (!is.null(adminBound)) {
    if (!inherits(adminBound, "sf")) stop("adminBound must be sf.", call. = FALSE)
    if (!adminID %in% names(adminBound)) stop("adminBound missing adminID.", call. = FALSE)
    if (sf::st_crs(adminBound) != sf::st_crs(cloud_sf)) {
      stop("adminBound and cloud_sf must share CRS.", call. = FALSE)
    }
    
    r <- suppressWarnings(raster::mask(r, methods::as(adminBound, "Spatial")))
  }
  
  # Renormalize to sum to 1 over non-NA
  vals <- raster::values(r)
  ok <- !is.na(vals)
  s <- sum(vals[ok])
  if (is.finite(s) && s > 0) {
    vals[ok] <- vals[ok] / s
    raster::values(r) <- vals
  } else {
    raster::values(r) <- NA_real_
  }
  
  r
}

#' Convert a probability raster to weighted points
#'
#' @description
#' Converts a RasterLayer of probability mass into an sf point set located at cell centers,
#' with a weight column equal to the cell probability.
#'
#' @param prob_raster RasterLayer with probability mass per cell (sums to 1).
#' @param weightsCol Output weight column name (default "layer").
#' @param drop_zeros If TRUE, drops zero and NA cells.
#'
#' @return \code{sf} POINT with one row per nonzero cell.
#' @export
pb_grid_to_points <- function(prob_raster,
                              weightsCol = "layer",
                              drop_zeros = TRUE) {
  
  if (!inherits(prob_raster, "RasterLayer")) stop("prob_raster must be a RasterLayer.", call. = FALSE)
  
  pts <- raster::rasterToPoints(prob_raster, spatial = TRUE)
  if (!inherits(pts, "SpatialPointsDataFrame")) stop("Unexpected rasterToPoints output.", call. = FALSE)
  
  pts <- sf::st_as_sf(pts)
  if (ncol(sf::st_drop_geometry(pts)) != 1) {
    stop("Probability raster must have exactly one value layer.", call. = FALSE)
  }
  
  names(pts)[1] <- weightsCol
  
  if (isTRUE(drop_zeros)) {
    v <- sf::st_drop_geometry(pts)[[weightsCol]]
    pts <- pts[!is.na(v) & v > 0, ]
  }
  
  # Ensure weights sum to 1 (numerical safety)
  v <- sf::st_drop_geometry(pts)[[weightsCol]]
  s <- sum(v, na.rm = TRUE)
  if (is.finite(s) && s > 0) {
    pts[[weightsCol]] <- v / s
  }
  
  pts
}
