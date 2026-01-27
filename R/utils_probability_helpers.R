#' @keywords internal
pb__weights_raster_to_points <- function(w_raster,
                                         weight_name = "w",
                                         drop_zeros = TRUE,
                                         min_weight = 0) {
  # Convert a RasterLayer of weights to an sf POINT object with a weight column.
  # Uses cell centers as candidate locations.
  #
  # w_raster: RasterLayer with weights (probability mass or density*cell_area).
  # drop_zeros: if TRUE, keep only strictly positive (or > min_weight) weights.
  #
  # Returns: sf POINT with columns: x, y, w
  
  stopifnot(inherits(w_raster, "RasterLayer"))
  
  # raster::rasterToPoints is fast and stable for this use.
  pts_mat <- raster::rasterToPoints(w_raster, spatial = FALSE) # columns: x, y, layer
  if (!is.matrix(pts_mat) || nrow(pts_mat) == 0) {
    # Return empty sf with correct columns
    empty <- sf::st_as_sf(
      data.frame(x = numeric(0), y = numeric(0), w = numeric(0)),
      coords = c("x", "y"),
      crs = raster::crs(w_raster)
    )
    names(empty)[names(empty) == "w"] <- weight_name
    return(empty)
  }
  
  colnames(pts_mat) <- c("x", "y", "w")
  
  if (drop_zeros) {
    pts_mat <- pts_mat[is.finite(pts_mat[, "w"]) & pts_mat[, "w"] > min_weight, , drop = FALSE]
  } else {
    pts_mat <- pts_mat[is.finite(pts_mat[, "w"]), , drop = FALSE]
  }
  
  pts_df <- as.data.frame(pts_mat)
  names(pts_df)[names(pts_df) == "w"] <- weight_name
  
  sf::st_as_sf(pts_df, coords = c("x", "y"), crs = raster::crs(w_raster))
}

#' @keywords internal
pb__weighted_quantile <- function(x, w, probs = c(0.5)) {
  # Weighted quantiles with deterministic behavior.
  # x: numeric vector
  # w: nonnegative numeric weights (not necessarily normalized)
  # probs: numeric in [0,1]
  stopifnot(length(x) == length(w))
  ok <- is.finite(x) & is.finite(w) & w > 0
  x <- x[ok]; w <- w[ok]
  
  if (length(x) == 0) return(stats::setNames(rep(NA_real_, length(probs)), paste0("q", probs)))
  
  o <- order(x)
  x <- x[o]; w <- w[o]
  cw <- cumsum(w)
  W <- sum(w)
  if (W <= 0) return(stats::setNames(rep(NA_real_, length(probs)), paste0("q", probs)))
  
  # First x where cumulative weight >= p*W
  q <- vapply(probs, function(p) {
    x[which(cw >= p * W)[1]]
  }, numeric(1))
  
  stats::setNames(q, paste0("q", probs))
}
