#' @keywords internal
pb__weights_raster_to_points <- function(w_raster,
                                         weight_name = "w",
                                         drop_zeros = TRUE,
                                         min_weight = 0) {
  # Convert a RasterLayer of weights to an sf POINT object with a weight column.
  # Uses cell centers as candidate locations.
  #
  # Returns: sf POINT with columns: <weight_name> plus geometry
  
  stopifnot(inherits(w_raster, "RasterLayer"))
  
  # Fast path: if all values are NA, return empty sf
  v <- raster::values(w_raster)
  if (is.null(v) || length(v) == 0 || all(!is.finite(v))) {
    empty <- sf::st_as_sf(
      data.frame(
        x = numeric(0),
        y = numeric(0),
        tmp_w = numeric(0)
      ),
      coords = c("x", "y"),
      crs = raster::crs(w_raster)
    )
    names(empty)[names(empty) == "tmp_w"] <- weight_name
    return(empty)
  }
  
  # rasterToPoints is stable and returns cell centers (x,y) + layer value.
  # NOTE: it will materialize points for all non-NA cells.
  pts_mat <- raster::rasterToPoints(w_raster, spatial = FALSE)
  
  if (!is.matrix(pts_mat) || nrow(pts_mat) == 0) {
    empty <- sf::st_as_sf(
      data.frame(
        x = numeric(0),
        y = numeric(0),
        tmp_w = numeric(0)
      ),
      coords = c("x", "y"),
      crs = raster::crs(w_raster)
    )
    names(empty)[names(empty) == "tmp_w"] <- weight_name
    return(empty)
  }
  
  colnames(pts_mat) <- c("x", "y", "tmp_w")
  
  ok <- is.finite(pts_mat[, "tmp_w"])
  if (drop_zeros) {
    ok <- ok & (pts_mat[, "tmp_w"] > min_weight)
  }
  pts_mat <- pts_mat[ok, , drop = FALSE]
  
  if (nrow(pts_mat) == 0) {
    empty <- sf::st_as_sf(
      data.frame(
        x = numeric(0),
        y = numeric(0),
        tmp_w = numeric(0)
      ),
      coords = c("x", "y"),
      crs = raster::crs(w_raster)
    )
    names(empty)[names(empty) == "tmp_w"] <- weight_name
    return(empty)
  }
  
  pts_df <- as.data.frame(pts_mat)
  names(pts_df)[names(pts_df) == "tmp_w"] <- weight_name
  
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
