# ---- probability rasterization + conversion helpers (internal) ----

#' Rasterize a centered displacement kernel around a point
#'
#' Internal helper used by join_* functions.
#'
#' @param densityBuffer list(weightedCircle, cellMeters, radiusMeters)
#' @param initialCoords numeric length-2 (x, y) in projected CRS units (meters)
#' @param inputCRS CRS object or WKT/PROJ string; passed to raster::crs()
#'
#' @return raster::RasterLayer with probability mass (sum ~ 1)
#' @keywords internal
rasterizeDisplacement <- function(densityBuffer, initialCoords, inputCRS = NULL) {
  
  if (is.null(densityBuffer) || !is.list(densityBuffer)) {
    stop("densityBuffer must be a list.", call. = FALSE)
  }
  req <- c("weightedCircle", "cellMeters", "radiusMeters")
  if (!all(req %in% names(densityBuffer))) {
    stop("densityBuffer must contain: weightedCircle, cellMeters, radiusMeters.", call. = FALSE)
  }
  
  m <- densityBuffer[["weightedCircle"]]
  cell <- densityBuffer[["cellMeters"]]
  
  if (!is.matrix(m)) stop("densityBuffer$weightedCircle must be a matrix.", call. = FALSE)
  if (!is.numeric(cell) || length(cell) != 1 || is.na(cell) || cell <= 0) {
    stop("densityBuffer$cellMeters must be a positive scalar.", call. = FALSE)
  }
  if (!is.numeric(initialCoords) || length(initialCoords) < 2 || anyNA(initialCoords[1:2])) {
    stop("initialCoords must be numeric length-2 c(x, y).", call. = FALSE)
  }
  
  nrows <- nrow(m)
  ncols <- ncol(m)
  
  x <- initialCoords[1]
  y <- initialCoords[2]
  
  half_x <- (ncols * cell) / 2
  half_y <- (nrows * cell) / 2
  
  r <- raster::raster(
    nrows = nrows,
    ncols = ncols,
    xmn = x - half_x,
    xmx = x + half_x,
    ymn = y - half_y,
    ymx = y + half_y
  )
  
  if (!is.null(inputCRS)) {
    raster::crs(r) <- inputCRS
  }
  
  # Fill raster with the kernel values
  # (t(m) mapping matches common raster<-matrix orientation expectations)
  r[] <- as.vector(t(m))
  
  # Normalize defensively
  s <- sum(r[], na.rm = TRUE)
  if (is.finite(s) && s > 0) r[] <- r[] / s
  
  r
}

#' Convert probability raster to sf points with weight column "layer"
#'
#' Internal helper used by join_* functions.
#'
#' @param probRaster raster::RasterLayer
#' @param weightsCol name for weights column (default "layer")
#'
#' @return sf POINT with weights column
#' @keywords internal
st_rasterAsPoint <- function(probRaster, weightsCol = "layer") {
  
  if (!inherits(probRaster, "RasterLayer")) {
    stop("probRaster must be a raster::RasterLayer.", call. = FALSE)
  }
  
  sp_pts <- raster::rasterToPoints(probRaster, spatial = TRUE)
  if (is.null(sp_pts) || nrow(sp_pts) == 0) {
    # Return empty sf with correct columns
    out <- sf::st_sf(data.frame(!!weightsCol := numeric(0)), geometry = sf::st_sfc(), crs = sf::st_crs(NA))
    return(out)
  }
  
  sf_pts <- sf::st_as_sf(sp_pts)
  
  # Ensure the weight column is named consistently
  dat_cols <- setdiff(names(sf_pts), attr(sf_pts, "sf_column"))
  if (!(weightsCol %in% names(sf_pts))) {
    # rasterToPoints usually yields a single data column (often named "layer")
    if (length(dat_cols) >= 1) {
      names(sf_pts)[names(sf_pts) == dat_cols[1]] <- weightsCol
    } else {
      stop("Could not identify raster value column to rename to weightsCol.", call. = FALSE)
    }
  }
  
  sf_pts
}
