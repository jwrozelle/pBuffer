#' DHS-style single-draw coordinate displacement
#'
#' Displaces point coordinates according to DHS-style public release rules:
#' urban clusters are displaced uniformly up to 2 km; rural clusters are displaced
#' uniformly up to 5 km, with a small fraction (default 1%) displaced up to 10 km.
#'
#' This function performs **one** displacement draw per input point. If you want
#' a probability buffer (many candidate locations with weights), wrap this in a
#' Monte Carlo routine (e.g., replicate draws per point and aggregate).
#'
#' Important: This function assumes the input geometry is in a **projected CRS
#' with meter units**. If you pass lon/lat coordinates, the displacement will be
#' wrong because offsets are computed in meters.
#'
#' @param sfDataset An `sf` object with POINT geometry. Must contain a column
#'   indicating urban/rural status.
#' @param useCRS Optional CRS override. May be an EPSG integer, a `crs` object,
#'   or `NULL`. If `NULL` (default), the output CRS is taken from `sfDataset`.
#' @param urban_col Name of the urban/rural indicator column (default `"URBAN_RURA"`).
#'   Values are expected to be `"U"` and `"R"`.
#' @param p_rural_hi Proportion of rural points assigned the larger maximum
#'   displacement (default `0.01`).
#' @param rural_max Maximum displacement in meters for most rural points (default `5000`).
#' @param rural_max_hi Maximum displacement in meters for the `p_rural_hi` rural points
#'   (default `10000`).
#' @param urban_max Maximum displacement in meters for urban points (default `2000`).
#' @param uniform_area Logical. If `FALSE` (default), distances are drawn uniform in
#'   radius: `r ~ U(0, max)`. If `TRUE`, distances are drawn uniform over area within
#'   the circle: `r ~ sqrt(U(0,1)) * max`. DHS documentation is commonly interpreted as
#'   random displacement within a circle; analysts differ on whether that implies
#'   uniform radius vs uniform area—choose explicitly.
#'
#' @return An `sf` object with displaced POINT geometry. The output includes:
#'   `originalX`, `originalY`, `angle_rad`, `m_displaced`, `random_meters`,
#'   `xOffset`, `yOffset` as columns, plus all original non-geometry attributes.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # dhs_sf is an sf POINT layer in a projected CRS (meters) with URBAN_RURA column
#' displaced <- pb_displace_once_dhs(dhs_sf)
#' }
pb_displace_once_dhs <- function(sfDataset,
                                 useCRS = NULL,
                                 urban_col = "URBAN_RURA",
                                 p_rural_hi = 0.01,
                                 rural_max = 5000,
                                 rural_max_hi = 10000,
                                 urban_max = 2000,
                                 uniform_area = FALSE) {
  
  if (!inherits(sfDataset, "sf")) {
    stop("sfDataset must be an sf object.")
  }
  if (!inherits(sf::st_geometry(sfDataset), "sfc_POINT")) {
    stop("sfDataset must have POINT geometry.")
  }
  if (!urban_col %in% names(sfDataset)) {
    stop("sfDataset must contain the urban/rural column: ", urban_col)
  }
  
  crs_in <- sf::st_crs(sfDataset)
  if (is.null(crs_in) && is.null(useCRS)) {
    stop("Input CRS is missing. Provide useCRS or set CRS on sfDataset.")
  }
  
  # Output CRS: default to input CRS unless overridden
  crs_out <- if (is.null(useCRS)) crs_in else sf::st_crs(useCRS)
  
  n <- nrow(sfDataset)
  if (n == 0) return(sfDataset)
  
  # Extract coordinates once (faster + avoids repeated st_coordinates calls)
  xy <- sf::st_coordinates(sfDataset)
  
  # Random direction in radians
  angle_rad <- stats::runif(n, 0, 2 * pi)
  
  # Determine max displacement per point
  ur <- sfDataset[[urban_col]]
  is_rural <- ur == "R"
  is_urban <- ur == "U"
  
  # Default max displacement: urban_max for urban, rural_max for rural, NA otherwise
  m_displaced <- rep(NA_real_, n)
  m_displaced[is_urban] <- urban_max
  m_displaced[is_rural] <- rural_max
  
  # Assign p_rural_hi of rural points to rural_max_hi
  # DHS protocol: ~1% of rural clusters, rounded down.
  rural_idx <- which(is_rural)
  n_rural <- length(rural_idx)
  
  k <- 0L
  if (n_rural > 0 && is.finite(p_rural_hi) && p_rural_hi > 0 && rural_max_hi > rural_max) {
    k <- as.integer(floor(n_rural * p_rural_hi))
    if (k >= 1L) {
      hi_idx <- sample(rural_idx, size = k, replace = FALSE)
      m_displaced[hi_idx] <- rural_max_hi
    }
  }
  
  
  # Draw displacement distance (meters)
  u <- stats::runif(n, 0, 1)
  if (uniform_area) {
    random_meters <- sqrt(u) * m_displaced
  } else {
    random_meters <- u * m_displaced
  }
  
  # Offsets: x = cos(theta)*r, y = sin(theta)*r
  xOffset <- cos(angle_rad) * random_meters
  yOffset <- sin(angle_rad) * random_meters
  
  # New coordinates
  newX <- xy[, 1] + xOffset
  newY <- xy[, 2] + yOffset
  
  # Build output as attributes + new geometry
  out <- sf::st_drop_geometry(sfDataset)
  out$angle_rad <- angle_rad
  out$m_displaced <- m_displaced
  out$random_meters <- random_meters
  out$xOffset <- xOffset
  out$yOffset <- yOffset
  out$originalX <- xy[, 1]
  out$originalY <- xy[, 2]
  
  out$newX <- newX
  out$newY <- newY
  
  # Recreate sf with displaced geometry
  out_sf <- sf::st_as_sf(out, coords = c("newX", "newY"), crs = crs_out)
  
  out_sf
}
