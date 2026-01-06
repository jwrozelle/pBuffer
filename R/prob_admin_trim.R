#' Trim/renormalize a probability raster to the admin polygon of a point
#'
#' Internal helper: given a probability raster centered on a point, mask it to the
#' admin polygon that contains that point and renormalize.
#'
#' @param probRaster raster::RasterLayer probability surface (values sum ~1 before trim)
#' @param point_sf single-row sf POINT (the community/cluster)
#' @param adminBound sf polygons
#' @param adminID admin polygon id column in adminBound
#'
#' @return raster::RasterLayer with values summing to 1.
#' @keywords internal
pb_trim_probRaster_to_point_admin <- function(probRaster,
                                              point_sf,
                                              adminBound,
                                              adminID = "ID_2") {
  
  pb_check_sf(point_sf, "point_sf")
  pb_check_geom_type(point_sf, "POINT", "point_sf")
  if (nrow(point_sf) != 1) stop("point_sf must have exactly 1 row.", call. = FALSE)
  
  pb_check_sf(adminBound, "adminBound")
  pb_check_cols(adminBound, adminID, "adminBound")
  
  # Identify containing polygon
  hit <- suppressWarnings(sf::st_join(point_sf, adminBound, join = sf::st_within, left = TRUE))
  pid <- sf::st_drop_geometry(hit)[[adminID]][1]
  
  if (is.na(pid)) {
    stop("Point does not fall within any adminBound polygon; cannot trim.", call. = FALSE)
  }
  
  poly <- adminBound[adminBound[[adminID]] == pid, ]
  if (nrow(poly) != 1) stop("Expected exactly 1 matching admin polygon.", call. = FALSE)
  
  trimProbBuff(probRaster = probRaster, adminBound = poly)
}
