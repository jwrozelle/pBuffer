#' Link each point to its nearest feature
#'
#' Creates a deterministic nearest-neighbor linkage from a point layer (e.g., DHS clusters)
#' to a feature layer (e.g., facilities). For each point, the function identifies the
#' nearest feature (by planar distance in the current CRS) and returns a two-column
#' key pairing the point ID to the nearest feature ID.
#'
#' Notes:
#' - For meaningful distances, inputs should be in a **projected CRS**. If `hf_sf` is in a
#'   different CRS, it will be transformed to match `dhs_sf`.
#' - Nearest-feature lookup uses `sf::st_nearest_feature()`, which operates on geometries;
#'   `hf_sf` may be points, lines, or polygons.
#'
#' @param dhs_sf An `sf` object with POINT geometry representing origins (e.g., DHS clusters).
#' @param hf_sf An `sf` object representing candidate destination features (e.g., facilities).
#' @param dhsID Name of the unique ID column in `dhs_sf` (default `"DHSID"`).
#' @param hfID Name of the unique ID column in `hf_sf` (default `"SPAID"`).
#' @param out_hf_id Name of the output column containing the linked feature ID
#'   (default `"nearHFID"`).
#'
#' @return A data.frame with columns `dhsID` and `out_hf_id`. One row per input point.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' links <- pb_nearest_feature_link(dhs_sf, hf_sf, dhsID = "DHSID", hfID = "SPAID")
#' head(links)
#' }
pb_nearest_feature_link <- function(dhs_sf,
                                    hf_sf,
                                    dhsID = "DHSID",
                                    hfID = "SPAID",
                                    out_hf_id = "nearHFID") {
  
  if (!inherits(dhs_sf, "sf")) stop("dhs_sf must be an sf object.")
  if (!inherits(hf_sf, "sf")) stop("hf_sf must be an sf object.")
  if (!dhsID %in% names(dhs_sf)) stop("dhsID column not found in dhs_sf: ", dhsID)
  if (!hfID %in% names(hf_sf)) stop("hfID column not found in hf_sf: ", hfID)
  
  # Require CRS info
  crs_dhs <- sf::st_crs(dhs_sf)
  crs_hf  <- sf::st_crs(hf_sf)
  if (is.null(crs_dhs)) stop("dhs_sf has no CRS; set it before linking.")
  if (is.null(crs_hf))  stop("hf_sf has no CRS; set it before linking.")
  
  # Prefer POINT origins; fail loudly if not (your workflow assumes point clusters)
  if (!inherits(sf::st_geometry(dhs_sf), "sfc_POINT")) {
    stop("dhs_sf must have POINT geometry.")
  }
  
  # Transform hf to dhs CRS if needed
  if (crs_hf != crs_dhs) {
    hf_sf <- sf::st_transform(hf_sf, crs_dhs)
  }
  
  # Nearest feature index (row number in hf_sf for each dhs point)
  idx <- sf::st_nearest_feature(dhs_sf, hf_sf)
  
  # Build output keypair
  out <- data.frame(
    dhs_id = sf::st_drop_geometry(dhs_sf)[[dhsID]],
    hf_id  = sf::st_drop_geometry(hf_sf)[[hfID]][idx],
    stringsAsFactors = FALSE
  )
  
  names(out) <- c(dhsID, out_hf_id)
  out
}
