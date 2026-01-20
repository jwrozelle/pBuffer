#' Assign administrative polygon IDs to points
#'
#' Joins point locations to a polygon layer using a within-join, returning a
#' lightweight data.frame of point IDs and polygon IDs (geometry removed).
#'
#' This is a helper for DHS-style displacement protocols that require displaced
#' points to remain within the same administrative unit as the original point.
#'
#' @param point_sf An `sf` object with POINT geometry.
#' @param poly_sf An `sf` object with POLYGON/MULTIPOLYGON geometry.
#' @param point_id Name of the unique ID column in `point_sf` (default `"DHSID"`).
#' @param poly_id Name of the unique ID column in `poly_sf` (default `"ID_2"`).
#' @param newsuffix Optional suffix to append to `poly_id` in the output column name
#'   (default `""`). Useful for creating parallel columns like `"ID_2.disp"`.
#'
#' @return A data.frame with two columns: `point_id` and `paste0(poly_id, newsuffix)`.
#'   Geometry is dropped.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' admins <- pb_assign_admin(dhs_sf, adm2_sf, point_id = "DHSID", poly_id = "ID_2")
#' }
pb_assign_admin <- function(point_sf,
                            poly_sf,
                            point_id = "DHSID",
                            poly_id = "ID_2",
                            newsuffix = "") {
  
  pb_check_sf(point_sf, "point_sf")
  pb_check_geom_type(point_sf, "POINT", "point_sf")
  pb_check_cols(point_sf, point_id, "point_sf")
  pb_check_unique_id(sf::st_drop_geometry(point_sf), point_id, "point_sf")
  
  pb_check_sf(poly_sf, "poly_sf")
  pb_check_cols(poly_sf, poly_id, "poly_sf")
  pb_check_unique_id(sf::st_drop_geometry(poly_sf), poly_id, "poly_sf")
  
  
  joined <- sf::st_join(
    point_sf,
    poly_sf,
    join = sf::st_within,
    suffix = c("", newsuffix),
    left = TRUE
  )
  
  out_cols <- c(point_id, paste0(poly_id, newsuffix))
  if (!all(out_cols %in% names(joined))) {
    stop("Expected columns not found after join: ", paste(out_cols, collapse = ", "))
  }
  
  sf::st_set_geometry(joined[, out_cols, drop = FALSE], NULL)
}


#' DHS displacement constrained to administrative boundaries (protocol)
#'
#' @description
#' Implements a DHS-style displacement protocol where displaced points must remain
#' within the same administrative polygon as their original locations.
#'
#' The protocol is:
#' 1) Assign each original point to an admin polygon (within-join).
#' 2) Displace points once according to DHS rules.
#' 3) Re-assign displaced points to polygons.
#' 4) Any points that changed polygons are re-displaced (re-drawn) until all points
#'    remain within their original polygon.
#'
#' @param point_sf An `sf` object with POINT geometry.
#' @param poly_sf An `sf` object with POLYGON/MULTIPOLYGON geometry.
#' @param point_id Name of unique ID column in `point_sf` (default `"DHSID"`).
#' @param poly_id Name of polygon ID column in `poly_sf` (default `"ID_2"`).
#' @param useCRS Optional CRS override passed through to `pb_displace_once_dhs()`.
#' @param urban_col Name of urban/rural indicator column (default `"URBAN_RURA"`).
#' @param p_rural_hi Proportion of rural points displaced to `rural_max_hi`.
#' @param rural_max Max rural displacement in meters.
#' @param rural_max_hi Max rural displacement in meters for the tail fraction.
#' @param urban_max Max urban displacement in meters.
#' @param uniform_area If TRUE, sample displacement radius uniformly over area.
#' @param verbose If TRUE, prints progress messages.
#'
#' @return An `sf` object of displaced points, with the original admin ID retained.
#' @export
pb_displace_protocol_dhs <- function(point_sf,
                                     poly_sf,
                                     point_id = "DHSID",
                                     poly_id = "ID_2",
                                     useCRS = NULL,
                                     urban_col = "URBAN_RURA",
                                     p_rural_hi = 0.01,
                                     rural_max = 5000,
                                     rural_max_hi = 10000,
                                     urban_max = 2000,
                                     uniform_area = FALSE,
                                     verbose = TRUE) {
  
  # ---- validation ----
  pb_check_sf(point_sf, "point_sf")
  pb_check_geom_type(point_sf, "POINT", "point_sf")
  pb_check_cols(point_sf, point_id, "point_sf")
  pb_check_unique_id(sf::st_drop_geometry(point_sf), point_id, "point_sf")
  pb_check_crs(point_sf, "point_sf")
  pb_check_projected(point_sf, "point_sf")
  
  pb_check_sf(poly_sf, "poly_sf")
  pb_check_cols(poly_sf, poly_id, "poly_sf")
  pb_check_unique_id(sf::st_drop_geometry(poly_sf), poly_id, "poly_sf")
  
  # ---- original admin IDs ----
  original_admins <- pb_assign_admin(
    point_sf,
    poly_sf,
    point_id = point_id,
    poly_id = poly_id
  )
  
  # Attach original admin IDs to points (preserve sf)
  dataset.sf <- merge(point_sf, original_admins, by = point_id)
  
  # Storage for accepted displaced points (keep sf class)
  final_displacement <- dataset.sf[0, ]
  
  # Everything starts mismatched until proven otherwise
  mismatched.sf <- dataset.sf
  mismatches <- nrow(mismatched.sf)
  iter <- 0L
  
  # ---- protocol loop ----
  while (mismatches > 0L) {
    iter <- iter + 1L
    if (verbose) message(sprintf("Protocol iteration %d; redrawing %d points.", iter, mismatches))
    
    # Draw a displacement for currently mismatched points
    displaced <- pb_displace_once_dhs(
      sfDataset    = mismatched.sf,
      useCRS       = useCRS,          # FIX: respect user argument
      urban_col    = urban_col,
      p_rural_hi   = p_rural_hi,
      rural_max    = rural_max,
      rural_max_hi = rural_max_hi,
      urban_max    = urban_max,
      uniform_area = uniform_area     # FIX: respect user argument
    )
    
    # Re-assign displaced points to polygons
    displaced_admins <- pb_assign_admin(
      displaced,
      poly_sf,
      point_id = point_id,
      poly_id  = poly_id,
      newsuffix = ".disp"
    )
    
    # Merge the displaced admin IDs onto displaced points
    displaced <- merge(displaced, displaced_admins, by = point_id)
    
    # Determine which points stayed in the same admin unit
    orig_col <- poly_id
    disp_col <- paste0(poly_id, ".disp")
    
    # If any displaced points fail assignment (e.g., topological edge cases),
    # treat them as mismatches so they get redrawn.
    same_admin <- !is.na(displaced[[disp_col]]) & displaced[[orig_col]] == displaced[[disp_col]]
    
    matching_admins.sf <- displaced[same_admin, , drop = FALSE]
    mismatched.sf      <- displaced[!same_admin, , drop = FALSE]
    
    mismatches <- nrow(mismatched.sf)
    
    # Accumulate the accepted points
    if (nrow(matching_admins.sf) > 0) {
      final_displacement <- rbind(final_displacement, matching_admins.sf)
    }
    
    if (mismatches == 0L && verbose) {
      message(sprintf("Completed displacement of %d clusters.", nrow(final_displacement)))
    }
  }
  
  # ---- return ----
  # FIX: remove invalid rowname assignment; if rownames are desired, base them on the output ID.
  rownames(final_displacement) <- paste0("disp_", final_displacement[[point_id]])
  
  final_displacement
}
