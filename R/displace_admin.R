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


#' DHS-protocol displacement constrained to original administrative unit
#'
#' Applies DHS-style displacement (urban up to 2 km; rural up to 5 km with a small
#' fraction allowed up to 10 km) and **repeats the displacement draw** until every
#' displaced point falls within the same administrative polygon as its original
#' location.
#'
#' Internally, this function:
#' 1. assigns each original point to an admin polygon via `pb_assign_admin()`;
#' 2. displaces points once via `pb_displace_once_dhs()`;
#' 3. re-assigns displaced points to polygons;
#' 4. re-draws displacement only for points whose polygon changed (or became NA),
#'    repeating until all points match.
#'
#' @param point_sf An `sf` object with POINT geometry representing original locations.
#' @param poly_sf An `sf` object with POLYGON/MULTIPOLYGON geometry used to constrain
#'   displacement (e.g., admin2 boundaries).
#' @param point_id Name of the unique ID column in `point_sf` (default `"DHSID"`).
#' @param poly_id Name of the unique ID column in `poly_sf` (default `"ID_2"`).
#' @param verbose Logical; if TRUE, prints iteration progress messages (default TRUE).
#'
#' @return An `sf` object with displaced POINT geometry and original attributes, plus:
#'   - the original admin ID column (`poly_id`)
#'   - the displaced admin ID column (`paste0(poly_id, ".disp")`)
#'   - displacement diagnostic columns created by `pb_displace_once_dhs()` (e.g.,
#'     `originalX`, `originalY`, offsets, etc.)
#'
#' @export
#'
#' @examples
#' \dontrun{
#' displaced <- pb_displace_protocol_dhs(
#'   point_sf = dhs_sf,
#'   poly_sf = adm2_sf,
#'   point_id = "DHSID",
#'   poly_id = "ID_2"
#' )
#' }
pb_displace_protocol_dhs <- function(point_sf,
                                     poly_sf,
                                     point_id = "DHSID",
                                     poly_id = "ID_2",
                                     verbose = TRUE) {
  
  pb_check_sf(point_sf, "point_sf")
  pb_check_geom_type(point_sf, "POINT", "point_sf")
  pb_check_cols(point_sf, point_id, "point_sf")
  pb_check_unique_id(sf::st_drop_geometry(point_sf), point_id, "point_sf")
  pb_check_crs(point_sf, "point_sf")
  pb_check_projected(point_sf, "point_sf")
  
  pb_check_sf(poly_sf, "poly_sf")
  pb_check_cols(poly_sf, poly_id, "poly_sf")
  pb_check_unique_id(sf::st_drop_geometry(poly_sf), poly_id, "poly_sf")
  
  
  # Original admin IDs
  original_admins <- pb_assign_admin(point_sf, poly_sf, point_id = point_id, poly_id = poly_id)
  
  # Attach original admin IDs to points
  dataset.sf <- merge(point_sf, original_admins, by = point_id)
  
  # Storage for accepted displaced points (KEEP sf class!)
  final_displacement <- dataset.sf[0, ]
  
  i <- 0L
  mismatches <- nrow(dataset.sf)
  mismatched.sf <- NULL
  
  while (mismatches > 0) {
    i <- i + 1L
    
    # If we have mismatches from last iteration, reset geometry back to original coordinates
    if (!is.null(mismatched.sf) && nrow(mismatched.sf) > 0) {
      tmp <- sf::st_set_geometry(mismatched.sf, NULL)
      # pb_displace_once_dhs() expects geometry at original coords; reconstruct from originalX/Y
      dataset.sf <- sf::st_as_sf(tmp, coords = c("originalX", "originalY"), crs = crs_in)
    }
    
    if (verbose) {
      message(sprintf("Starting iteration %d with %d observations needing displacement", i, nrow(dataset.sf)))
    }
    
    # Displace once (uses DHS-style rules)
    displaced <- pb_displace_once_dhs(dataset.sf, useCRS = NULL)
    
    if (verbose) message(sprintf("Displaced %d clusters.", nrow(displaced)))
    
    displaced_admins <- pb_assign_admin(
      displaced,
      poly_sf,
      point_id = point_id,
      poly_id = poly_id,
      newsuffix = ".disp"
    )
    
    displaced <- merge(
      displaced,
      displaced_admins,
      by = point_id,
      all.x = TRUE,
      all.y = FALSE
    )
    
    # Identify mismatches:
    # mismatch if original polygon is non-missing and displaced polygon differs or becomes NA
    disp_col <- paste0(poly_id, ".disp")
    
    displaced$mismat <- 0L
    displaced$mismat <- ifelse(!is.na(displaced[[poly_id]]) &
                                 (is.na(displaced[[disp_col]]) | displaced[[disp_col]] != displaced[[poly_id]]),
                               1L, displaced$mismat)
    
    # Extract mismatched
    mismatched.sf <- displaced[displaced$mismat == 1L, ]
    
    # Drop the displaced polygon ID column so it can be re-created cleanly next loop
    if (nrow(mismatched.sf) > 0) {
      mismatched.sf[[disp_col]] <- NULL
      mismatched.sf$mismat <- NULL
    }
    
    mismatches <- nrow(mismatched.sf)
    
    # Keep the displacements that are in the same polygon
    matching_admins.sf <- displaced[displaced$mismat == 0L, ]
    matching_admins.sf$mismat <- NULL
    
    if (verbose) {
      message(sprintf("There are %d added to the final_displacement.", nrow(matching_admins.sf)))
    }
    
    if (nrow(matching_admins.sf) > 0) {
      final_displacement <- rbind(final_displacement, matching_admins.sf)
    }
    
    if (mismatches == 0L && verbose) {
      message(sprintf("Completed displacement of %d clusters.", nrow(final_displacement)))
    }
  }
  
  final_displacement
}