# ======================================================================
# join_polygon.R
#
# Probability-aware polygon joins:
#   - pb_polyProbsJoin(): for each displaced point, return ALL polygons with
#     their probability mass under the uncertainty kernel.
#   - pb_polyJoin(): return ONLY the most likely polygon (MAP) per displaced point.
#
# Interpretation:
#   For point i, candidate locations k have weights w_k (sum to 1). We assign
#   each candidate to a polygon, then aggregate weights by polygon id:
#     P(poly=j | i) = sum_{k in j} w_k
# ======================================================================

#' Probability-weighted polygon membership from a density buffer
#'
#' @description
#' For each displaced point, applies a `densityBuffer` (via `pb_apply_densityBuffer()`)
#' to generate a probability surface of plausible true locations. Candidate locations
#' (cell centers) are then joined to polygons and weights are aggregated to yield
#' a probability distribution over polygon membership.
#'
#' This is the core building block for:
#' - admin-area membership probabilities
#' - Thiessen polygon probabilities (probabilistic "nearest facility")
#'
#' @param point_sf `sf` POINT layer of displaced points (one geometry per ID).
#' @param poly_sf `sf` POLYGON/MULTIPOLYGON layer.
#' @param densityBuffer A density buffer object created by `pb_densityBuffer()`.
#' @param templateRaster Optional `RasterLayer` used as a template for rasterization
#'   (CRS/resolution/extent). If `NULL`, the function will try to use the
#'   densityBuffer's own raster metadata; if that is insufficient in your current
#'   implementation of `pb_apply_densityBuffer()`, pass a template explicitly.
#' @param point_id Name of unique ID field in `point_sf`.
#' @param poly_id Name of unique ID field in `poly_sf`.
#' @param weight_name Output column name for probability mass (default `"p"`).
#' @param min_weight Drop candidate cells with weights <= min_weight (default 0).
#' @param adminBound Optional admin polygons used to trim the kernel to the polygon
#'   containing each displaced point.
#' @param adminID ID column in `adminBound`.
#'
#' @returns A list with:
#'   - `probs_long`: data.frame with columns `{point_id, poly_id, p}`
#'   - `poly_map`:  data.frame with one row per point: `{point_id, poly_id, p}`
#'
#' @export
pb_polyProbsJoin <- function(point_sf,
                             poly_sf,
                             densityBuffer,
                             templateRaster = NULL,
                             point_id = "DHSID",
                             poly_id = "ID_2",
                             weight_name = "p",
                             min_weight = 0,
                             adminBound = NULL,
                             adminID = "ID_2") {
  
  stopifnot(inherits(point_sf, "sf"), inherits(poly_sf, "sf"))
  stopifnot(sf::st_geometry_type(point_sf, by_geometry = FALSE) %in% c("POINT", "MULTIPOINT"))
  
  # Ensure IDs exist
  if (!point_id %in% names(point_sf)) stop(sprintf("point_id '%s' not found in point_sf", point_id), call. = FALSE)
  if (!poly_id %in% names(poly_sf)) stop(sprintf("poly_id '%s' not found in poly_sf", poly_id), call. = FALSE)
  
  # CRS safety
  if (sf::st_crs(point_sf) != sf::st_crs(poly_sf)) {
    stop("point_sf and poly_sf must have the same CRS.", call. = FALSE)
  }
  
  # Enforce uniqueness (one row per displaced point)
  point_ids <- point_sf[[point_id]]
  if (anyDuplicated(point_ids)) {
    stop("point_id must be unique in point_sf (one row per displaced point).", call. = FALSE)
  }
  
  # Apply density buffer to each point -> weight raster(s)
  w_out <- pb_apply_densityBuffer(
    point_sf = point_sf,
    densityBuffer = densityBuffer,
    templateRaster = templateRaster,
    adminBound = adminBound,
    adminID = adminID
  )
  
  # Accept either:
  # - list of rasters (multi-point)
  # - single RasterLayer (single-point)
  if (inherits(w_out, "RasterLayer")) {
    if (nrow(point_sf) != 1) {
      stop("pb_apply_densityBuffer returned a single RasterLayer but point_sf has >1 row.", call. = FALSE)
    }
    w_list <- list(w_out)
  } else if (is.list(w_out)) {
    w_list <- w_out
  } else {
    stop("pb_apply_densityBuffer() must return a RasterLayer (nrow==1) or a list of RasterLayer/NULL.", call. = FALSE)
  }
  
  if (length(w_list) != nrow(point_sf)) {
    stop("pb_apply_densityBuffer() output length must equal nrow(point_sf).", call. = FALSE)
  }
  
  # Aggregate by polygon
  probs_accum <- vector("list", length(w_list))
  
  # Keep only needed polygon ID column to reduce attribute bloat
  poly_min <- poly_sf[, poly_id, drop = FALSE]
  
  for (i in seq_along(w_list)) {
    
    w_r <- w_list[[i]]
    
    # If kernel failed / trimmed to nothing, return NA row
    if (is.null(w_r)) {
      probs_accum[[i]] <- setNames(
        data.frame(
          point_sf[[point_id]][i],
          NA,
          NA_real_,
          stringsAsFactors = FALSE
        ),
        c(point_id, poly_id, weight_name)
      )
      next
    }
    
    cand_sf <- pb__weights_raster_to_points(
      w_r,
      weight_name = weight_name,
      drop_zeros  = TRUE,
      min_weight  = min_weight
    )
    
    if (!inherits(cand_sf, "sf") || nrow(cand_sf) == 0) {
      probs_accum[[i]] <- setNames(
        data.frame(
          point_sf[[point_id]][i],
          NA,
          NA_real_,
          stringsAsFactors = FALSE
        ),
        c(point_id, poly_id, weight_name)
      )
      next
    }
    
    # CRS safety (cand_sf should already be same CRS, but be explicit)
    if (sf::st_crs(cand_sf) != sf::st_crs(poly_min)) {
      sf::st_crs(cand_sf) <- sf::st_crs(poly_min)
    }
    
    joined <- sf::st_join(cand_sf, poly_min, join = sf::st_within, left = FALSE)
    
    if (nrow(joined) == 0) {
      probs_accum[[i]] <- setNames(
        data.frame(
          point_sf[[point_id]][i],
          NA,
          NA_real_,
          stringsAsFactors = FALSE
        ),
        c(point_id, poly_id, weight_name)
      )
      next
    }
    
    df <- sf::st_drop_geometry(joined)
    
    # Attach point_id explicitly
    df[[point_id]] <- point_sf[[point_id]][i]
    
    # Sum weights by polygon id
    agg <- stats::aggregate(
      x  = df[[weight_name]],
      by = list(point = df[[point_id]], poly = df[[poly_id]]),
      FUN = sum
    )
    
    names(agg) <- c(point_id, poly_id, weight_name)
    probs_accum[[i]] <- agg
  }
  
  probs_long <- do.call(rbind, probs_accum)
  
  # MAP polygon per point (ties broken deterministically by first after ordering)
  poly_map_list <- vector("list", length = nrow(point_sf))
  names(poly_map_list) <- as.character(point_sf[[point_id]])
  
  split_list <- split(probs_long, probs_long[[point_id]])
  
  for (pid in names(poly_map_list)) {
    
    d <- split_list[[pid]]
    
    if (is.null(d) || nrow(d) == 0 || all(!is.finite(d[[weight_name]]))) {
      poly_map_list[[pid]] <- setNames(
        data.frame(pid, NA, NA_real_, stringsAsFactors = FALSE),
        c(point_id, poly_id, weight_name)
      )
      next
    }
    
    d <- d[order(d[[weight_name]], decreasing = TRUE), , drop = FALSE]
    poly_map_list[[pid]] <- d[1, c(point_id, poly_id, weight_name), drop = FALSE]
  }
  
  poly_map <- do.call(rbind, poly_map_list)
  rownames(poly_map) <- NULL
  
  list(
    probs_long = probs_long,
    poly_map   = poly_map
  )
}


#' Most likely polygon under positional uncertainty (MAP)
#'
#' @description
#' Wrapper around `pb_polyProbsJoin()` that returns, for each displaced point,
#' the polygon with maximum estimated probability mass.
#'
#' @inheritParams pb_polyProbsJoin
#'
#' @return A data.frame with columns:
#'   - displaced.id
#'   - polygons.id (MAP polygon; NA if none)
#'   - p_hat (probability of MAP; NA if none)
#' @export
pb_polyJoin <- function(displaced.sf,
                        polygons.sf,
                        displaced.id = "DHSID",
                        polygons.id = "SPAID",
                        densityBuffer = NULL,
                        adminBound = NULL,
                        adminID = "ID_2",
                        templateRaster = NULL,
                        min_weight = 0) {
  
  if (is.null(densityBuffer)) {
    stop("densityBuffer must be provided.", call. = FALSE)
  }
  
  res <- pb_polyProbsJoin(
    point_sf       = displaced.sf,
    poly_sf        = polygons.sf,
    densityBuffer  = densityBuffer,
    templateRaster = templateRaster,
    point_id       = displaced.id,
    poly_id        = polygons.id,
    weight_name    = "p",
    min_weight     = min_weight,
    adminBound     = adminBound,
    adminID        = adminID
  )
  
  out <- res$poly_map
  names(out)[names(out) == "p"] <- "p_hat"
  
  # Ensure every displaced point appears exactly once, preserving input order
  all_ids <- displaced.sf[[displaced.id]]
  
  # Add missing rows if any (should be rare, but keep it robust)
  if (!all(all_ids %in% out[[displaced.id]])) {
    missing_ids <- setdiff(all_ids, out[[displaced.id]])
    if (length(missing_ids) > 0) {
      add <- setNames(
        data.frame(
          missing_ids,
          NA,
          NA_real_,
          stringsAsFactors = FALSE
        ),
        c(displaced.id, polygons.id, "p_hat")
      )
      out <- rbind(out, add)
    }
  }
  
  out <- out[match(all_ids, out[[displaced.id]]), , drop = FALSE]
  rownames(out) <- NULL
  out
}
