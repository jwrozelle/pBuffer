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
                             min_weight = 0) {
  
  stopifnot(inherits(point_sf, "sf"), inherits(poly_sf, "sf"))
  stopifnot(sf::st_geometry_type(point_sf, by_geometry = FALSE) %in% c("POINT", "MULTIPOINT"))
  
  # Ensure IDs exist
  if (!point_id %in% names(point_sf)) stop(sprintf("point_id '%s' not found in point_sf", point_id))
  if (!poly_id %in% names(poly_sf)) stop(sprintf("poly_id '%s' not found in poly_sf", poly_id))
  
  # Defensive: keep only needed columns
  point_ids <- point_sf[[point_id]]
  if (anyDuplicated(point_ids)) stop("point_id must be unique in point_sf (one row per displaced point).")
  
  # Apply density buffer to each point -> weight raster (one per point)
  # NOTE: this assumes pb_apply_densityBuffer can accept templateRaster;
  # if your existing signature differs, adapt the call here but keep output a RasterLayer.
  w_list <- pb_apply_densityBuffer(
    point_sf = point_sf,
    densityBuffer = densityBuffer,
    templateRaster = templateRaster
  )
  
  # Expect: list of RasterLayer, aligned to template, one per point in input order.
  if (!is.list(w_list) || length(w_list) != nrow(point_sf)) {
    stop("pb_apply_densityBuffer() must return a list of RasterLayer, length == nrow(point_sf).")
  }
  
  # Aggregate by polygon
  probs_accum <- vector("list", length(w_list))
  
  for (i in seq_along(w_list)) {
    w_r <- w_list[[i]]
    cand_sf <- pb__weights_raster_to_points(w_r, weight_name = weight_name, drop_zeros = TRUE, min_weight = min_weight)
    
    if (nrow(cand_sf) == 0) {
      probs_accum[[i]] <- data.frame(
        point_id = point_sf[[point_id]][i],
        poly_id  = NA,
        p        = NA_real_
      )
      next
    }
    
    # Join candidate points to polygons
    # Keep poly_id only to avoid inflating attributes
    poly_min <- poly_sf[, poly_id, drop = FALSE]
    
    joined <- sf::st_join(cand_sf, poly_min, join = sf::st_within, left = FALSE)
    if (nrow(joined) == 0) {
      probs_accum[[i]] <- data.frame(
        point_id = point_sf[[point_id]][i],
        poly_id  = NA,
        p        = NA_real_
      )
      next
    }
    
    # Aggregate weights by polygon
    df <- sf::st_drop_geometry(joined)
    # Standardize column names in output
    df[[point_id]] <- point_sf[[point_id]][i]
    
    # sum weights by polygon id
    agg <- stats::aggregate(df[[weight_name]],
                            by = list(point = df[[point_id]], poly = df[[poly_id]]),
                            FUN = sum)
    
    names(agg) <- c(point_id, poly_id, weight_name)
    probs_accum[[i]] <- agg
  }
  
  probs_long <- do.call(rbind, probs_accum)
  
  # MAP polygon for each point
  # (ties broken by first encountered; deterministic order of aggregate)
  poly_map <- do.call(rbind, lapply(split(probs_long, probs_long[[point_id]]), function(d) {
    if (nrow(d) == 0 || all(!is.finite(d[[weight_name]]))) {
      return(data.frame(d[1, point_id, drop = FALSE], poly_id = NA, p = NA_real_))
    }
    d <- d[order(d[[weight_name]], decreasing = TRUE), , drop = FALSE]
    out <- d[1, c(point_id, poly_id, weight_name), drop = FALSE]
    names(out)[names(out) == weight_name] <- "p"
    out
  }))
  
  # Normalize output column names a bit
  names(poly_map)[names(poly_map) == "p"] <- weight_name
  
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
                        adminID = NULL) {
  
  probs <- pb_polyProbsJoin(
    displaced.sf = displaced.sf,
    polygons.sf  = polygons.sf,
    displaced.id = displaced.id,
    polygons.id  = polygons.id,
    densityBuffer = densityBuffer,
    adminBound = adminBound,
    adminID = adminID,
    drop_zero = TRUE
  )
  
  # If a point has no polygon mass, return NA row for that ID
  if (nrow(probs) == 0) {
    out <- data.frame(
      id = displaced.sf[[displaced.id]],
      poly = NA_character_,
      p_hat = NA_real_
    )
    names(out) <- c(displaced.id, polygons.id, "p_hat")
    return(out)
  }
  
  # Pick MAP per displaced.id
  # (Base R split/apply to avoid extra deps)
  split_idx <- split(seq_len(nrow(probs)), probs[[displaced.id]])
  
  out_list <- lapply(names(split_idx), function(id_i) {
    rows <- split_idx[[id_i]]
    j <- rows[which.max(probs$p_hat[rows])]
    data.frame(
      id = probs[[displaced.id]][j],
      poly = probs[[polygons.id]][j],
      p_hat = probs$p_hat[j]
    )
  })
  
  out <- do.call(rbind, out_list)
  names(out) <- c(displaced.id, polygons.id, "p_hat")
  
  # Ensure every displaced point appears exactly once (even if none)
  missing_ids <- setdiff(displaced.sf[[displaced.id]], out[[displaced.id]])
  if (length(missing_ids) > 0) {
    out <- rbind(out, data.frame(
      tmp_id = missing_ids,
      tmp_poly = NA_character_,
      tmp_p = NA_real_
    ) |> `names<-`(c(displaced.id, polygons.id, "p_hat")))
  }
  
  # Preserve input order
  out <- out[match(displaced.sf[[displaced.id]], out[[displaced.id]]), , drop = FALSE]
  
  rownames(out) <- NULL
  out
}
