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

#' Polygon membership probabilities under positional uncertainty
#'
#' @description
#' Returns a discrete distribution of polygon membership probabilities for each
#' displaced point. Probability mass is computed by applying a `densityBuffer`
#' at each displaced point, assigning each candidate location to a polygon, and
#' summing candidate weights by polygon ID.
#'
#' If `densityBuffer` is NULL, this performs a conventional point-in-polygon join
#' and returns probability 1 for the polygon containing the point (or no rows if
#' none).
#'
#' @param displaced.sf `sf` POINT object of displaced locations (one per ID).
#' @param polygons.sf `sf` POLYGON/MULTIPOLYGON object to join against.
#' @param displaced.id Unique ID column in `displaced.sf`.
#' @param polygons.id Unique ID column in `polygons.sf`.
#' @param densityBuffer Optional pb_densityBuffer object (from `pb_densityBuffer()`).
#' @param adminBound Optional admin polygons used to trim kernels (passed to
#'   `pb_apply_densityBuffer()`).
#' @param adminID Polygon ID column in `adminBound` (required if `adminBound` provided).
#' @param drop_zero If TRUE, drop polygons with probability 0 (default TRUE).
#'
#' @return A data.frame with columns:
#'   - displaced.id
#'   - polygons.id
#'   - p_hat (probability mass; sums to ~1 per displaced.id if candidates cover space)
#' @export
pb_polyProbsJoin <- function(displaced.sf,
                             polygons.sf,
                             displaced.id = "DHSID",
                             polygons.id = "SPAID",
                             densityBuffer = NULL,
                             adminBound = NULL,
                             adminID = NULL,
                             drop_zero = TRUE) {
  
  # ---- validation ----
  pb_check_sf(displaced.sf, "displaced.sf")
  pb_check_geom_type(displaced.sf, "POINT", "displaced.sf")
  pb_check_cols(displaced.sf, displaced.id, "displaced.sf")
  pb_check_unique_id(sf::st_drop_geometry(displaced.sf), displaced.id, "displaced.sf")
  
  pb_check_sf(polygons.sf, "polygons.sf")
  pb_check_cols(polygons.sf, polygons.id, "polygons.sf")
  pb_check_unique_id(sf::st_drop_geometry(polygons.sf), polygons.id, "polygons.sf")
  
  if (sf::st_crs(displaced.sf) != sf::st_crs(polygons.sf)) {
    stop("CRS mismatch: displaced.sf and polygons.sf must have the same CRS.", call. = FALSE)
  }
  
  ids <- displaced.sf[[displaced.id]]
  
  # ---- conventional path (no densityBuffer) ----
  if (is.null(densityBuffer)) {
    
    joined <- suppressWarnings(sf::st_join(
      displaced.sf[, displaced.id, drop = FALSE],
      polygons.sf[, polygons.id, drop = FALSE],
      join = sf::st_within,
      left = FALSE
    ))
    
    if (nrow(joined) == 0) {
      out <- data.frame(
        displaced_id = character(0),
        polygon_id   = character(0),
        p_hat        = numeric(0)
      )
      names(out) <- c(displaced.id, polygons.id, "p_hat")
      return(out)
    }
    
    out <- sf::st_drop_geometry(joined)
    
    # Handle multi-hit (point on boundary or overlapping polygons):
    # keep first, but warn because conventional join can't disambiguate.
    dup <- out[[displaced.id]][duplicated(out[[displaced.id]])]
    if (length(dup) > 0) {
      warning("Some points intersect multiple polygons; keeping the first match per point.", call. = FALSE)
      out <- out[!duplicated(out[[displaced.id]]), , drop = FALSE]
    }
    
    out$p_hat <- 1
    return(as.data.frame(out))
  }
  
  # ---- probability-aware path ----
  res_list <- lapply(ids, function(id_i) {
    
    pt <- displaced.sf[displaced.sf[[displaced.id]] == id_i, , drop = FALSE]
    
    cand <- pb_apply_densityBuffer(
      displaced.sf = pt,
      densityBuffer = densityBuffer,
      adminBound = adminBound,
      adminID = adminID
    )
    
    pb_check_cols(cand, "layer", "candidate points")
    w <- pb_normalize_weights(cand[["layer"]], "candidate weights")
    
    # Assign candidates to polygons (within join)
    inter <- suppressWarnings(sf::st_join(
      cand[, "layer", drop = FALSE],
      polygons.sf[, polygons.id, drop = FALSE],
      join = sf::st_within,
      left = FALSE
    ))
    
    if (nrow(inter) == 0) {
      # No candidates fell in any polygon
      return(data.frame(
        id = id_i,
        poly = character(0),
        p_hat = numeric(0)
      ))
    }
    
    dt <- sf::st_drop_geometry(inter)
    
    # Aggregate weight mass by polygon id
    # Use base aggregate to avoid extra deps
    mass <- stats::aggregate(
      x  = dt[["layer"]],
      by = list(poly = dt[[polygons.id]]),
      FUN = sum, na.rm = TRUE
    )
    names(mass)[2] <- "p_hat"
    
    # Normalize defensively (pb_apply_densityBuffer should already do this)
    mass$p_hat <- pb_normalize_weights(mass$p_hat, "polygon probability mass")
    
    if (drop_zero) mass <- mass[mass$p_hat > 0, , drop = FALSE]
    
    data.frame(
      id = id_i,
      poly = as.character(mass$poly),
      p_hat = as.numeric(mass$p_hat)
    )
  })
  
  out <- do.call(rbind, res_list)
  
  # Standardize names
  if (nrow(out) == 0) {
    out <- data.frame(
      displaced_id = character(0),
      polygon_id   = character(0),
      p_hat        = numeric(0)
    )
    names(out) <- c(displaced.id, polygons.id, "p_hat")
    return(out)
  }
  
  names(out) <- c(displaced.id, polygons.id, "p_hat")
  out
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
