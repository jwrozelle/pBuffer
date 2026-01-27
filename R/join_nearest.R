#' Probability-aware point-to-nearest-feature join from a density buffer
#'
#' @description
#' For each displaced point, applies a `densityBuffer` to obtain a probability
#' surface of plausible true locations. Each candidate location deterministically
#' maps to a nearest feature (by Euclidean distance in the projected CRS).
#' We then aggregate candidate weights by feature to obtain a probability
#' distribution over "true nearest feature".
#'
#' Also returns probability-weighted distance summaries and a lightweight
#' distance distribution (binned probability mass by distance).
#'
#' @param point_sf `sf` POINT layer of displaced points (one geometry per ID).
#' @param feature_sf `sf` POINT layer of candidate nearest features (e.g., facilities).
#' @param densityBuffer A density buffer object from `pb_densityBuffer()`.
#' @param templateRaster Optional `RasterLayer` template for `pb_apply_densityBuffer()`.
#' @param point_id Unique ID in `point_sf`.
#' @param feature_id Unique ID in `feature_sf`.
#' @param weight_name Name for probability mass column (default `"p"`).
#' @param min_weight Drop candidate cells with weights <= min_weight.
#' @param dist_quantiles Numeric vector of quantiles to report (default c(.05,.25,.5,.75,.95)).
#' @param dist_breaks Numeric vector of distance breaks (in CRS units, typically meters)
#'   used to build a binned probability mass distribution. If NULL, no distribution is returned.
#'
#' @returns A list with:
#'   - `probs_long`: (point_id, feature_id, p)
#'   - `map_nearest`: (point_id, feature_id, p)
#'   - `dist_summary`: (point_id, w_mean_dist, w_median_dist, q05.., q95..)
#'   - `dist_distribution`: (point_id, bin_lo, bin_hi, p) if dist_breaks supplied; else NULL
#'
#' @export
pb_pointProbJoin <- function(point_sf,
                             feature_sf,
                             densityBuffer,
                             templateRaster = NULL,
                             point_id = "DHSID",
                             feature_id = "SPAid",
                             weight_name = "p",
                             min_weight = 0,
                             dist_quantiles = c(.05, .25, .5, .75, .95),
                             dist_breaks = seq(0, 20000, by = 250)) {
  
  stopifnot(inherits(point_sf, "sf"), inherits(feature_sf, "sf"))
  if (!point_id %in% names(point_sf)) stop(sprintf("point_id '%s' not found in point_sf", point_id))
  if (!feature_id %in% names(feature_sf)) stop(sprintf("feature_id '%s' not found in feature_sf", feature_id))
  
  # Enforce one row per displaced point
  if (anyDuplicated(point_sf[[point_id]])) stop("point_id must be unique in point_sf (one row per displaced point).")
  if (anyDuplicated(feature_sf[[feature_id]])) stop("feature_id must be unique in feature_sf.")
  
  # CRS sanity: distances assume projected CRS with linear units
  if (sf::st_is_longlat(point_sf) || sf::st_is_longlat(feature_sf)) {
    stop("pb_pointProbJoin requires projected CRS (not lon/lat). Reproject first.")
  }
  
  # Apply density buffer to each point -> weight raster per point
  w_list <- pb_apply_densityBuffer(
    point_sf = point_sf,
    densityBuffer = densityBuffer,
    templateRaster = templateRaster
  )
  if (!is.list(w_list) || length(w_list) != nrow(point_sf)) {
    stop("pb_apply_densityBuffer() must return a list of RasterLayer, length == nrow(point_sf).")
  }
  
  # Prep minimal feature layer for nearest lookup
  feat_min <- feature_sf[, feature_id, drop = FALSE]
  
  probs_accum <- vector("list", length(w_list))
  dist_sum_accum <- vector("list", length(w_list))
  dist_dist_accum <- if (!is.null(dist_breaks)) vector("list", length(w_list)) else NULL
  
  for (i in seq_along(w_list)) {
    pid <- point_sf[[point_id]][i]
    w_r <- w_list[[i]]
    
    cand_sf <- pb__weights_raster_to_points(w_r, weight_name = weight_name, drop_zeros = TRUE, min_weight = min_weight)
    if (nrow(cand_sf) == 0) {
      probs_accum[[i]] <- data.frame(point_id = pid, feature_id = NA, p = NA_real_)
      dist_sum_accum[[i]] <- data.frame(point_id = pid, w_mean_dist = NA_real_, w_median_dist = NA_real_)
      if (!is.null(dist_dist_accum)) {
        dist_dist_accum[[i]] <- data.frame(point_id = pid, bin_lo = NA_real_, bin_hi = NA_real_, p = NA_real_)
      }
      next
    }
    
    # Nearest feature index for each candidate
    nn_idx <- sf::st_nearest_feature(cand_sf, feat_min)
    nn_feat <- feat_min[[feature_id]][nn_idx]
    
    # Distance from candidate to its nearest feature (units preserved)
    # st_distance returns units; convert to numeric in CRS units (typically meters)
    d <- sf::st_distance(sf::st_geometry(cand_sf), sf::st_geometry(feat_min)[nn_idx], by_element = TRUE)
    d_num <- as.numeric(d)
    
    w <- cand_sf[[weight_name]]
    
    # Aggregate probability by nearest feature
    agg <- stats::aggregate(w,
                            by = list(point = rep(pid, length(nn_feat)), feat = nn_feat),
                            FUN = sum)
    names(agg) <- c(point_id, feature_id, weight_name)
    probs_accum[[i]] <- agg
    
    # Distance summaries (weighted mean + weighted quantiles)
    w_mean <- sum(d_num * w) / sum(w)
    qs <- pb__weighted_quantile(d_num, w, probs = dist_quantiles)
    
    # Weighted median is quantile at 0.5 if present, else compute directly
    w_median <- if ("q0.5" %in% names(qs)) unname(qs[["q0.5"]]) else unname(pb__weighted_quantile(d_num, w, probs = 0.5)[[1]])
    
    dist_row <- data.frame(
      point_id = pid,
      w_mean_dist = w_mean,
      w_median_dist = w_median
    )
    # Append named quantiles as columns
    for (nm in names(qs)) dist_row[[nm]] <- unname(qs[[nm]])
    dist_sum_accum[[i]] <- dist_row
    
    # Optional: binned distance distribution (probability mass per bin)
    if (!is.null(dist_breaks)) {
      # Ensure breaks cover max distance; if not, extend last break
      if (max(d_num, na.rm = TRUE) > max(dist_breaks)) {
        dist_breaks2 <- c(dist_breaks, max(d_num, na.rm = TRUE))
      } else {
        dist_breaks2 <- dist_breaks
      }
      
      bin <- cut(d_num, breaks = dist_breaks2, include.lowest = TRUE, right = FALSE)
      # Sum probability mass by bin
      dist_df <- data.frame(bin = bin, w = w)
      dist_agg <- stats::aggregate(dist_df$w, by = list(bin = dist_df$bin), FUN = sum)
      names(dist_agg) <- c("bin", weight_name)
      
      # Parse bin interval into lo/hi for clean downstream plotting
      # Format like "[0,250)" (depends on locale); handle robustly
      bin_chr <- as.character(dist_agg$bin)
      # Extract numeric bounds
      lo <- as.numeric(sub("^\\[|\\(|,.*$", "", bin_chr))
      hi <- as.numeric(sub("^.*,(.*)\\)|\\]$", "\\1", bin_chr))
      
      dist_dist_accum[[i]] <- data.frame(
        point_id = pid,
        bin_lo = lo,
        bin_hi = hi,
        p = dist_agg[[weight_name]]
      )
    }
  }
  
  probs_long <- do.call(rbind, probs_accum)
  
  # MAP nearest feature per point
  map_nearest <- do.call(rbind, lapply(split(probs_long, probs_long[[point_id]]), function(d) {
    if (nrow(d) == 0 || all(!is.finite(d[[weight_name]]))) {
      return(data.frame(d[1, point_id, drop = FALSE], feature_id = NA, p = NA_real_))
    }
    d <- d[order(d[[weight_name]], decreasing = TRUE), , drop = FALSE]
    out <- d[1, c(point_id, feature_id, weight_name), drop = FALSE]
    names(out)[names(out) == weight_name] <- "p"
    out
  }))
  names(map_nearest)[names(map_nearest) == "p"] <- weight_name
  
  dist_summary <- do.call(rbind, dist_sum_accum)
  
  dist_distribution <- if (!is.null(dist_dist_accum)) do.call(rbind, dist_dist_accum) else NULL
  
  list(
    probs_long = probs_long,
    map_nearest = map_nearest,
    dist_summary = dist_summary,
    dist_distribution = dist_distribution
  )
}

