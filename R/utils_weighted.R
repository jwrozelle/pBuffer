# ---- internal probability-weighted aggregation helpers (not exported) ----

#' Probability-weighted aggregation of polygon attributes given a probability surface
#'
#' @description
#' Internal helper that computes probability-weighted summaries of polygon-level attributes
#' (e.g., facility readiness scores or other metrics) for a single displaced location.
#'
#' The function intersects a set of probability-weighted geometries (typically points or
#' small polygons representing probability mass) with a polygon layer, aggregates the
#' probability mass to the polygon level, normalizes weights to sum to one, and then
#' returns weighted means of the requested polygon attributes.
#'
#' @details
#' **Intended use-case:** You have a displaced community location, you have already
#' generated a probability surface around that location (e.g., via `pb_Density()` /
#' `pb_integratedDensity()`), and you want a single best estimate of polygon attributes
#' under that uncertainty (plus the polygon-level weights used to construct it).
#'
#' **Weighting logic:** Let \eqn{w_i} be the total probability mass intersecting polygon
#' \eqn{i}. The function normalizes these so that \eqn{\sum_i w_i = 1}, then for each
#' metric \eqn{m}, returns \eqn{\sum_i w_i m_i}. Normalization is deliberate: intersections
#' can drop mass when buffers extend beyond the polygon set, when polygons do not fully
#' cover the support of the probability surface, or when geometries are invalid.
#'
#' **Failure modes:** The function errors if there are no intersections between the
#' probability layer and polygons (i.e., no basis to allocate mass).
#'
#' @param bufferPoints.sf An `sf` object representing the probability surface for a single
#' displaced observation. Must contain a numeric column given by `weightVar` (default:
#' `"layer"`) that represents probability mass or relative weights.
#'
#' @param polygons.sf An `sf` polygon layer containing polygon identifiers (`poly_id`) and
#' polygon-level attributes to be weighted (given by `metrics`).
#'
#' @param poly_id Name of the polygon identifier column in `polygons.sf`. Defaults to
#' `"SPAID"`. The identifier is used only to aggregate weights and join polygon attributes.
#'
#' @param weightVar Name of the weight/probability column in `bufferPoints.sf`. Defaults
#' to `"layer"`.
#'
#' @param metrics Character vector of column names in `polygons.sf` giving the polygon-level
#' attributes to aggregate (e.g., readiness scores).
#'
#' @return A list with two elements:
#' \describe{
#'   \item{weightedMetrics.df}{A one-row `data.frame` with columns named by `metrics`,
#'   containing probability-weighted means of each metric.}
#'   \item{hfWeights}{A `data.frame` with columns `poly_id` (as named in `polygons.sf`) and
#'   `weight`, giving the normalized weight assigned to each polygon.}
#' }
#'
#' @keywords internal
pb_weighted_metrics <- function(
    bufferPoints.sf,
    polygons.sf,
    poly_id = "SPAID",
    weightVar = "layer",
    metrics
) {
  
  pb_check_sf(bufferPoints.sf, "bufferPoints.sf")
  pb_check_sf(polygons.sf, "polygons.sf")
  
  pb_check_cols(bufferPoints.sf, weightVar, "bufferPoints.sf")
  pb_check_cols(polygons.sf, poly_id, "polygons.sf")
  pb_check_cols(polygons.sf, metrics, "polygons.sf")
  
  
  if (length(metrics) == 0) {
    stop("metrics must contain at least one column name.", call. = FALSE)
  }
  metrics <- unique(metrics)
  
  
  # Intersect; keep warning suppression local
  intersection <- suppressWarnings(sf::st_intersection(
    polygons.sf[, c(poly_id, metrics)],
    bufferPoints.sf[, c(weightVar)]
  ))
  
  if (nrow(intersection) == 0) {
    stop("No intersections between bufferPoints.sf and polygons.sf; cannot compute weighted metrics.", call. = FALSE)
  }
  
  # Aggregate probability mass to polygon id
  int_result <- intersection |>
    sf::st_drop_geometry() |>
    dplyr::group_by(.data[[poly_id]]) |>
    dplyr::summarise(weight = sum(.data[[weightVar]], na.rm = TRUE), .groups = "drop")
  
  # Normalize weights to sum to 1 (robust to partial overlap / dropped pieces)
  if (any(int_result$weight < 0, na.rm = TRUE)) {
    stop("intersection weights contain negative values; cannot normalize.", call. = FALSE)
  }
  int_result$weight <- pb_normalize_weights(int_result$weight, "intersection weights")
  
  includedScores <- merge(
    int_result,
    sf::st_drop_geometry(polygons.sf[, c(poly_id, metrics)]),
    by = poly_id,
    all.x = TRUE,
    all.y = FALSE
  )
  
  # Weighted means; remove NAs at the metric-level
  est <- vapply(
    metrics,
    function(m) stats::weighted.mean(includedScores[[m]], w = includedScores$weight, na.rm = TRUE),
    numeric(1)
  )
  
  weightedMetrics.df <- as.data.frame(t(est), stringsAsFactors = FALSE)
  colnames(weightedMetrics.df) <- metrics
  
  list(weightedMetrics.df = weightedMetrics.df, hfWeights = int_result)
}