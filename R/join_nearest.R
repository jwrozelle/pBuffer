#' Probability-aware nearest-point distance join under displacement uncertainty
#'
#' @description
#' For each displaced point, this function builds the set of plausible true
#' locations implied by \code{densityBuffer}, treats those as weighted candidate
#' \code{sf} points, and then assigns each candidate location to its nearest
#' point or multipoint feature.
#'
#' The function then summarizes that nearest-feature assignment in three ways:
#' \enumerate{
#'   \item an optional long table giving the displaced-point ID, nearest-feature
#'   ID, nearest distance in meters, and the probability mass for that
#'   feature-distance pair;
#'   \item a one-row summary per displaced point giving the probability-weighted
#'   median nearest distance; and
#'   \item a one-row summary per displaced point giving the most likely nearest
#'   distance after optional rounding or bucketing with \code{round()}.
#' }
#'
#' A few implementation choices are deliberate here.
#'
#' First, I am not routing this through an intermediate raster join anymore.
#' I reconstruct the candidate support directly as weighted \code{sf} points,
#' which is both clearer and usually lighter on memory.
#'
#' Second, the long table can get large if you keep exact distances. So I let
#' you control whether it is returned at all, whether it is also streamed to
#' disk chunk-by-chunk, and whether distances are rounded before aggregation.
#' If you want a stable, interpretable mode and a smaller long table, it is a
#' good idea to set \code{distance_digits} to something like \code{-1} or
#' \code{-2} so the mode is calculated on 10 m or 100 m buckets rather than on
#' essentially exact distances.
#'
#' @returns A list with three elements:
#' \describe{
#'   \item{nearest_distance_probabilities}{A data frame with four columns:
#'   displaced-point ID, nearest-feature ID, nearest distance in meters, and the
#'   probability for that feature-distance pair. If \code{return_long = FALSE},
#'   this element is returned as \code{NULL}. If \code{long_file} is supplied,
#'   this same table is also written to disk in chunks.}
#'   \item{median_distance}{A data frame with two columns: displaced-point ID
#'   and the probability-weighted median nearest distance in meters.}
#'   \item{max_likely_distance}{A data frame with three columns:
#'   displaced-point ID, the maximum-likelihood nearest distance in meters after
#'   optional rounding via \code{distance_digits}, and the probability mass of
#'   that distance bucket.}
#' }
#'
#' @param point_sf An \code{sf} object with POINT geometries for displaced
#' locations.
#' @param feature_sf An \code{sf} object with POINT or MULTIPOINT geometries for
#' the candidate nearest features.
#' @param densityBuffer A density-buffer object created by
#' \code{pb_densityBuffer()}, \code{pb_integratedDensity()},
#' \code{pb_Density()}, \code{pb_mcDensity()}, or
#' \code{pb_mcDensity_dhs()}.
#' @param templateRaster Kept only for backward compatibility. It is not used by
#' this implementation because the nearest-feature join is carried out directly
#' on candidate \code{sf} points.
#' @param point_id Column name giving the unique displaced-point ID.
#' @param feature_id Column name giving the unique nearest-feature ID.
#' @param weight_name Retained only for backward compatibility with older calls
#' to \code{pb_pointProbJoin()}. It is ignored by this implementation, which
#' always returns a probability column named \code{probability}.
#' @param adminBound Optional administrative boundary used to trim plausible true
#' locations before assigning nearest features. If multiple admin polygons are
#' supplied, the displaced point is first matched to its containing admin
#' polygon and only that polygon is used for trimming.
#' @param adminID Column name giving the admin ID in \code{adminBound}. This is
#' only required when \code{adminBound} has more than one row.
#' @param min_weight Drop candidate true-location support points with weights
#' less than or equal to this threshold before assigning nearest features, then
#' renormalize the remaining weights so the resulting probabilities still sum to
#' 1.
#' @param dist_quantiles Retained only for backward compatibility with older
#' calls. It is ignored by this implementation.
#' @param dist_breaks Retained only for backward compatibility with older calls.
#' It is ignored by this implementation.
#' @param distance_digits Optional single numeric value passed to
#' \code{round(x, digits = distance_digits)} when aggregating nearest distances
#' for the long output and for the maximum-likelihood distance. Use
#' \code{NULL} to keep exact distances, positive values for decimal places,
#' \code{0} for whole meters, and negative values for coarser buckets such as
#' 10 m or 100 m.
#' @param return_long Logical. If \code{TRUE}, return the long
#' feature-distance-probability table in memory. If \code{FALSE}, skip storing
#' it in memory.
#' @param long_file Optional file path for writing the long table to disk in
#' chunked CSV form. This is useful when the long output is too large to keep in
#' memory. If supplied together with \code{return_long = TRUE}, the long table
#' is both written to disk and returned in memory.
#' @param overwrite_long_file Logical. If \code{TRUE} and \code{long_file}
#' already exists, overwrite it. If \code{FALSE}, error instead.
#' @param chunk_size Optional number of displaced points to process per chunk.
#' If \code{NULL}, the function chooses a default that is conservative when the
#' long table is requested and otherwise processes everything in one chunk.
#' @param n.cores Number of workers to use when parallel execution is requested.
#' @param parallel Logical. If \code{TRUE}, evaluate displaced points in
#' parallel with \code{future.apply}. If \code{FALSE}, evaluation is serial
#' unless \code{n.cores > 1}.
#' @param future.seed Seed handling passed to
#' \code{future.apply::future_lapply()}.
#' @param future.scheduling Scheduling argument passed to
#' \code{future.apply::future_lapply()}.
#'
#' @author J.W. Rozelle
#'
#' @export pb_pointProbJoin
#' @examples
#' # coming soon!
pb_pointProbJoin <- function(point_sf,
                             feature_sf,
                             densityBuffer,
                             templateRaster = NULL,
                             point_id = "DHSID",
                             feature_id = "SPAid",
                             weight_name = "p",
                             min_weight = 0,
                             dist_quantiles = c(0.05, 0.25, 0.5, 0.75, 0.95),
                             dist_breaks = seq(0, 20000, by = 250),
                             adminBound = NULL,
                             adminID = NULL,
                             distance_digits = NULL,
                             return_long = TRUE,
                             long_file = NULL,
                             overwrite_long_file = FALSE,
                             chunk_size = NULL,
                             n.cores = 1,
                             parallel = FALSE,
                             future.seed = TRUE,
                             future.scheduling = 1) {
  
  # Future-me: the displaced layer has to be clean POINT geometry with one row
  # per ID, otherwise the probability support I build below is not well-defined.
  pb_check_sf(point_sf, "point_sf")
  pb_check_geom_type(point_sf, "POINT", "point_sf")
  pb_check_crs(point_sf, "point_sf")
  pb_check_projected(point_sf, "point_sf")
  pb_check_cols(point_sf, point_id, "point_sf")
  pb_check_unique_id(sf::st_drop_geometry(point_sf), point_id, "point_sf")
  
  # Future-me: the nearest-feature layer can be POINT or MULTIPOINT, but it
  # still needs a clean CRS and a unique feature ID so the output is interpretable.
  pb_check_sf(feature_sf, "feature_sf")
  pb_check_crs(feature_sf, "feature_sf")
  pb_check_cols(feature_sf, feature_id, "feature_sf")
  pb_check_unique_id(sf::st_drop_geometry(feature_sf), feature_id, "feature_sf")
  
  feature_geom <- unique(as.character(sf::st_geometry_type(feature_sf)))
  if (!all(feature_geom %in% c("POINT", "MULTIPOINT"))) {
    stop(
      "feature_sf must have POINT or MULTIPOINT geometries.",
      call. = FALSE
    )
  }
  
  if (nrow(feature_sf) == 0) {
    stop("feature_sf must contain at least one feature.", call. = FALSE)
  }
  
  pb_check_pos_scalar(n.cores, "n.cores")
  
  if (!is.numeric(min_weight) || length(min_weight) != 1 || is.na(min_weight) || min_weight < 0) {
    stop("min_weight must be a single non-negative number.", call. = FALSE)
  }
  
  if (!is.null(distance_digits)) {
    if (!is.numeric(distance_digits) || length(distance_digits) != 1 || is.na(distance_digits) || !is.finite(distance_digits)) {
      stop("distance_digits must be NULL or a single finite numeric value.", call. = FALSE)
    }
  }
  
  if (!is.logical(return_long) || length(return_long) != 1 || is.na(return_long)) {
    stop("return_long must be TRUE or FALSE.", call. = FALSE)
  }
  
  if (!is.null(long_file)) {
    if (!is.character(long_file) || length(long_file) != 1 || is.na(long_file) || nchar(long_file) == 0) {
      stop("long_file must be NULL or a single non-empty file path.", call. = FALSE)
    }
  }
  
  if (!is.logical(overwrite_long_file) || length(overwrite_long_file) != 1 || is.na(overwrite_long_file)) {
    stop("overwrite_long_file must be TRUE or FALSE.", call. = FALSE)
  }
  
  if (!is.null(chunk_size)) {
    if (!is.numeric(chunk_size) || length(chunk_size) != 1 || is.na(chunk_size) || chunk_size <= 0) {
      stop("chunk_size must be NULL or a single positive number.", call. = FALSE)
    }
    chunk_size <- as.integer(ceiling(chunk_size))
  }
  
  # Future-me: templateRaster stays in the signature only so older calling code
  # does not break. I am intentionally not using it in this implementation.
  
  # Future-me: weight_name, dist_quantiles, and dist_breaks are also still in
  # the signature so old calling code does not fail on unused arguments. The new
  # return object does not use them, because this rewrite is focused on the
  # nearest-feature probability table plus median and modal distances only.
  
  if (is.null(densityBuffer)) {
    stop(
      "densityBuffer must be provided (for example from pb_densityBuffer(), ",
      "pb_integratedDensity(), pb_Density(), pb_mcDensity(), or pb_mcDensity_dhs()).",
      call. = FALSE
    )
  }
  
  # Future-me: I want this to work with both the newer pb_densityBuffer objects
  # and the older list-style density buffers that still exist elsewhere here.
  density_is_new <- inherits(densityBuffer, "pb_densityBuffer")
  density_is_legacy <- is.list(densityBuffer) &&
    all(c("weightedCircle", "cellMeters", "radiusMeters") %in% names(densityBuffer))
  
  if (!density_is_new && !density_is_legacy) {
    stop(
      "densityBuffer must be either a pb_densityBuffer object or a legacy list ",
      "with weightedCircle, cellMeters, and radiusMeters.",
      call. = FALSE
    )
  }
  
  # Future-me: all distance work has to happen in one projected CRS. I use the
  # displaced-point CRS as the target and move the other layers onto it first.
  target_crs <- sf::st_crs(point_sf)
  if (sf::st_crs(feature_sf) != target_crs) {
    feature_sf <- sf::st_transform(feature_sf, target_crs)
  }
  
  if (!is.null(adminBound)) {
    pb_check_sf(adminBound, "adminBound")
    pb_check_crs(adminBound, "adminBound")
    
    if (sf::st_crs(adminBound) != target_crs) {
      adminBound <- sf::st_transform(adminBound, target_crs)
    }
    
    if (nrow(adminBound) > 1) {
      if (is.null(adminID)) {
        stop(
          "adminID must be provided when adminBound has more than one feature.",
          call. = FALSE
        )
      }
      pb_check_cols(adminBound, adminID, "adminBound")
      pb_check_unique_id(sf::st_drop_geometry(adminBound), adminID, "adminBound")
      adminBound$adminID_use <- sf::st_drop_geometry(adminBound)[[adminID]]
    } else {
      adminBound$adminID_use <- if (!is.null(adminID) && adminID %in% names(adminBound)) {
        sf::st_drop_geometry(adminBound)[[adminID]]
      } else {
        1L
      }
    }
  }
  
  # Future-me: I only carry the feature ID and geometry here because that is all
  # the nearest-distance calculation actually needs.
  feature_min <- feature_sf[, feature_id, drop = FALSE]
  feature_ids <- feature_min[[feature_id]]
  
  if (!is.null(long_file)) {
    long_dir <- dirname(long_file)
    if (!dir.exists(long_dir)) {
      dir.create(long_dir, recursive = TRUE, showWarnings = FALSE)
    }
    
    if (file.exists(long_file)) {
      if (isTRUE(overwrite_long_file)) {
        file.remove(long_file)
      } else {
        stop(
          "long_file already exists. Set overwrite_long_file = TRUE to overwrite it.",
          call. = FALSE
        )
      }
    }
  }
  
  # Future-me: if the caller did not specify a chunk size, I use a conservative
  # default whenever the long table is being kept or written, because that is the
  # piece most likely to chew through RAM.
  if (is.null(chunk_size)) {
    chunk_size <- if (isTRUE(return_long) || !is.null(long_file)) 100L else max(1L, nrow(point_sf))
  }
  
  # Future-me: make empty outputs with the right column names up front so I do
  # not have to special-case rbind on empty results later.
  empty_long <- function() {
    out <- data.frame(
      tmp_point_id = point_sf[[point_id]][FALSE],
      tmp_feature_id = feature_ids[FALSE],
      distance_m = numeric(0),
      probability = numeric(0),
      stringsAsFactors = FALSE
    )
    names(out)[1] <- point_id
    names(out)[2] <- feature_id
    out
  }
  
  empty_median <- function() {
    out <- data.frame(
      tmp_point_id = point_sf[[point_id]][FALSE],
      median_distance_m = numeric(0),
      stringsAsFactors = FALSE
    )
    names(out)[1] <- point_id
    out
  }
  
  empty_mode <- function() {
    out <- data.frame(
      tmp_point_id = point_sf[[point_id]][FALSE],
      max_likely_distance_m = numeric(0),
      probability = numeric(0),
      stringsAsFactors = FALSE
    )
    names(out)[1] <- point_id
    out
  }
  
  make_long_na <- function(id_value) {
    out <- data.frame(
      tmp_point_id = id_value,
      tmp_feature_id = feature_ids[NA_integer_][1],
      distance_m = NA_real_,
      probability = NA_real_,
      stringsAsFactors = FALSE
    )
    names(out)[1] <- point_id
    names(out)[2] <- feature_id
    out
  }
  
  make_median_na <- function(id_value) {
    out <- data.frame(
      tmp_point_id = id_value,
      median_distance_m = NA_real_,
      stringsAsFactors = FALSE
    )
    names(out)[1] <- point_id
    out
  }
  
  make_mode_na <- function(id_value) {
    out <- data.frame(
      tmp_point_id = id_value,
      max_likely_distance_m = NA_real_,
      probability = NA_real_,
      stringsAsFactors = FALSE
    )
    names(out)[1] <- point_id
    out
  }
  
  matrix_to_candidate_points <- function(weight_matrix, cell_m, x0, y0, weight_scale = 1) {
    # Future-me: this recreates the same cell centers I would have gotten from
    # the raster path, but I keep everything as plain coordinates until I really
    # need sf geometry for the nearest-feature lookup.
    if (!is.matrix(weight_matrix)) {
      weight_matrix <- as.matrix(weight_matrix)
    }
    
    n_rows <- nrow(weight_matrix)
    n_cols <- ncol(weight_matrix)
    
    x_centers <- x0 - (n_cols / 2 - 0.5) * cell_m + (seq_len(n_cols) - 1) * cell_m
    y_centers <- y0 + (n_rows / 2 - 0.5) * cell_m - (seq_len(n_rows) - 1) * cell_m
    
    grid_index <- expand.grid(
      col = seq_len(n_cols),
      row = seq_len(n_rows)
    )
    
    out <- data.frame(
      x = x_centers[grid_index$col],
      y = y_centers[grid_index$row],
      weight = as.vector(t(weight_matrix)) * weight_scale,
      stringsAsFactors = FALSE
    )
    
    out <- out[is.finite(out$weight) & !is.na(out$weight) & out$weight > 0, , drop = FALSE]
    out
  }
  
  choose_density_components <- function(single_point) {
    # Future-me: newer density buffers can carry mixtures, while the legacy ones
    # just carry one weighted circle. I normalize those two paths here so the
    # rest of the function does not care which object came in.
    if (!density_is_new) {
      return(list(list(weight = 1, radius_m = densityBuffer$radiusMeters)))
    }
    
    mix <- densityBuffer$mixture
    if (is.null(mix) || !isTRUE(mix$enabled)) {
      return(list(list(weight = 1, radius_m = max(densityBuffer$radii_m))))
    }
    
    urban_col <- densityBuffer$urban_col %||% "URBAN_RURA"
    urban_value <- densityBuffer$urban_value %||% "U"
    
    is_urban <- FALSE
    if (urban_col %in% names(single_point)) {
      urban_raw <- as.character(sf::st_drop_geometry(single_point)[[urban_col]][1])
      is_urban <- isTRUE(!is.na(urban_raw) && urban_raw == urban_value)
    }
    
    applies_to <- mix$applies_to %||% "rural_only"
    use_mix <- if (identical(applies_to, "global")) TRUE else !is_urban
    
    if (use_mix) {
      return(Map(
        function(w, r) list(weight = w, radius_m = r),
        mix$weights,
        mix$radii_m
      ))
    }
    
    # Future-me: if the mixture is rural-only and there is a dedicated urban
    # radius stored separately, I use that instead of accidentally falling back
    # to the biggest radius in the object.
    non_mix_radii <- setdiff(densityBuffer$radii_m, mix$radii_m)
    fallback_radius <- if (length(non_mix_radii) > 0) {
      max(non_mix_radii)
    } else {
      min(densityBuffer$radii_m)
    }
    
    list(list(weight = 1, radius_m = fallback_radius))
  }
  
  get_trim_polygon <- function(single_point) {
    if (is.null(adminBound)) {
      return(NULL)
    }
    
    if (nrow(adminBound) == 1) {
      return(adminBound)
    }
    
    point_min <- single_point[, character(0), drop = FALSE]
    joined_admin <- suppressWarnings(
      sf::st_join(
        point_min,
        adminBound[, "adminID_use", drop = FALSE],
        join = sf::st_intersects,
        left = TRUE
      )
    )
    
    admin_value <- sf::st_drop_geometry(joined_admin)$adminID_use[1]
    if (is.na(admin_value)) {
      return(NULL)
    }
    
    adminBound[adminBound$adminID_use == admin_value, , drop = FALSE]
  }
  
  build_candidate_support <- function(single_point) {
    xy <- suppressWarnings(sf::st_coordinates(single_point))
    if (nrow(xy) != 1 || anyNA(xy[1, c("X", "Y")])) {
      return(NULL)
    }
    
    x0 <- xy[1, "X"]
    y0 <- xy[1, "Y"]
    
    components <- choose_density_components(single_point)
    candidate_parts <- lapply(components, function(comp) {
      if (density_is_new) {
        radius_key <- as.character(comp$radius_m)
        if (!radius_key %in% names(densityBuffer$kernels)) {
          stop(
            "densityBuffer is missing a stored kernel for radius ",
            comp$radius_m,
            ".",
            call. = FALSE
          )
        }
        
        matrix_to_candidate_points(
          weight_matrix = densityBuffer$kernels[[radius_key]],
          cell_m = densityBuffer$cell_m,
          x0 = x0,
          y0 = y0,
          weight_scale = comp$weight
        )
      } else {
        matrix_to_candidate_points(
          weight_matrix = densityBuffer$weightedCircle,
          cell_m = densityBuffer$cellMeters,
          x0 = x0,
          y0 = y0,
          weight_scale = comp$weight
        )
      }
    })
    
    candidate_df <- do.call(rbind, candidate_parts)
    if (is.null(candidate_df) || nrow(candidate_df) == 0) {
      return(NULL)
    }
    
    # Future-me: mixture components can overlap on the same support point, so I
    # add those weights together before I do trimming or thresholding.
    candidate_df <- stats::aggregate(
      candidate_df$weight,
      by = list(x = candidate_df$x, y = candidate_df$y),
      FUN = sum
    )
    names(candidate_df) <- c("x", "y", "weight")
    
    candidate_sf <- sf::st_as_sf(candidate_df, coords = c("x", "y"), crs = sf::st_crs(single_point))
    
    trim_poly <- get_trim_polygon(single_point)
    if (!is.null(trim_poly)) {
      keep_idx <- lengths(sf::st_intersects(candidate_sf, trim_poly)) > 0
      candidate_sf <- candidate_sf[keep_idx, , drop = FALSE]
    }
    
    if (nrow(candidate_sf) == 0) {
      return(NULL)
    }
    
    if (min_weight > 0) {
      candidate_sf <- candidate_sf[candidate_sf$weight > min_weight, , drop = FALSE]
    }
    
    if (nrow(candidate_sf) == 0) {
      return(NULL)
    }
    
    # Future-me: once support points have been trimmed and tiny weights possibly
    # dropped, I renormalize so everything downstream stays on a true probability scale.
    candidate_sf$weight <- candidate_sf$weight / sum(candidate_sf$weight)
    candidate_sf
  }
  
  maybe_round_distance <- function(x) {
    if (is.null(distance_digits)) {
      return(x)
    }
    round(x, digits = distance_digits)
  }
  
  one_point <- function(row_index) {
    row_point_id <- point_sf[[point_id]][row_index]
    single_point <- point_sf[row_index, , drop = FALSE]
    
    # Future-me: if the source geometry is empty or malformed, I return an
    # explicit missing result instead of letting that silently vanish later.
    if (isTRUE(sf::st_is_empty(single_point))) {
      return(list(
        nearest_distance_probabilities = if (isTRUE(return_long) || !is.null(long_file)) make_long_na(row_point_id) else NULL,
        median_distance = make_median_na(row_point_id),
        max_likely_distance = make_mode_na(row_point_id)
      ))
    }
    
    candidate_sf <- build_candidate_support(single_point)
    
    if (is.null(candidate_sf) || nrow(candidate_sf) == 0) {
      return(list(
        nearest_distance_probabilities = if (isTRUE(return_long) || !is.null(long_file)) make_long_na(row_point_id) else NULL,
        median_distance = make_median_na(row_point_id),
        max_likely_distance = make_mode_na(row_point_id)
      ))
    }
    
    # Future-me: I only ever need the nearest feature for each candidate point,
    # not a full distance matrix. That keeps this substantially lighter than the
    # old all-by-all approach.
    nearest_idx <- sf::st_nearest_feature(candidate_sf, feature_min)
    nearest_feature_id <- feature_ids[nearest_idx]
    
    nearest_distance <- sf::st_distance(
      sf::st_geometry(candidate_sf),
      sf::st_geometry(feature_min)[nearest_idx],
      by_element = TRUE
    )
    nearest_distance_m <- as.numeric(nearest_distance)
    weight_vec <- candidate_sf$weight
    
    # Future-me: the weighted median should be based on the raw distances, not
    # the rounded buckets, because this is meant to be the probability-weighted
    # central tendency and not a histogram summary.
    median_df <- data.frame(
      tmp_point_id = row_point_id,
      median_distance_m = as.numeric(
        pb__weighted_quantile(
          x = nearest_distance_m,
          w = weight_vec,
          probs = 0.5
        )[[1]]
      ),
      stringsAsFactors = FALSE
    )
    names(median_df)[1] <- point_id
    
    distance_bucket_m <- maybe_round_distance(nearest_distance_m)
    
    # Future-me: the mode is intentionally calculated on the rounded distance
    # buckets if distance_digits is supplied. That makes the result much more
    # stable and much easier to interpret than an exact-distance mode.
    mode_agg <- stats::aggregate(
      weight_vec,
      by = list(distance_m = distance_bucket_m),
      FUN = sum
    )
    names(mode_agg) <- c("distance_m", "probability")
    
    mode_agg <- mode_agg[
      order(-mode_agg$probability, mode_agg$distance_m),
      ,
      drop = FALSE
    ]
    
    # Future-me: if two distance buckets tie in probability, I keep the smaller
    # distance so the tie-break is deterministic and conservative.
    mode_df <- data.frame(
      tmp_point_id = row_point_id,
      max_likely_distance_m = mode_agg$distance_m[1],
      probability = mode_agg$probability[1],
      stringsAsFactors = FALSE
    )
    names(mode_df)[1] <- point_id
    
    long_df <- NULL
    if (isTRUE(return_long) || !is.null(long_file)) {
      long_df <- data.frame(
        tmp_feature_id = nearest_feature_id,
        distance_m = distance_bucket_m,
        probability = weight_vec,
        stringsAsFactors = FALSE
      )
      
      # Future-me: I aggregate here by feature-distance pair so the long output
      # is the actual probability mass for each combination, not one row per
      # support cell. That alone can save a lot of storage.
      long_df <- stats::aggregate(
        long_df$probability,
        by = list(
          tmp_feature_id = long_df$tmp_feature_id,
          distance_m = long_df$distance_m
        ),
        FUN = sum
      )
      names(long_df) <- c(feature_id, "distance_m", "probability")
      long_df[[point_id]] <- row_point_id
      long_df <- long_df[, c(point_id, feature_id, "distance_m", "probability"), drop = FALSE]
      long_df <- long_df[
        order(-long_df$probability, long_df$distance_m),
        ,
        drop = FALSE
      ]
      rownames(long_df) <- NULL
    }
    
    list(
      nearest_distance_probabilities = long_df,
      median_distance = median_df,
      max_likely_distance = mode_df
    )
  }
  
  if (nrow(point_sf) == 0) {
    out <- list(
      nearest_distance_probabilities = if (isTRUE(return_long)) empty_long() else NULL,
      median_distance = empty_median(),
      max_likely_distance = empty_mode()
    )
    if (!is.null(long_file)) {
      attr(out, "nearest_distance_probabilities_file") <- long_file
    }
    return(out)
  }
  
  row_index <- seq_len(nrow(point_sf))
  chunk_ids <- split(row_index, ceiling(seq_along(row_index) / chunk_size))
  use_parallel <- isTRUE(parallel) || n.cores > 1
  
  if (use_parallel) {
    # Future-me: the work is embarrassingly parallel across displaced points, so
    # I let future handle the split but restore the incoming plan on the way out.
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    future::plan(future::multisession, workers = n.cores)
  }
  
  long_chunks <- if (isTRUE(return_long)) list() else NULL
  median_chunks <- list()
  mode_chunks <- list()
  wrote_long_header <- FALSE
  
  for (chunk_i in seq_along(chunk_ids)) {
    this_chunk <- chunk_ids[[chunk_i]]
    
    chunk_results <- if (use_parallel) {
      future.apply::future_lapply(
        X = this_chunk,
        FUN = one_point,
        future.seed = future.seed,
        future.scheduling = future.scheduling
      )
    } else {
      lapply(this_chunk, one_point)
    }
    
    chunk_median <- do.call(
      rbind,
      lapply(chunk_results, `[[`, "median_distance")
    )
    chunk_mode <- do.call(
      rbind,
      lapply(chunk_results, `[[`, "max_likely_distance")
    )
    
    rownames(chunk_median) <- NULL
    rownames(chunk_mode) <- NULL
    
    median_chunks[[length(median_chunks) + 1L]] <- chunk_median
    mode_chunks[[length(mode_chunks) + 1L]] <- chunk_mode
    
    if (isTRUE(return_long) || !is.null(long_file)) {
      chunk_long_list <- lapply(chunk_results, `[[`, "nearest_distance_probabilities")
      chunk_long_list <- Filter(Negate(is.null), chunk_long_list)
      
      if (length(chunk_long_list) > 0) {
        chunk_long <- do.call(rbind, chunk_long_list)
        rownames(chunk_long) <- NULL
      } else {
        chunk_long <- empty_long()
      }
      
      if (!is.null(long_file) && nrow(chunk_long) > 0) {
        utils::write.table(
          x = chunk_long,
          file = long_file,
          sep = ",",
          row.names = FALSE,
          col.names = !wrote_long_header,
          append = wrote_long_header,
          quote = TRUE,
          qmethod = "double"
        )
        wrote_long_header <- TRUE
      }
      
      if (isTRUE(return_long)) {
        long_chunks[[length(long_chunks) + 1L]] <- chunk_long
      }
    }
    
    # Future-me: on very large jobs this helps release chunk-local objects before
    # the next chunk is built, especially when the long output is being written.
    rm(chunk_results)
  }
  
  median_distance <- if (length(median_chunks) > 0) {
    out <- do.call(rbind, median_chunks)
    rownames(out) <- NULL
    out
  } else {
    empty_median()
  }
  
  max_likely_distance <- if (length(mode_chunks) > 0) {
    out <- do.call(rbind, mode_chunks)
    rownames(out) <- NULL
    out
  } else {
    empty_mode()
  }
  
  nearest_distance_probabilities <- if (isTRUE(return_long)) {
    if (length(long_chunks) > 0) {
      out <- do.call(rbind, long_chunks)
      rownames(out) <- NULL
      out
    } else {
      empty_long()
    }
  } else {
    NULL
  }
  
  out <- list(
    nearest_distance_probabilities = nearest_distance_probabilities,
    median_distance = median_distance,
    max_likely_distance = max_likely_distance
  )
  
  if (!is.null(long_file)) {
    attr(out, "nearest_distance_probabilities_file") <- long_file
  }
  
  out
}
