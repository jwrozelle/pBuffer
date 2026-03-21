#' Probability-aware nearest-point distance join under displacement uncertainty
#'
#' @description
#' For each displaced point, this function builds the set of plausible true
#' locations implied by a `pb_densityBuffer()` object, treats those as weighted
#' candidate `sf` points, and then assigns each candidate location to its
#' nearest point or multipoint feature.
#'
#' The function then summarizes that nearest-feature assignment in three ways:
#' \enumerate{
#'   \item an optional long table giving the displaced-point ID, nearest-feature
#'   ID, nearest distance in meters, and the probability mass for that
#'   feature-distance pair;
#'   \item a one-row summary per displaced point giving the probability-weighted
#'   median nearest distance; and
#'   \item a one-row summary per displaced point giving the most likely nearest
#'   distance after optional rounding or bucketing with `round()`.
#' }
#'
#' This rewrite is intentionally direct:
#' \itemize{
#'   \item no raster path,
#'   \item no `templateRaster`,
#'   \item no legacy density-buffer objects,
#'   \item and no splitting `sf` objects themselves for parallel work.
#' }
#'
#' @param point_sf An `sf` POINT layer of displaced points.
#' @param feature_sf An `sf` POINT or MULTIPOINT layer of candidate nearest
#'   features.
#' @param densityBuffer A `pb_densityBuffer` object created by
#'   `pb_densityBuffer()`.
#' @param point_id Column name giving the unique displaced-point ID.
#' @param feature_id Column name giving the unique nearest-feature ID.
#' @param min_weight Drop candidate true-location support points with weights
#'   less than or equal to this threshold before assigning nearest features,
#'   then renormalize the remaining weights so the probabilities still sum to 1.
#' @param adminBound Optional administrative boundary used to trim plausible
#'   true locations before assigning nearest features.
#' @param adminID Column name giving the admin ID in `adminBound`. If this
#'   column is also present in `point_sf`, the function uses it directly to
#'   subset the matching admin polygon for each displaced point. Otherwise, if
#'   `adminBound` has more than one row, the function falls back to a spatial
#'   match between the displaced point and `adminBound`.
#' @param distance_digits Optional single numeric value passed to
#'   `round(x, digits = distance_digits)` when aggregating nearest distances for
#'   the long output and for the maximum-likelihood distance. Use `NULL` to keep
#'   exact distances, positive values for decimal places, `0` for whole meters,
#'   and negative values for coarser buckets such as 10 m or 100 m.
#' @param return_long Logical. If `TRUE`, return the long
#'   feature-distance-probability table in memory. If `FALSE`, skip storing it
#'   in memory.
#' @param long_file Optional file path for writing the long table to disk in
#'   chunked CSV form. This is useful when the long output is too large to keep
#'   in memory. If supplied together with `return_long = TRUE`, the long table
#'   is both written to disk and returned in memory.
#' @param overwrite_long_file Logical. If `TRUE` and `long_file` already exists,
#'   overwrite it. If `FALSE`, error instead.
#' @param chunk_size Optional number of displaced points to process per chunk.
#'   If `NULL`, the function chooses a conservative default whenever the long
#'   table is being kept or written, and otherwise processes everything in one
#'   chunk.
#' @param n.cores Number of workers to use when parallel execution is requested.
#' @param parallel Logical. If `TRUE`, evaluate displaced points in parallel
#'   with `future.apply`. If `FALSE`, evaluation is serial.
#' @param future.seed Seed handling passed to `future.apply::future_lapply()`.
#' @param future.scheduling Scheduling argument passed to
#'   `future.apply::future_lapply()`.
#'
#' @returns A list with three elements:
#' \describe{
#'   \item{nearest_distance_probabilities}{A data frame with four columns:
#'   displaced-point ID, nearest-feature ID, nearest distance in meters, and the
#'   probability for that feature-distance pair. If `return_long = FALSE`, this
#'   element is returned as `NULL`. If `long_file` is supplied, this same table
#'   is also written to disk in chunks.}
#'   \item{median_distance}{A data frame with two columns: displaced-point ID
#'   and the probability-weighted median nearest distance in meters.}
#'   \item{max_likely_distance}{A data frame with three columns:
#'   displaced-point ID, the maximum-likelihood nearest distance in meters after
#'   optional rounding via `distance_digits`, and the probability mass of that
#'   distance bucket.}
#' }
#'
#' @author J.W. Rozelle
#'
#' @export
pb_pointProbJoin <- function(point_sf,
                             feature_sf,
                             densityBuffer,
                             point_id = "DHSID",
                             feature_id = "SPAid",
                             min_weight = 0,
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
  
  # ---- helper: rebuild sf explicitly if the geometry pointer got messy ----
  repair_sf <- function(x) {
    if (!inherits(x, "sf")) {
      stop("Expected an sf object.", call. = FALSE)
    }
    
    geom <- sf::st_geometry(x)
    geom_col <- attr(x, "sf_column")
    
    if (is.null(geom_col) || !geom_col %in% names(x)) {
      geom_col <- "geometry"
    }
    
    x_df <- sf::st_drop_geometry(x)
    x_df[[geom_col]] <- geom
    
    x <- sf::st_as_sf(
      x_df,
      sf_column_name = geom_col,
      crs = sf::st_crs(geom)
    )
    
    sf::st_geometry(x) <- geom_col
    x
  }
  
  # ---- helper: keep only the columns I actually need and preserve sf ----
  keep_sf_cols <- function(x, cols) {
    x <- repair_sf(x)
    
    missing_cols <- setdiff(cols, names(x))
    if (length(missing_cols) > 0) {
      stop(
        "These required columns are missing: ",
        paste(missing_cols, collapse = ", "),
        call. = FALSE
      )
    }
    
    geom_col <- attr(x, "sf_column")
    x <- x[, unique(c(cols, geom_col)), drop = FALSE]
    sf::st_geometry(x) <- geom_col
    x
  }
  
  # ---- helper: simple weighted median on the raw distances ----
  weighted_median <- function(x, w) {
    keep <- is.finite(x) & is.finite(w) & w > 0
    x <- x[keep]
    w <- w[keep]
    
    if (length(x) == 0) {
      return(NA_real_)
    }
    
    ord <- order(x)
    x <- x[ord]
    w <- w[ord]
    
    cw <- cumsum(w) / sum(w)
    x[which(cw >= 0.5)[1]]
  }
  
  # ---- helper: precompute offset support for one kernel matrix ----
  build_kernel_support <- function(weight_matrix, cell_m) {
    # Future me:
    # I only want to compute the origin-centered support grid once per stored
    # kernel radius. After that, each displaced point is just x0 + dx, y0 + dy.
    if (!is.matrix(weight_matrix)) {
      weight_matrix <- as.matrix(weight_matrix)
    }
    
    n_rows <- nrow(weight_matrix)
    n_cols <- ncol(weight_matrix)
    
    x_centers <- - (n_cols / 2 - 0.5) * cell_m + (seq_len(n_cols) - 1) * cell_m
    y_centers <-   (n_rows / 2 - 0.5) * cell_m - (seq_len(n_rows) - 1) * cell_m
    
    grid_index <- expand.grid(
      col = seq_len(n_cols),
      row = seq_len(n_rows)
    )
    
    out <- data.frame(
      dx = x_centers[grid_index$col],
      dy = y_centers[grid_index$row],
      weight = as.vector(t(weight_matrix)),
      stringsAsFactors = FALSE
    )
    
    out <- out[is.finite(out$weight) & !is.na(out$weight) & out$weight > 0, , drop = FALSE]
    
    if (nrow(out) == 0) {
      return(out)
    }
    
    # Future me:
    # I normalize the stored support now so each kernel is cleanly on a
    # probability scale before I apply any mixture weight or admin trimming.
    out$weight <- out$weight / sum(out$weight)
    rownames(out) <- NULL
    out
  }
  
  # ---- helper: figure out which radius / weights apply to this point ----
  choose_density_components <- function(single_point) {
    mix <- densityBuffer$mixture
    
    # Future me:
    # For the clean non-mixture path, I want exactly one radius in the object.
    if (is.null(mix) || !isTRUE(mix$enabled)) {
      if (length(densityBuffer$radii_m) != 1) {
        stop(
          "When densityBuffer$mixture is not enabled, densityBuffer should contain ",
          "exactly one radius. Build separate buffers explicitly if you want ",
          "different urban and rural radii.",
          call. = FALSE
        )
      }
      
      return(list(
        list(weight = 1, radius_m = densityBuffer$radii_m[1])
      ))
    }
    
    urban_col <- if (!is.null(densityBuffer$urban_col)) densityBuffer$urban_col else "URBAN_RURA"
    urban_value <- if (!is.null(densityBuffer$urban_value)) densityBuffer$urban_value else "U"
    
    is_urban <- FALSE
    if (urban_col %in% names(single_point)) {
      urban_raw <- as.character(sf::st_drop_geometry(single_point)[[urban_col]][1])
      is_urban <- !is.na(urban_raw) && identical(urban_raw, urban_value)
    }
    
    applies_to <- if (!is.null(mix$applies_to)) mix$applies_to else "rural_only"
    use_mix <- if (identical(applies_to, "global")) TRUE else !is_urban
    
    if (isTRUE(use_mix)) {
      return(Map(
        function(w, r) list(weight = w, radius_m = r),
        mix$weights,
        mix$radii_m
      ))
    }
    
    # Future me:
    # If the mixture only applies to rural points, the object needs to contain
    # exactly one non-mixture radius to use for urban points.
    non_mix_radii <- setdiff(densityBuffer$radii_m, mix$radii_m)
    
    if (length(non_mix_radii) != 1) {
      stop(
        "This densityBuffer does not contain a single non-mixture fallback radius ",
        "for non-mixture points. Build separate urban and rural density buffers ",
        "explicitly instead.",
        call. = FALSE
      )
    }
    
    list(
      list(weight = 1, radius_m = non_mix_radii[1])
    )
  }
  
  # ---- helper: get the admin polygon for this displaced point ----
  get_trim_polygon <- function(single_point) {
    if (is.null(admin_min)) {
      return(NULL)
    }
    
    if (nrow(admin_min) == 1) {
      return(admin_min)
    }
    
    # Future me:
    # Fast path first: if the point layer already carries adminID, use it.
    if (!is.null(adminID) && adminID %in% names(single_point)) {
      this_admin_val <- sf::st_drop_geometry(single_point)[[adminID]][1]
      
      if (!is.na(this_admin_val)) {
        this_admin <- admin_min[admin_min[[adminID]] == this_admin_val, , drop = FALSE]
        if (nrow(this_admin) > 0) {
          return(repair_sf(this_admin))
        }
      }
    }
    
    # Future me:
    # Fallback path: if the point layer does not carry adminID, spatially match
    # the displaced point to the admin polygon that contains / touches it.
    point_stub <- single_point[, character(0), drop = FALSE]
    
    joined <- suppressWarnings(
      sf::st_join(
        point_stub,
        admin_min[, adminID, drop = FALSE],
        join = sf::st_intersects,
        left = TRUE
      )
    )
    
    this_admin_val <- sf::st_drop_geometry(joined)[[adminID]][1]
    
    if (is.na(this_admin_val)) {
      return(NULL)
    }
    
    this_admin <- admin_min[admin_min[[adminID]] == this_admin_val, , drop = FALSE]
    if (nrow(this_admin) == 0) {
      return(NULL)
    }
    
    repair_sf(this_admin)
  }
  
  # ---- helper: build the weighted candidate support for one displaced point ----
  build_candidate_support <- function(single_point) {
    xy <- suppressWarnings(sf::st_coordinates(single_point))
    
    if (nrow(xy) != 1 || anyNA(xy[1, c("X", "Y")])) {
      return(NULL)
    }
    
    x0 <- xy[1, "X"]
    y0 <- xy[1, "Y"]
    
    components <- choose_density_components(single_point)
    
    candidate_parts <- lapply(components, function(comp) {
      radius_key <- as.character(comp$radius_m)
      
      if (!radius_key %in% names(kernel_supports)) {
        stop(
          "densityBuffer is missing a stored kernel for radius ",
          comp$radius_m,
          ".",
          call. = FALSE
        )
      }
      
      support_df <- kernel_supports[[radius_key]]
      
      data.frame(
        x = x0 + support_df$dx,
        y = y0 + support_df$dy,
        weight = support_df$weight * comp$weight,
        stringsAsFactors = FALSE
      )
    })
    
    candidate_df <- do.call(rbind, candidate_parts)
    
    if (is.null(candidate_df) || nrow(candidate_df) == 0) {
      return(NULL)
    }
    
    # Future me:
    # Mixture components can overlap on the same support point, so I collapse
    # those duplicates before I do any trimming.
    candidate_df <- stats::aggregate(
      candidate_df$weight,
      by = list(x = candidate_df$x, y = candidate_df$y),
      FUN = sum
    )
    names(candidate_df) <- c("x", "y", "weight")
    
    candidate_sf <- sf::st_as_sf(
      candidate_df,
      coords = c("x", "y"),
      crs = sf::st_crs(single_point)
    )
    
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
    
    # Future me:
    # After admin trimming and min_weight filtering, I need to renormalize or
    # the downstream probabilities will no longer sum to 1.
    candidate_sf$weight <- candidate_sf$weight / sum(candidate_sf$weight)
    repair_sf(candidate_sf)
  }
  
  # ---- helper: optional distance bucketing for the long table and mode ----
  maybe_round_distance <- function(x) {
    if (is.null(distance_digits)) {
      return(x)
    }
    round(x, digits = distance_digits)
  }
  
  # ---- helper: empty outputs with the right column names ----
  empty_long <- function() {
    out <- data.frame(
      tmp_point_id = point_min[[point_id]][FALSE],
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
      tmp_point_id = point_min[[point_id]][FALSE],
      median_distance_m = numeric(0),
      stringsAsFactors = FALSE
    )
    names(out)[1] <- point_id
    out
  }
  
  empty_mode <- function() {
    out <- data.frame(
      tmp_point_id = point_min[[point_id]][FALSE],
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
  
  # ---- helper: evaluate one displaced point ----
  one_point <- function(row_index) {
    row_point_id <- point_min[[point_id]][row_index]
    single_point <- repair_sf(point_min[row_index, , drop = FALSE])
    
    # Future me:
    # If the geometry is empty or malformed, I want an explicit NA result rather
    # than silently dropping the observation.
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
    
    # Future me:
    # I only ever need the nearest feature for each candidate point, not a full
    # all-by-all distance matrix. That is the main reason this stays manageable.
    nearest_idx <- sf::st_nearest_feature(candidate_sf, feature_min)
    nearest_feature_id <- feature_ids[nearest_idx]
    
    nearest_distance_m <- as.numeric(
      sf::st_distance(
        candidate_sf,
        feature_min[nearest_idx, , drop = FALSE],
        by_element = TRUE
      )
    )
    
    weight_vec <- candidate_sf$weight
    
    median_df <- data.frame(
      tmp_point_id = row_point_id,
      median_distance_m = weighted_median(nearest_distance_m, weight_vec),
      stringsAsFactors = FALSE
    )
    names(median_df)[1] <- point_id
    
    distance_bucket_m <- maybe_round_distance(nearest_distance_m)
    
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
    
    # Future me:
    # If two distance buckets tie in probability, I keep the smaller distance so
    # the tie-break is deterministic and conservative.
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
      
      # Future me:
      # I aggregate by feature-distance pair so the long output is actual mass
      # on each combination, not one row per support cell.
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
  
  # ---- validate the main inputs ----
  point_sf <- repair_sf(point_sf)
  feature_sf <- repair_sf(feature_sf)
  
  pb_check_sf(point_sf, "point_sf")
  pb_check_geom_type(point_sf, "POINT", "point_sf")
  pb_check_crs(point_sf, "point_sf")
  pb_check_projected(point_sf, "point_sf")
  pb_check_cols(point_sf, point_id, "point_sf")
  pb_check_unique_id(sf::st_drop_geometry(point_sf), point_id, "point_sf")
  
  pb_check_sf(feature_sf, "feature_sf")
  pb_check_crs(feature_sf, "feature_sf")
  pb_check_cols(feature_sf, feature_id, "feature_sf")
  pb_check_unique_id(sf::st_drop_geometry(feature_sf), feature_id, "feature_sf")
  
  feature_geom <- unique(as.character(sf::st_geometry_type(feature_sf)))
  if (!all(feature_geom %in% c("POINT", "MULTIPOINT"))) {
    stop("feature_sf must have POINT or MULTIPOINT geometries.", call. = FALSE)
  }
  
  if (nrow(feature_sf) == 0) {
    stop("feature_sf must contain at least one feature.", call. = FALSE)
  }
  
  if (!inherits(densityBuffer, "pb_densityBuffer")) {
    stop(
      "densityBuffer must be a pb_densityBuffer object created by pb_densityBuffer().",
      call. = FALSE
    )
  }
  
  if (is.null(densityBuffer$kernels) || !is.list(densityBuffer$kernels) || length(densityBuffer$kernels) == 0) {
    stop("densityBuffer$kernels is missing or empty.", call. = FALSE)
  }
  
  if (is.null(densityBuffer$cell_m) || !is.numeric(densityBuffer$cell_m) ||
      length(densityBuffer$cell_m) != 1 || is.na(densityBuffer$cell_m) || densityBuffer$cell_m <= 0) {
    stop("densityBuffer$cell_m must be a single positive number.", call. = FALSE)
  }
  
  pb_check_pos_scalar(n.cores, "n.cores")
  
  if (!is.numeric(min_weight) || length(min_weight) != 1 || is.na(min_weight) || min_weight < 0) {
    stop("min_weight must be a single non-negative number.", call. = FALSE)
  }
  
  if (!is.null(distance_digits)) {
    if (!is.numeric(distance_digits) || length(distance_digits) != 1 ||
        is.na(distance_digits) || !is.finite(distance_digits)) {
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
    if (!is.numeric(chunk_size) || length(chunk_size) != 1 ||
        is.na(chunk_size) || chunk_size <= 0) {
      stop("chunk_size must be NULL or a single positive number.", call. = FALSE)
    }
    chunk_size <- as.integer(ceiling(chunk_size))
  }
  
  if (!is.logical(parallel) || length(parallel) != 1 || is.na(parallel)) {
    stop("parallel must be TRUE or FALSE.", call. = FALSE)
  }
  
  # ---- bring everything into the displaced-point CRS ----
  target_crs <- sf::st_crs(point_sf)
  
  if (sf::st_crs(feature_sf) != target_crs) {
    feature_sf <- sf::st_transform(feature_sf, target_crs)
    feature_sf <- repair_sf(feature_sf)
  }
  
  admin_min <- NULL
  
  if (!is.null(adminBound)) {
    adminBound <- repair_sf(adminBound)
    
    pb_check_sf(adminBound, "adminBound")
    pb_check_crs(adminBound, "adminBound")
    
    if (sf::st_crs(adminBound) != target_crs) {
      adminBound <- sf::st_transform(adminBound, target_crs)
      adminBound <- repair_sf(adminBound)
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
      admin_min <- keep_sf_cols(adminBound, adminID)
    } else {
      if (!is.null(adminID) && adminID %in% names(adminBound)) {
        admin_min <- keep_sf_cols(adminBound, adminID)
      } else {
        admin_min <- repair_sf(adminBound)
      }
    }
  }
  
  # ---- keep only the displaced-point columns I actually need ----
  point_keep <- point_id
  
  urban_col <- if (!is.null(densityBuffer$urban_col)) densityBuffer$urban_col else NULL
  if (!is.null(urban_col) && urban_col %in% names(point_sf)) {
    point_keep <- c(point_keep, urban_col)
  }
  
  if (!is.null(admin_min) && !is.null(adminID) && adminID %in% names(point_sf)) {
    point_keep <- c(point_keep, adminID)
  }
  
  point_min <- keep_sf_cols(point_sf, unique(point_keep))
  feature_min <- keep_sf_cols(feature_sf, feature_id)
  feature_ids <- feature_min[[feature_id]]
  
  # ---- precompute support tables once per kernel radius ----
  kernel_supports <- lapply(
    densityBuffer$kernels,
    build_kernel_support,
    cell_m = densityBuffer$cell_m
  )
  
  # ---- if the caller did not set chunk_size, choose a conservative default ----
  if (is.null(chunk_size)) {
    chunk_size <- if (isTRUE(return_long) || !is.null(long_file)) {
      100L
    } else {
      max(1L, nrow(point_min))
    }
  }
  
  # ---- prepare the optional on-disk long output ----
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
  
  # ---- quick exit for an empty point layer ----
  if (nrow(point_min) == 0) {
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
  
  # ---- chunk the row indices, not the sf object itself ----
  # Future me:
  # This is the big geometry-stability choice. I am not splitting point_sf into
  # a bunch of mini sf objects and shipping those around. I chunk row indices,
  # then subset inside the worker / loop.
  row_index <- seq_len(nrow(point_min))
  chunk_ids <- split(row_index, ceiling(seq_along(row_index) / chunk_size))
  use_parallel <- isTRUE(parallel) && n.cores > 1
  
  if (use_parallel) {
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
    
    # Future me:
    # On very large jobs this helps release chunk-local objects before the next
    # chunk is built.
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