#' @name pb_pointProbJoin
#' @rdname pb_pointProbJoin
#' @title pb_pointProbJoin
#'
#' @description
#' For each displaced point, apply a spatial uncertainty kernel, identify the
#' nearest feature to each candidate true location, and summarise the resulting
#' probability distribution of nearest-feature IDs and nearest distances.
#'
#' This version is written to be memory-conscious and robust to sf geometry
#' issues in parallel workflows. In particular, it avoids splitting sf objects
#' directly for parallel processing, because that can leave the sf geometry
#' pointer in a broken state on workers.
#'
#' @param point_sf sf POINT object of displaced points.
#' @param feature_sf sf POINT or MULTIPOINT object of candidate nearest features.
#' @param densityBuffer A pb_densityBuffer object, or legacy density buffer list.
#' @param point_id Character scalar naming the unique displaced-point ID column.
#' @param feature_id Character scalar naming the unique feature ID column.
#' @param min_weight Optional numeric threshold; support cells below this
#'   probability are dropped before nearest-feature assignment.
#' @param adminBound Optional sf polygon layer used to constrain candidate true
#'   locations to the corresponding administrative unit.
#' @param adminID Character scalar naming the admin ID column in adminBound and
#'   point_sf (when adminBound is used).
#' @param distance_digits Optional integer passed to round() for the
#'   max-likely-distance summary and, if returned, the long nearest-distance
#'   table. Use negative values for coarser bucketing (e.g. -1 = nearest 10 m).
#' @param return_long Logical; if TRUE, return the long nearest-distance
#'   probability table in memory.
#' @param long_file Optional file path. If provided, the long nearest-distance
#'   table is written incrementally to CSV.
#' @param overwrite_long_file Logical; if TRUE and long_file exists, overwrite it.
#' @param chunk_size Integer number of displaced points to process per chunk.
#' @param n.cores Integer number of workers for optional parallel processing.
#' @param parallel Logical; if TRUE, process chunks in parallel.
#' @param weight_name Deprecated/ignored, retained for compatibility.
#' @param dist_quantiles Deprecated/ignored, retained for compatibility.
#' @param dist_breaks Deprecated/ignored, retained for compatibility.
#'
#' @returns A list with three items:
#' \describe{
#'   \item{nearest_distance_probabilities}{A data frame with four columns:
#'   displaced point ID, nearest feature ID, nearest distance in meters, and
#'   probability. Returned as NULL when return_long = FALSE.}
#'   \item{median_distance}{A data frame with two columns: displaced point ID
#'   and the probability-weighted median nearest distance in meters.}
#'   \item{max_likely_distance}{A data frame with three columns: displaced point
#'   ID, the maximum-likelihood nearest distance bucket, and the probability
#'   mass for that bucket.}
#' }
#'
#' @export
pb_pointProbJoin <- function(
    point_sf,
    feature_sf,
    densityBuffer,
    point_id = "DHSID",
    feature_id = "facID",
    min_weight = 0,
    adminBound = NULL,
    adminID = NULL,
    distance_digits = NULL,
    return_long = TRUE,
    long_file = NULL,
    overwrite_long_file = FALSE,
    chunk_size = 100,
    n.cores = 1,
    parallel = FALSE,
    weight_name = NULL,
    dist_quantiles = NULL,
    dist_breaks = NULL
) {
  
  # Future me:
  # These old arguments are kept only so older code does not fail.
  # They are intentionally ignored in this rewrite.
  invisible(weight_name)
  invisible(dist_quantiles)
  invisible(dist_breaks)
  
  # ---- helper: rebuild sf explicitly ------------------------------------
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
  
  # ---- helper: keep only needed columns while preserving sf -------------
  keep_sf_cols <- function(x, cols) {
    geom_col <- attr(x, "sf_column")
    if (is.null(geom_col) || !geom_col %in% names(x)) {
      x <- repair_sf(x)
      geom_col <- attr(x, "sf_column")
    }
    
    miss <- setdiff(cols, names(x))
    if (length(miss) > 0) {
      stop(
        "These required columns are missing: ",
        paste(miss, collapse = ", "),
        call. = FALSE
      )
    }
    
    x <- x[, unique(c(cols, geom_col)), drop = FALSE]
    sf::st_geometry(x) <- geom_col
    x
  }
  
  # ---- helper: weighted median ------------------------------------------
  weighted_median <- function(x, w) {
    keep <- !(is.na(x) | is.na(w) | w <= 0)
    x <- x[keep]
    w <- w[keep]
    
    if (length(x) == 0) {
      return(NA_real_)
    }
    
    ord <- order(x)
    x <- x[ord]
    w <- w[ord]
    
    w <- w / sum(w)
    cw <- cumsum(w)
    
    x[which(cw >= 0.5)[1]]
  }
  
  # ---- helper: empty outputs --------------------------------------------
  empty_output <- function(return_long, point_id, feature_id) {
    long_df <- NULL
    if (isTRUE(return_long)) {
      long_df <- data.frame(
        tmp_point = character(0),
        tmp_feature = character(0),
        distance_m = numeric(0),
        probability = numeric(0),
        stringsAsFactors = FALSE
      )
      names(long_df)[1] <- point_id
      names(long_df)[2] <- feature_id
    }
    
    median_df <- data.frame(
      tmp_point = character(0),
      median_distance_m = numeric(0),
      stringsAsFactors = FALSE
    )
    names(median_df)[1] <- point_id
    
    mode_df <- data.frame(
      tmp_point = character(0),
      max_likely_distance_m = numeric(0),
      probability = numeric(0),
      stringsAsFactors = FALSE
    )
    names(mode_df)[1] <- point_id
    
    list(
      nearest_distance_probabilities = long_df,
      median_distance = median_df,
      max_likely_distance = mode_df
    )
  }
  
  # ---- helper: chunk one displaced point --------------------------------
  process_one_point <- function(i, point_sf, feature_sf, densityBuffer,
                                point_id, feature_id, min_weight,
                                adminBound, adminID, distance_digits,
                                return_long) {
    
    pt <- point_sf[i, , drop = FALSE]
    pt <- repair_sf(pt)
    
    this_point_id <- pt[[point_id]][1]
    
    # Future me:
    # If adminBound is supplied, I trim candidate true locations to the
    # corresponding admin polygon for this displaced point.
    this_admin <- NULL
    if (!is.null(adminBound)) {
      if (is.null(adminID)) {
        stop("adminID must be supplied when adminBound is used.", call. = FALSE)
      }
      this_admin_val <- pt[[adminID]][1]
      this_admin <- adminBound[adminBound[[adminID]] == this_admin_val, , drop = FALSE]
      if (nrow(this_admin) == 0) {
        this_admin <- NULL
      } else {
        this_admin <- repair_sf(this_admin)
      }
    }
    
    # Future me:
    # pb_apply_densityBuffer() returns the gridded probability surface for this
    # displaced point after optional admin trimming.
    prob_raster <- pBuffer:::pb_apply_densityBuffer(
      point_sf = pt,
      densityBuffer = densityBuffer,
      templateRaster = NULL,
      adminBound = this_admin,
      adminID = adminID
    )
    
    if (is.null(prob_raster)) {
      return(list(
        long = NULL,
        median = data.frame(
          tmp_point = this_point_id,
          median_distance_m = NA_real_,
          stringsAsFactors = FALSE
        ),
        mode = data.frame(
          tmp_point = this_point_id,
          max_likely_distance_m = NA_real_,
          probability = NA_real_,
          stringsAsFactors = FALSE
        )
      ))
    }
    
    cand_pts <- pBuffer::pb_grid_to_points(
      prob_raster = prob_raster,
      weightsCol = "p",
      drop_zeros = TRUE
    )
    
    if (is.null(cand_pts) || nrow(cand_pts) == 0) {
      return(list(
        long = NULL,
        median = data.frame(
          tmp_point = this_point_id,
          median_distance_m = NA_real_,
          stringsAsFactors = FALSE
        ),
        mode = data.frame(
          tmp_point = this_point_id,
          max_likely_distance_m = NA_real_,
          probability = NA_real_,
          stringsAsFactors = FALSE
        )
      ))
    }
    
    cand_pts <- repair_sf(cand_pts)
    
    if (!"p" %in% names(cand_pts)) {
      stop("Candidate points are missing probability column 'p'.", call. = FALSE)
    }
    
    cand_pts <- cand_pts[cand_pts$p > min_weight, , drop = FALSE]
    cand_pts <- repair_sf(cand_pts)
    
    if (nrow(cand_pts) == 0) {
      return(list(
        long = NULL,
        median = data.frame(
          tmp_point = this_point_id,
          median_distance_m = NA_real_,
          stringsAsFactors = FALSE
        ),
        mode = data.frame(
          tmp_point = this_point_id,
          max_likely_distance_m = NA_real_,
          probability = NA_real_,
          stringsAsFactors = FALSE
        )
      ))
    }
    
    # Future me:
    # I find the nearest feature for each candidate true location, then compute
    # the nearest distance in meters for that candidate location.
    nearest_idx <- sf::st_nearest_feature(cand_pts, feature_sf)
    
    nearest_feat <- feature_sf[nearest_idx, , drop = FALSE]
    nearest_feat <- repair_sf(nearest_feat)
    
    dist_m <- as.numeric(sf::st_distance(cand_pts, nearest_feat, by_element = TRUE))
    
    long_df <- data.frame(
      tmp_point = rep(this_point_id, length(nearest_idx)),
      tmp_feature = nearest_feat[[feature_id]],
      distance_m = dist_m,
      probability = cand_pts$p,
      stringsAsFactors = FALSE
    )
    
    names(long_df)[1] <- point_id
    names(long_df)[2] <- feature_id
    
    # Future me:
    # I renormalize here in case support cells were dropped by min_weight.
    if (sum(long_df$probability, na.rm = TRUE) > 0) {
      long_df$probability <- long_df$probability / sum(long_df$probability, na.rm = TRUE)
    }
    
    median_df <- data.frame(
      tmp_point = this_point_id,
      median_distance_m = weighted_median(long_df$distance_m, long_df$probability),
      stringsAsFactors = FALSE
    )
    names(median_df)[1] <- point_id
    
    mode_dist <- long_df$distance_m
    if (!is.null(distance_digits)) {
      mode_dist <- round(mode_dist, digits = distance_digits)
    }
    
    mode_tab <- data.frame(
      distance_bucket_m = mode_dist,
      probability = long_df$probability
    ) |>
      dplyr::group_by(distance_bucket_m) |>
      dplyr::summarise(probability = sum(probability, na.rm = TRUE), .groups = "drop") |>
      dplyr::arrange(dplyr::desc(probability), distance_bucket_m)
    
    mode_df <- data.frame(
      tmp_point = this_point_id,
      max_likely_distance_m = mode_tab$distance_bucket_m[1],
      probability = mode_tab$probability[1],
      stringsAsFactors = FALSE
    )
    names(mode_df)[1] <- point_id
    
    if (!isTRUE(return_long)) {
      long_df <- NULL
    } else {
      # Future me:
      # Aggregate repeated feature-distance combinations so the long table is
      # much smaller than one row per support cell.
      long_df <- long_df |>
        dplyr::group_by(.data[[point_id]], .data[[feature_id]], distance_m) |>
        dplyr::summarise(probability = sum(probability, na.rm = TRUE), .groups = "drop")
      
      if (!is.null(distance_digits)) {
        long_df$distance_m <- round(long_df$distance_m, digits = distance_digits)
        long_df <- long_df |>
          dplyr::group_by(.data[[point_id]], .data[[feature_id]], distance_m) |>
          dplyr::summarise(probability = sum(probability, na.rm = TRUE), .groups = "drop")
      }
    }
    
    list(
      long = long_df,
      median = median_df,
      mode = mode_df
    )
  }
  
  # ---- helper: write one chunk of long output ---------------------------
  append_long_csv <- function(df, file_path, append = TRUE) {
    if (is.null(df) || nrow(df) == 0) {
      return(invisible(NULL))
    }
    
    dir.create(dirname(file_path), recursive = TRUE, showWarnings = FALSE)
    
    utils::write.table(
      df,
      file = file_path,
      sep = ",",
      row.names = FALSE,
      col.names = !append || !file.exists(file_path),
      append = append && file.exists(file_path),
      qmethod = "double"
    )
    
    invisible(NULL)
  }
  
  # ---- validate and repair inputs ---------------------------------------
  point_sf <- repair_sf(point_sf)
  feature_sf <- repair_sf(feature_sf)
  
  if (!is.null(adminBound)) {
    adminBound <- repair_sf(adminBound)
  }
  
  if (!point_id %in% names(point_sf)) {
    stop(point_id, " is missing from point_sf.", call. = FALSE)
  }
  
  if (!feature_id %in% names(feature_sf)) {
    stop(feature_id, " is missing from feature_sf.", call. = FALSE)
  }
  
  if (anyDuplicated(point_sf[[point_id]]) > 0) {
    stop(point_id, " must uniquely identify rows in point_sf.", call. = FALSE)
  }
  
  if (anyDuplicated(feature_sf[[feature_id]]) > 0) {
    stop(feature_id, " must uniquely identify rows in feature_sf.", call. = FALSE)
  }
  
  if (!is.null(adminBound)) {
    if (is.null(adminID)) {
      stop("adminID must be supplied when adminBound is used.", call. = FALSE)
    }
    if (!adminID %in% names(adminBound)) {
      stop(adminID, " is missing from adminBound.", call. = FALSE)
    }
    if (!adminID %in% names(point_sf)) {
      stop(adminID, " is missing from point_sf.", call. = FALSE)
    }
    adminBound <- keep_sf_cols(adminBound, adminID)
  }
  
  point_keep <- c(point_id, if (!is.null(adminBound)) adminID else NULL)
  point_sf <- keep_sf_cols(point_sf, point_keep)
  feature_sf <- keep_sf_cols(feature_sf, feature_id)
  
  if (!is.null(long_file) && file.exists(long_file)) {
    if (isTRUE(overwrite_long_file)) {
      file.remove(long_file)
    } else {
      stop("long_file already exists; set overwrite_long_file = TRUE.", call. = FALSE)
    }
  }
  
  if (nrow(point_sf) == 0) {
    return(empty_output(
      return_long = return_long,
      point_id = point_id,
      feature_id = feature_id
    ))
  }
  
  # ---- build index chunks, not sf-object chunks -------------------------
  idx <- seq_len(nrow(point_sf))
  chunk_ids <- ceiling(idx / chunk_size)
  idx_chunks <- split(idx, chunk_ids)
  
  # Future me:
  # This is the key robustness change: I chunk by row indices, then subset
  # inside the worker. I do not split the sf object itself and ship those
  # split objects around.
  process_chunk <- function(idx_chunk) {
    chunk_res <- lapply(
      idx_chunk,
      FUN = process_one_point,
      point_sf = point_sf,
      feature_sf = feature_sf,
      densityBuffer = densityBuffer,
      point_id = point_id,
      feature_id = feature_id,
      min_weight = min_weight,
      adminBound = adminBound,
      adminID = adminID,
      distance_digits = distance_digits,
      return_long = return_long
    )
    
    list(
      long = dplyr::bind_rows(lapply(chunk_res, `[[`, "long")),
      median = dplyr::bind_rows(lapply(chunk_res, `[[`, "median")),
      mode = dplyr::bind_rows(lapply(chunk_res, `[[`, "mode"))
    )
  }
  
  if (isTRUE(parallel) && n.cores > 1) {
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    
    future::plan(future::multisession, workers = n.cores)
    
    chunk_out <- future.apply::future_lapply(
      X = idx_chunks,
      FUN = process_chunk
    )
  } else {
    chunk_out <- lapply(
      X = idx_chunks,
      FUN = process_chunk
    )
  }
  
  if (!is.null(long_file)) {
    for (j in seq_along(chunk_out)) {
      append_long_csv(
        df = chunk_out[[j]]$long,
        file_path = long_file,
        append = TRUE
      )
      chunk_out[[j]]$long <- NULL
    }
  }
  
  out <- list(
    nearest_distance_probabilities = if (isTRUE(return_long) && is.null(long_file)) {
      dplyr::bind_rows(lapply(chunk_out, `[[`, "long"))
    } else {
      NULL
    },
    median_distance = dplyr::bind_rows(lapply(chunk_out, `[[`, "median")),
    max_likely_distance = dplyr::bind_rows(lapply(chunk_out, `[[`, "mode"))
  )
  
  out
}
