#' @name pb_radiusLimitJoin
#' @rdname pb_radiusLimitJoin
#' @title  pb_radiusLimitJoin
#'
#' @description For each displaced point, turns a density buffer into a weighted
#' set of plausible true locations, counts nearby point features within a fixed
#' radius at each candidate location, and then summarizes the resulting
#' probability distribution.
#' 
#' @returns A list with three data frames:
#' \describe{
#'   \item{feature_probabilities}{Columns for displaced point ID, joined point
#'   ID, and the probability that the joined point is within \code{radiusLength}
#'   of the true location.}
#'   \item{count_distribution}{Columns for displaced point ID, the number of
#'   joined points within \code{radiusLength}, and the probability of that count.}
#'   \item{median_count}{Columns for displaced point ID and the weighted median
#'   count across all candidate true locations.}
#' }
#'
#'
#' @param displaced.sf Displaced coordinates in an sf object type.
#' @param features2count.sf sf point feature that is counted.
#' @param limitDist a numeric vector to determine the distance from a displaced community within which the number of features2count.sf will be counted.
#' @param displaced.id a character vector of the column name specifying the unique ID of the displaced communities
#' @param features2count.id a charactor vector of the column name specifying the unique ID of the point features being connected
#' @param metrics Character vector of column names containing the numeric metrics you wish to calculate.
#' @param densityBuffer A density buffer list object created with the pb_integratedDens function.
#' @param adminBound (Optional) the administrative boundary that circumscribes the displacement.
#' @param adminID The unique ID for the adminBound features.
#' @param weightedMedian Boolean value, if true, the weighted median will be calculated. If false, the weighted mean will be calculated. Defaults to true.
#' @param probCuts A numeric vector containing values  between 0 and 1, these are percentiles to report the credible interval of distances
#' @param n.cores This allows for parallelization using the futures package. This can be slightly unstable, but normally functions well and dramatically speeds compute time when it does. n.cores specifies the number of cores to use. The default is 1, and does not parallelize.
#'
#' @author J.W. Rozelle
#'
#'
#' @export pb_radiusLimitJoin
#' @examples
#'
#' # coming soon!
#' 
#' 
#' 
#' 


# polygons.sf
# features2count.sf <- htSPAge_t2.sf
# n.cores <- 18
# densityBuffer <- pb_integratedDensity(boundaries = 5e3)


pb_radiusLimitJoin <- function(displaced.sf, 
                               features2count.sf, 
                               limitDist = 5e3,
                               displaced.id = "DHSID", 
                               features2count.id = "SPAID",
                               metrics = c("sri_score", 
                                           "sri_basicamenities", 
                                           "sri_basicequip", 
                                           "sri_diagcapacity", 
                                           "sri_infprev", 
                                           "sri_med"),
                               densityBuffer = NULL, 
                               adminBound = NULL, 
                               adminID = NULL,
                               weightedMedian = TRUE,
                               probCuts = c(0.1, 0.9),
                               n.cores = 1) {
  
  # --- validation (centralized) ---
  pb_check_sf(displaced.sf, "displaced.sf")
  pb_check_geom_type(displaced.sf, "POINT", "displaced.sf")
  pb_check_crs(displaced.sf, "displaced.sf")
  pb_check_projected(displaced.sf, "displaced.sf")
  
  pb_check_sf(features2count.sf, "features2count.sf")
  pb_check_crs(features2count.sf, "features2count.sf")
  
  pb_check_cols(displaced.sf, displaced.id, "displaced.sf")
  pb_check_unique_id(sf::st_drop_geometry(displaced.sf), displaced.id, "displaced.sf")
  
  if (!is.null(adminBound)) {
    pb_check_sf(adminBound, "adminBound")
    if (is.null(adminID)) stop("adminID must be provided when adminBound is provided.", call. = FALSE)
    pb_check_cols(adminBound, adminID, "adminBound")
    pb_check_unique_id(sf::st_drop_geometry(adminBound), adminID, "adminBound")
  }
  
  if (is.null(densityBuffer)) {
    stop("densityBuffer must be provided (e.g., pb_integratedDensity(...) or pb_mcDensity(...)).", call. = FALSE)
  }
  
  pb_check_pos_scalar(limitDist, "limitDist")
  pb_check_prob_vec(probCuts, "probCuts")
  
  if (!all(metrics %in% names(features2count.sf))) {
    stop(
      "One or more column names specified in metrics cannot be found in features2count.sf: ",
      paste(setdiff(metrics, names(features2count.sf)), collapse = ", "),
      call. = FALSE
    )
  }
  
  # standard internal column names
  displaced.sf$DHSID <- displaced.sf[[displaced.id]]
  features2count.sf$SPAID_use <- features2count.sf[[features2count.id]]
  
  if (!is.null(adminBound)) {
    adminBound$adminID <- sf::st_drop_geometry(adminBound)[[adminID]]
  }
  
  # worker function (used by both branches)
  one_comm <- function(xrow) {
    
    rowDHSID <- xrow[[displaced.id]]
    singleComm <- displaced.sf[displaced.sf[[displaced.id]] == rowDHSID, ]
    
    # guard: missing/zero geometry (your convention)
    xy <- sf::st_coordinates(singleComm)
    if (nrow(xy) < 1 || is.na(xy[1,2]) || xy[1,2] == 0) return(NA)
    
    crs2use <- raster::crs(singleComm)
    
    # 1) build centered probability raster
    singleDens.raster <- rasterizeDisplacement(
      densityBuffer,
      initialCoords = xy[1, 1:2],
      inputCRS = crs2use
    )
    
    # 2) optionally admin-trim + renormalize
    if (!is.null(adminBound)) {
      singleDens.raster <- pb_trim_probRaster_to_point_admin(
        probRaster = singleDens.raster,
        point_sf   = singleComm,
        adminBound = adminBound,
        adminID    = "adminID"
      )
    }
    
    # 3) raster -> weighted points
    singleDens.sf <- st_rasterAsPoint(singleDens.raster)
    singleDens.sf$dispid <- row.names(singleDens.sf)
    rm(singleDens.raster)
    
    # 4) distances from weighted points to all facilities
    distance.mtx <- sf::st_distance(singleDens.sf, features2count.sf) |> as.matrix()
    
    row.names(distance.mtx) <- singleDens.sf$dispid
    colnames(distance.mtx) <- features2count.sf$SPAID_use
    
    distance.df <- as.data.frame(distance.mtx)
    distance.df$DHSID <- rowDHSID
    distance.df$dispid <- row.names(distance.mtx)
    rm(distance.mtx)
    
    singleDens.sf <- merge(singleDens.sf, distance.df, by = "dispid")
    singleDens.df <- sf::st_drop_geometry(singleDens.sf)
    rm(singleDens.sf)
    
    # 5) summarize per facility: median distance + prob within limit + CrI
    featureDistances.list <- lapply(features2count.sf$SPAID_use, function(fid) {
      
      singleFeatureDistance.df <- singleDens.df[c(fid, "layer")]
      singleFeatureDistance.df$distance <- as.numeric(singleFeatureDistance.df[[fid]])
      
      if (isTRUE(weightedMedian)) {
        likelyDist <- spatstat.geom::weighted.median(singleFeatureDistance.df$distance,
                                                     singleFeatureDistance.df$layer)
      } else {
        # weighted mean
        likelyDist <- stats::weighted.mean(singleFeatureDistance.df$distance,
                                           singleFeatureDistance.df$layer, na.rm = TRUE)
      }
      
      distCrI <- spatstat.geom::weighted.quantile(
        singleFeatureDistance.df$distance,
        singleFeatureDistance.df$layer,
        probs = probCuts
      )
      names(distCrI) <- paste0("limit_", names(distCrI))
      distCrI.df <- as.data.frame(t(distCrI))
      
      probInLimit <- sum(
        ifelse(singleFeatureDistance.df$distance <= limitDist &
                 !is.na(singleFeatureDistance.df$distance) &
                 !is.na(singleFeatureDistance.df$layer),
               singleFeatureDistance.df$layer,
               0),
        na.rm = TRUE
      )
      
      out <- data.frame(likelyDist = likelyDist, probInLimit = probInLimit)
      out[[features2count.id]] <- fid
      cbind(out, distCrI.df)
    })
    
    hf_weightedDist.df <- iotools::fdrbind(featureDistances.list)
    hf_weightedDist.df <- dplyr::filter(hf_weightedDist.df, probInLimit > 0)
    
    if (!is.null(metrics)) {
      hf_weightedDist.df <- merge(
        hf_weightedDist.df,
        sf::st_drop_geometry(features2count.sf)[c(features2count.id, metrics)],
        by = features2count.id,
        all.x = TRUE, all.y = FALSE
      )
    }
    
    hf_weightedDist.df
  }
  
  # run
  if (n.cores > 1) {
    cl <- parallel::makeCluster(n.cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    future::plan(future::cluster, workers = cl)
    
    numFeatures <- future.apply::future_apply(sf::st_drop_geometry(displaced.sf), 1, one_comm)
  } else {
    numFeatures <- apply(sf::st_drop_geometry(displaced.sf), 1, one_comm)
  }
  
  iotools::fdrbind(numFeatures)
}






#' @name pb_countFeatures
#' @rdname pb_countFeatures
#' @title pb_countFeatures
#'
#' @description
#' For each displaced point, this function builds the set of plausible true
#' locations implied by \code{densityBuffer}, keeps those locations as an
#' \code{sf} point layer, and then counts how many joined point features fall
#' within \code{radiusLength} of each plausible true location.
#'
#' In other words, this is a point-to-point join. I am not extracting raster
#' values here, and I am not summarizing a raster surface as the substantive
#' output. The only thing the function returns is the probability-weighted join
#' result for point or multipoint features.
#'
#' The function returns three related summaries for each displaced point:
#' \enumerate{
#'   \item the probability that each joined feature falls within the requested
#'   radius of the true location;
#'   \item the probability distribution over the number of joined features that
#'   fall within that radius; and
#'   \item the weighted median of that count distribution.
#' }
#'
#' @returns A list with three data frames:
#' \describe{
#'   \item{feature_probabilities}{A data frame with three columns: displaced-point
#'   ID, joined-point ID, and the probability that the joined feature is within
#'   \code{radiusLength} of the true location. If no joined feature has positive
#'   probability mass for a displaced point, that displaced point simply
#'   contributes zero rows to this table.}
#'   \item{count_distribution}{A data frame with three columns: displaced-point ID,
#'   count of joined features within \code{radiusLength}, and the probability of
#'   observing that count. The probabilities sum to 1 within each displaced
#'   point after any trimming or weight thresholding.}
#'   \item{median_count}{A data frame with two columns: displaced-point ID and the
#'   weighted median count of joined features within \code{radiusLength}.}
#' }
#'
#' @param displaced.sf An \code{sf} object with POINT geometries for displaced
#' cluster locations.
#' @param features2count.sf An \code{sf} object with POINT or MULTIPOINT
#' geometries that should be counted within \code{radiusLength}.
#' @param radiusLength A single numeric value giving the search radius, in the
#' linear units of the projected CRS (normally meters).
#' @param displaced.id Column name giving the unique displaced-point ID.
#' @param features2count.id Column name giving the unique joined-feature ID.
#' @param densityBuffer A density-buffer object created by
#' \code{pb_densityBuffer()}, \code{pb_integratedDensity()},
#' \code{pb_Density()}, \code{pb_mcDensity()}, or
#' \code{pb_mcDensity_dhs()}.
#' @param templateRaster Kept only for backward compatibility. It is not used by
#' this implementation because the join is carried out directly on candidate
#' \code{sf} points rather than through an intermediate raster object.
#' @param adminBound Optional administrative boundary used to trim plausible true
#' locations before counting. If multiple admin polygons are supplied, the
#' displaced point is first matched to its containing admin polygon and only that
#' polygon is used for trimming.
#' @param adminID Column name giving the admin ID in \code{adminBound}. This is
#' only required when \code{adminBound} has more than one row.
#' @param min_weight Drop candidate true-location support points with weights less
#' than or equal to this threshold before counting, then renormalize the
#' remaining weights so the resulting probabilities still sum to 1.
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
#' @export pb_countFeatures
#' @examples
#' # coming soon!
pb_countFeatures <- function(displaced.sf,
                             features2count.sf,
                             radiusLength = 1e4,
                             displaced.id = "DHSID",
                             features2count.id = "SPAID",
                             densityBuffer = NULL,
                             templateRaster = NULL,
                             adminBound = NULL,
                             adminID = NULL,
                             min_weight = 0,
                             n.cores = 1,
                             parallel = FALSE,
                             future.seed = TRUE,
                             future.scheduling = 1) {
  
  # Future-me: this join only makes sense for displaced POINT locations.
  # I am being strict here on purpose so I do not accidentally let one DHSID
  # sneak in with multiple source coordinates and then silently use the first.
  pb_check_sf(displaced.sf, "displaced.sf")
  pb_check_geom_type(displaced.sf, "POINT", "displaced.sf")
  pb_check_crs(displaced.sf, "displaced.sf")
  pb_check_projected(displaced.sf, "displaced.sf")
  pb_check_cols(displaced.sf, displaced.id, "displaced.sf")
  pb_check_unique_id(sf::st_drop_geometry(displaced.sf), displaced.id, "displaced.sf")
  
  # Future-me: the counted layer really is the thing this function is about.
  # I only want point-like features here because I am counting whether a feature
  # is inside a radius, not measuring overlap with polygons or extracting
  # anything from rasters.
  pb_check_sf(features2count.sf, "features2count.sf")
  pb_check_crs(features2count.sf, "features2count.sf")
  pb_check_cols(features2count.sf, features2count.id, "features2count.sf")
  pb_check_unique_id(sf::st_drop_geometry(features2count.sf), features2count.id, "features2count.sf")
  
  feature_geom <- unique(as.character(sf::st_geometry_type(features2count.sf)))
  if (!all(feature_geom %in% c("POINT", "MULTIPOINT"))) {
    stop(
      "features2count.sf must have POINT or MULTIPOINT geometries.",
      call. = FALSE
    )
  }
  
  pb_check_pos_scalar(radiusLength, "radiusLength")
  pb_check_pos_scalar(n.cores, "n.cores")
  if (!is.numeric(min_weight) || length(min_weight) != 1 || is.na(min_weight) || min_weight < 0) {
    stop("min_weight must be a single non-negative number.", call. = FALSE)
  }
  
  # Future-me: I am keeping templateRaster in the signature so old code does not
  # break, but I am intentionally not using it anymore.
  
  if (is.null(densityBuffer)) {
    stop(
      "densityBuffer must be provided (for example from pb_densityBuffer(), ",
      "pb_integratedDensity(), pb_Density(), pb_mcDensity(), or pb_mcDensity_dhs()).",
      call. = FALSE
    )
  }
  
  # Future-me: the newer package code uses class 'pb_densityBuffer', but a lot
  # of older helpers in this same source tree still return the old-style list
  # with weightedCircle/cellMeters/radiusMeters. I want pb_countFeatures to work
  # with both so it is actually drop-in safe inside this package.
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
  
  # Future-me: the join has to happen in one projected CRS. If the counted layer
  # or the admin boundary is in something else, I just move them onto the CRS of
  # the displaced points before I do any distance work.
  target_crs <- sf::st_crs(displaced.sf)
  if (sf::st_crs(features2count.sf) != target_crs) {
    features2count.sf <- sf::st_transform(features2count.sf, target_crs)
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
  
  features2count_ids <- features2count.sf[[features2count.id]]
  
  make_feature_prob_na <- function(id_value) {
    out <- data.frame(
      tmp_displaced_id = id_value,
      tmp_feature_id = features2count_ids[NA_integer_][1],
      probability = NA_real_,
      stringsAsFactors = FALSE
    )
    names(out)[1] <- displaced.id
    names(out)[2] <- features2count.id
    out
  }
  
  make_feature_prob_empty <- function(id_value) {
    out <- data.frame(
      tmp_displaced_id = id_value[FALSE],
      tmp_feature_id = features2count_ids[FALSE],
      probability = numeric(0),
      stringsAsFactors = FALSE
    )
    names(out)[1] <- displaced.id
    names(out)[2] <- features2count.id
    out
  }
  
  make_count_na <- function(id_value) {
    out <- data.frame(
      tmp_displaced_id = id_value,
      point_count = NA_integer_,
      probability = NA_real_,
      stringsAsFactors = FALSE
    )
    names(out)[1] <- displaced.id
    out
  }
  
  make_median_na <- function(id_value) {
    out <- data.frame(
      tmp_displaced_id = id_value,
      median_count = NA_real_,
      stringsAsFactors = FALSE
    )
    names(out)[1] <- displaced.id
    out
  }
  
  matrix_to_candidate_points <- function(weight_matrix, cell_m, x0, y0, weight_scale = 1) {
    # Future-me: this recreates the same cell centers I would have gotten from
    # the raster path, but keeps everything as sf points from start to finish.
    # Row 1 of the matrix is the top of the support, so y decreases with row.
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
  
  choose_density_components <- function(singleComm) {
    # Future-me: new-style density objects can store several kernels plus an
    # optional mixture. I am reconstructing the support points directly from the
    # stored matrices here instead of going through an intermediate raster.
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
    if (urban_col %in% names(singleComm)) {
      urban_raw <- as.character(sf::st_drop_geometry(singleComm)[[urban_col]][1])
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
    
    # Future-me: if a rural-only mixture object also carries a dedicated urban
    # radius, it should be in radii_m but outside mixture_radii_m. I use that
    # first. If it is missing, I fall back to the smallest stored radius instead
    # of the largest one so I do not accidentally balloon an urban kernel to a
    # rural maximum.
    non_mix_radii <- setdiff(densityBuffer$radii_m, mix$radii_m)
    fallback_radius <- if (length(non_mix_radii) > 0) {
      max(non_mix_radii)
    } else {
      min(densityBuffer$radii_m)
    }
    
    list(list(weight = 1, radius_m = fallback_radius))
  }
  
  get_trim_polygon <- function(singleComm) {
    if (is.null(adminBound)) {
      return(NULL)
    }
    
    if (nrow(adminBound) == 1) {
      return(adminBound)
    }
    
    point_min <- singleComm[, character(0), drop = FALSE]
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
  
  build_candidate_support <- function(singleComm) {
    xy <- suppressWarnings(sf::st_coordinates(singleComm))
    if (nrow(xy) != 1 || anyNA(xy[1, c("X", "Y")])) {
      return(NULL)
    }
    
    x0 <- xy[1, "X"]
    y0 <- xy[1, "Y"]
    
    components <- choose_density_components(singleComm)
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
    
    # Future-me: mixture components can land on the same support point. I need
    # to add those weights before I do any trimming or thresholding.
    candidate_df <- stats::aggregate(
      candidate_df$weight,
      by = list(x = candidate_df$x, y = candidate_df$y),
      FUN = sum
    )
    names(candidate_df) <- c("x", "y", "weight")
    
    candidate_sf <- sf::st_as_sf(candidate_df, coords = c("x", "y"), crs = sf::st_crs(singleComm))
    
    trim_poly <- get_trim_polygon(singleComm)
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
    
    # Future-me: once I have trimmed the support and possibly dropped tiny
    # weights, I have to renormalize so the outputs are still real probability
    # distributions that sum to one.
    candidate_sf$weight <- candidate_sf$weight / sum(candidate_sf$weight)
    candidate_sf
  }
  
  one_comm <- function(row_index) {
    rowDHSID <- displaced.sf[[displaced.id]][row_index]
    singleComm <- displaced.sf[row_index, , drop = FALSE]
    
    # Future-me: keep the bad-geometry escape hatch explicit. If a displaced
    # point is empty or malformed, I want a marked missing result, not a partial
    # object that looks real later on.
    if (isTRUE(sf::st_is_empty(singleComm))) {
      return(list(
        feature_probabilities = make_feature_prob_na(rowDHSID),
        count_distribution = make_count_na(rowDHSID),
        median_count = make_median_na(rowDHSID)
      ))
    }
    
    candidate_sf <- build_candidate_support(singleComm)
    
    if (is.null(candidate_sf) || nrow(candidate_sf) == 0) {
      return(list(
        feature_probabilities = make_feature_prob_na(rowDHSID),
        count_distribution = make_count_na(rowDHSID),
        median_count = make_median_na(rowDHSID)
      ))
    }
    
    # Future-me: this is the actual join. Every candidate true point asks which
    # joined features are inside the requested radius. After that, the summaries
    # are just weighted aggregations.
    within_list <- sf::st_is_within_distance(
      candidate_sf,
      features2count.sf,
      dist = radiusLength
    )
    
    candidate_df <- sf::st_drop_geometry(candidate_sf)
    point_counts <- lengths(within_list)
    
    feature_hit_list <- lapply(seq_along(within_list), function(i) {
      hit_idx <- within_list[[i]]
      if (length(hit_idx) == 0) {
        return(NULL)
      }
      
      data.frame(
        tmp_feature_id = features2count_ids[hit_idx],
        probability = candidate_df$weight[i],
        stringsAsFactors = FALSE
      )
    })
    feature_hit_list <- Filter(Negate(is.null), feature_hit_list)
    
    if (length(feature_hit_list) == 0) {
      feature_prob_df <- make_feature_prob_empty(rowDHSID)
    } else {
      feature_prob_df <- do.call(rbind, feature_hit_list)
      feature_prob_df <- stats::aggregate(
        feature_prob_df$probability,
        by = list(tmp_feature_id = feature_prob_df$tmp_feature_id),
        FUN = sum
      )
      names(feature_prob_df) <- c(features2count.id, "probability")
      feature_prob_df[[displaced.id]] <- rowDHSID
      feature_prob_df <- feature_prob_df[, c(displaced.id, features2count.id, "probability"), drop = FALSE]
      feature_prob_df <- feature_prob_df[
        order(feature_prob_df$probability, decreasing = TRUE),
        ,
        drop = FALSE
      ]
    }
    
    count_distribution <- stats::aggregate(
      candidate_df$weight,
      by = list(point_count = point_counts),
      FUN = sum
    )
    names(count_distribution) <- c("point_count", "probability")
    count_distribution[[displaced.id]] <- rowDHSID
    count_distribution <- count_distribution[, c(displaced.id, "point_count", "probability"), drop = FALSE]
    count_distribution$point_count <- as.integer(count_distribution$point_count)
    count_distribution <- count_distribution[
      order(count_distribution$point_count),
      ,
      drop = FALSE
    ]
    
    median_count <- data.frame(
      tmp_displaced_id = rowDHSID,
      median_count = as.numeric(
        pb__weighted_quantile(
          x = point_counts,
          w = candidate_df$weight,
          probs = 0.5
        )[[1]]
      ),
      stringsAsFactors = FALSE
    )
    names(median_count)[1] <- displaced.id
    
    list(
      feature_probabilities = feature_prob_df,
      count_distribution = count_distribution,
      median_count = median_count
    )
  }
  
  row_index <- seq_len(nrow(displaced.sf))
  use_parallel <- isTRUE(parallel) || n.cores > 1
  
  if (use_parallel) {
    # Future-me: this is embarrassingly parallel across displaced points, so I
    # let future.apply handle the split. I also restore the incoming plan on the
    # way out so this function does not leave side effects in the caller's
    # session.
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    future::plan(future::multisession, workers = n.cores)
    
    count_results <- future.apply::future_lapply(
      X = row_index,
      FUN = one_comm,
      future.seed = future.seed,
      future.scheduling = future.scheduling
    )
  } else {
    count_results <- lapply(row_index, one_comm)
  }
  
  feature_probabilities <- do.call(
    rbind,
    lapply(count_results, `[[`, "feature_probabilities")
  )
  count_distribution <- do.call(
    rbind,
    lapply(count_results, `[[`, "count_distribution")
  )
  median_count <- do.call(
    rbind,
    lapply(count_results, `[[`, "median_count")
  )
  
  rownames(feature_probabilities) <- NULL
  rownames(count_distribution) <- NULL
  rownames(median_count) <- NULL
  
  list(
    feature_probabilities = feature_probabilities,
    count_distribution = count_distribution,
    median_count = median_count
  )
}
