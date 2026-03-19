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
#' @title  pb_countFeatures
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
#' @param radiusLength a numeric vector to determine the distance from a displaced community within which the number of features2count.sf will be counted.
#' @param displaced.id Column name giving the unique displaced-point ID.
#' @param features2count.id Column name giving the unique joined-point ID.
#' @param densityBuffer A density buffer list object created with the pb_integratedDens function.
#' @param templateRaster Optional RasterLayer used to supply CRS/alignment when
#' applying the density buffer.
#' @param adminBound (Optional) the administrative boundary that circumscribes the displacement.
#' @param adminID The unique ID for the adminBound features, defaults to "ID_2".
#' @param min_weight Drop candidate true-location cells with weights less than or
#' equal to this threshold before counting.
#' @param n.cores Number of workers to use when parallelization is enabled.
#' @param parallel Logical; if TRUE, run displaced-point jobs with
#' \code{future.apply}. If FALSE, run serially unless \code{n.cores > 1}.
#' @param future.seed Seed handling passed to \code{future_lapply()} for
#' reproducible parallel execution. Defaults to TRUE.
#' @param future.scheduling Scheduling value passed to
#' \code{future.apply::future_lapply()}.
#' 
#' @author J.W. Rozelle
#'
#'
#' @export pb_countFeatures
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
  
  pb_check_sf(displaced.sf, "displaced.sf")
  pb_check_geom_type(displaced.sf, "POINT", "displaced.sf")
  pb_check_crs(displaced.sf, "displaced.sf")
  pb_check_projected(displaced.sf, "displaced.sf")
  
  pb_check_sf(features2count.sf, "features2count.sf")
  pb_check_crs(features2count.sf, "features2count.sf")
  
  pb_check_cols(displaced.sf, displaced.id, "displaced.sf")
  pb_check_cols(features2count.sf, features2count.id, "features2count.sf")
  pb_check_unique_id(sf::st_drop_geometry(displaced.sf), displaced.id, "displaced.sf")
  pb_check_unique_id(sf::st_drop_geometry(features2count.sf), features2count.id, "features2count.sf")
  
  pb_check_pos_scalar(radiusLength, "radiusLength")
  pb_check_pos_scalar(n.cores, "n.cores")
  
  if (is.null(densityBuffer)) {
    stop("densityBuffer must be provided (e.g., pb_integratedDensity(...) or pb_mcDensity(...)).", call. = FALSE)
  }
  if (!inherits(densityBuffer, "pb_densityBuffer")) {
    stop("densityBuffer must be a pb_densityBuffer object.", call. = FALSE)
  }
  if (is.null(templateRaster)) {
    templateRaster <- raster::raster(
      raster::extent(0, densityBuffer$cell_m, 0, densityBuffer$cell_m),
      res = densityBuffer$cell_m,
      crs = sf::st_crs(displaced.sf)$wkt
    )
  }
  if (!inherits(templateRaster, "RasterLayer")) {
    stop("templateRaster must be a raster::RasterLayer.", call. = FALSE)
  }
  
  if (!is.null(adminBound)) {
    pb_check_sf(adminBound, "adminBound")
    if (is.null(adminID)) stop("adminID must be provided when adminBound is provided.", call. = FALSE)
    pb_check_cols(adminBound, adminID, "adminBound")
    pb_check_unique_id(sf::st_drop_geometry(adminBound), adminID, "adminBound")
  }
  
  if (!is.null(adminBound)) {
    adminBound$adminID <- sf::st_drop_geometry(adminBound)[[adminID]]
  }
  
  one_comm <- function(row_index) {
    
    rowDHSID <- displaced.sf[[displaced.id]][row_index]
    singleComm <- displaced.sf[displaced.sf[[displaced.id]] == rowDHSID, ]
    
    # Future-me: keep the "bad geometry" escape hatch boring and explicit.
    # If the displaced point is missing, I want the downstream shapes to stay
    # predictable instead of exploding halfway through a batch.
    xy <- sf::st_coordinates(singleComm)
    if (nrow(xy) < 1 || anyNA(xy[1, c("X", "Y")])) {
      feature_probs <- data.frame(
        displaced_id_value = rowDHSID,
        feature_id_value = NA,
        probability = NA_real_
      )
      count_distribution <- data.frame(
        displaced_id_value = rowDHSID,
        point_count = NA_integer_,
        probability = NA_real_
      )
      median_count <- data.frame(
        displaced_id_value = rowDHSID,
        median_count = NA_real_
      )
      names(feature_probs)[1] <- displaced.id
      names(feature_probs)[2] <- features2count.id
      names(count_distribution)[1] <- displaced.id
      names(median_count)[1] <- displaced.id
      return(list(
        feature_probabilities = feature_probs,
        count_distribution = count_distribution,
        median_count = median_count
      ))
    }
    
    # Future-me: always go through pb_apply_densityBuffer() here so I get the
    # same behavior for integrated kernels and gridded MC kernels. If the MC
    # object was bucketed to save space, each grid point already carries the
    # probability mass I need.
    singleDens.raster <- pb_apply_densityBuffer(
      point_sf = singleComm,
      densityBuffer = densityBuffer,
      templateRaster = templateRaster,
      adminBound = adminBound,
      adminID = "adminID"
    )
    
    if (is.null(singleDens.raster)) {
      feature_probs <- data.frame(
        displaced_id_value = rowDHSID,
        feature_id_value = NA,
        probability = NA_real_
      )
      count_distribution <- data.frame(
        displaced_id_value = rowDHSID,
        point_count = NA_integer_,
        probability = NA_real_
      )
      median_count <- data.frame(
        displaced_id_value = rowDHSID,
        median_count = NA_real_
      )
      names(feature_probs)[1] <- displaced.id
      names(feature_probs)[2] <- features2count.id
      names(count_distribution)[1] <- displaced.id
      names(median_count)[1] <- displaced.id
      return(list(
        feature_probabilities = feature_probs,
        count_distribution = count_distribution,
        median_count = median_count
      ))
    }
    
    # Future-me: convert the probability raster to candidate true-point
    # locations. For gridded MC objects, this is the "bucketed simulation"
    # representation, so every point already has an associated probability.
    singleDens.sf <- pb_grid_to_points(
      prob_raster = singleDens.raster,
      weightsCol = "layer",
      drop_zeros = TRUE
    )
    singleDens.sf <- singleDens.sf[singleDens.sf$layer > min_weight, , drop = FALSE]
    
    if (nrow(singleDens.sf) == 0) {
      feature_probs <- data.frame(
        displaced_id_value = rowDHSID,
        feature_id_value = NA,
        probability = NA_real_
      )
      count_distribution <- data.frame(
        displaced_id_value = rowDHSID,
        point_count = NA_integer_,
        probability = NA_real_
      )
      median_count <- data.frame(
        displaced_id_value = rowDHSID,
        median_count = NA_real_
      )
      names(feature_probs)[1] <- displaced.id
      names(feature_probs)[2] <- features2count.id
      names(count_distribution)[1] <- displaced.id
      names(median_count)[1] <- displaced.id
      return(list(
        feature_probabilities = feature_probs,
        count_distribution = count_distribution,
        median_count = median_count
      ))
    }
    
    # Future-me: this is the core join. For each plausible true location, find
    # every feature that falls inside the requested radius. Because the
    # candidate locations are probability-weighted, all later summaries can just
    # sum those weights.
    within_list <- sf::st_is_within_distance(singleDens.sf, features2count.sf, dist = radiusLength)
    singleDens.df <- sf::st_drop_geometry(singleDens.sf)
    
    # Future-me: first output = per feature probability of being inside the
    # radius from the true point. If a feature shows up for multiple candidate
    # locations, I add those candidate probabilities together.
    feature_prob_list <- lapply(seq_along(within_list), function(idx) {
      hit_idx <- within_list[[idx]]
      if (length(hit_idx) == 0) return(NULL)
      data.frame(
        feature_id_value = features2count.sf[[features2count.id]][hit_idx],
        probability = singleDens.df$layer[idx],
        stringsAsFactors = FALSE
      )
    })
    feature_prob_df <- iotools::fdrbind(feature_prob_list)
    
    if (is.null(feature_prob_df) || nrow(feature_prob_df) == 0) {
      feature_prob_df <- data.frame(
        displaced_id_value = rowDHSID,
        feature_id_value = NA,
        probability = 0
      )
      names(feature_prob_df)[1] <- displaced.id
      names(feature_prob_df)[2] <- features2count.id
    } else {
      feature_prob_df <- stats::aggregate(
        feature_prob_df$probability,
        by = list(feature_id_value = feature_prob_df$feature_id_value),
        FUN = sum
      )
      names(feature_prob_df) <- c(features2count.id, "probability")
      feature_prob_df[[displaced.id]] <- rowDHSID
      feature_prob_df <- feature_prob_df[, c(displaced.id, features2count.id, "probability")]
    }
    
    # Future-me: second output = probability distribution over counts. Count how
    # many features each candidate true location can reach, then collapse those
    # counts by probability mass.
    singleDens.df$point_count <- lengths(within_list)
    count_distribution <- stats::aggregate(
      singleDens.df$layer,
      by = list(point_count = singleDens.df$point_count),
      FUN = sum
    )
    names(count_distribution) <- c("point_count", "probability")
    count_distribution[[displaced.id]] <- rowDHSID
    count_distribution <- count_distribution[, c(displaced.id, "point_count", "probability")]
    count_distribution <- count_distribution[order(count_distribution$point_count), , drop = FALSE]
    
    # Future-me: third output = weighted median of the count distribution. This
    # is my "most likely middle" count across all plausible true locations, not
    # the modal count.
    median_count <- data.frame(
      displaced_id_value = rowDHSID,
      median_count = as.numeric(
        pb__weighted_quantile(
          x = singleDens.df$point_count,
          w = singleDens.df$layer,
          probs = 0.5
        )[[1]]
      )
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
    # Future-me: the heavy lift is embarrassingly parallel across displaced
    # points, so expose a couple of future.* knobs instead of baking in a
    # single plan. That makes it easier to tune on big runs without rewriting
    # the function again later.
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    
    cl <- parallel::makeCluster(n.cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    future::plan(future::cluster, workers = cl)
    
    numFeatures <- future.apply::future_lapply(
      X = row_index,
      FUN = one_comm,
      future.seed = future.seed,
      future.scheduling = future.scheduling
    )
  } else {
    numFeatures <- lapply(row_index, one_comm)
  }
  
  feature_probabilities <- iotools::fdrbind(lapply(numFeatures, `[[`, "feature_probabilities"))
  count_distribution <- iotools::fdrbind(lapply(numFeatures, `[[`, "count_distribution"))
  median_count <- iotools::fdrbind(lapply(numFeatures, `[[`, "median_count"))
  
  list(
    feature_probabilities = feature_probabilities,
    count_distribution = count_distribution,
    median_count = median_count
  )
}
