#' @name pb_radiusLimitJoin
#' @rdname pb_radiusLimitJoin
#' @title  pb_radiusLimitJoin
#'
#' @description For each given displaced community, estimates the raster value based on the probability weighted mean of possible true locations
#' 
#' @returns a vector of the most likely count of point features within a defined radius of the displaced coordinates.
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
#' @description For each given displaced community, estimates the raster value based on the probability weighted mean of possible true locations
#' 
#' @returns a vector of the most likely count of point features within a defined radius of the displaced coordinates.
#'
#'
#' @param displaced.sf Displaced coordinates in an sf object type.
#' @param features2count.sf sf point feature that is counted.
#' @param radiusLength a numeric vector to determine the distance from a displaced community within which the number of features2count.sf will be counted.
#' @param densityBuffer A density buffer list object created with the pb_integratedDens function.
#' @param adminBound (Optional) the administrative boundary that circumscribes the displacement.
#' @param adminID The unique ID for the adminBound features, defaults to "ID_2".
#' @param n.cores This allows for parallelization using the futures package. This can be slightly unstable, but normally functions well and dramatically speeds compute time when it does. n.cores specifies the number of cores to use. The default is 1, and does not parallelize.
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
                             densityBuffer = NULL, 
                             adminBound = NULL, 
                             adminID = NULL,
                             n.cores = 1) {
  
  pb_check_sf(displaced.sf, "displaced.sf")
  pb_check_geom_type(displaced.sf, "POINT", "displaced.sf")
  pb_check_crs(displaced.sf, "displaced.sf")
  pb_check_projected(displaced.sf, "displaced.sf")
  
  pb_check_sf(features2count.sf, "features2count.sf")
  pb_check_crs(features2count.sf, "features2count.sf")
  
  pb_check_cols(displaced.sf, displaced.id, "displaced.sf")
  pb_check_unique_id(sf::st_drop_geometry(displaced.sf), displaced.id, "displaced.sf")
  
  pb_check_pos_scalar(radiusLength, "radiusLength")
  
  if (is.null(densityBuffer)) {
    stop("densityBuffer must be provided (e.g., pb_integratedDensity(...) or pb_mcDensity(...)).", call. = FALSE)
  }
  
  if (!is.null(adminBound)) {
    pb_check_sf(adminBound, "adminBound")
    if (is.null(adminID)) stop("adminID must be provided when adminBound is provided.", call. = FALSE)
    pb_check_cols(adminBound, adminID, "adminBound")
    pb_check_unique_id(sf::st_drop_geometry(adminBound), adminID, "adminBound")
  }
  
  displaced.sf$DHSID <- displaced.sf[[displaced.id]]
  
  if (!is.null(adminBound)) {
    adminBound$adminID <- sf::st_drop_geometry(adminBound)[[adminID]]
  }
  
  one_comm <- function(xrow) {
    
    rowDHSID <- xrow[[displaced.id]]
    singleComm <- displaced.sf[displaced.sf[[displaced.id]] == rowDHSID, ]
    
    xy <- sf::st_coordinates(singleComm)
    if (nrow(xy) < 1 || is.na(xy[1,2]) || xy[1,2] == 0) {
      return(list(ml_featureCount.df = NA, possibleNumFeatures.df = NA))
    }
    
    crs2use <- raster::crs(singleComm)
    
    singleDens.raster <- rasterizeDisplacement(
      densityBuffer,
      initialCoords = xy[1, 1:2],
      inputCRS = crs2use
    )
    
    if (!is.null(adminBound)) {
      singleDens.raster <- pb_trim_probRaster_to_point_admin(
        probRaster = singleDens.raster,
        point_sf   = singleComm,
        adminBound = adminBound,
        adminID    = "adminID"
      )
    }
    
    singleDens.sf <- st_rasterAsPoint(singleDens.raster)
    rm(singleDens.raster)
    
    singleBuffs.sf <- sf::st_buffer(singleDens.sf, radiusLength)
    singleBuffs.sf$pt_count <- lengths(sf::st_intersects(singleBuffs.sf, features2count.sf))
    
    numFeatures.df <- sf::st_drop_geometry(singleBuffs.sf) |>
      dplyr::group_by(pt_count) |>
      dplyr::summarise(weights = sum(layer, na.rm = TRUE), .groups = "drop") |>
      as.data.frame()
    
    ml_featureCount.df <- dplyr::filter(numFeatures.df, weights == max(weights, na.rm = TRUE)) |>
      as.data.frame()
    ml_featureCount.df[[displaced.id]] <- rowDHSID
    
    numFeatures.df[[displaced.id]] <- rowDHSID
    
    list(ml_featureCount.df = ml_featureCount.df, possibleNumFeatures.df = numFeatures.df)
  }
  
  if (n.cores > 1) {
    cl <- parallel::makeCluster(n.cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    future::plan(future::cluster, workers = cl)
    
    numFeatures <- future.apply::future_apply(sf::st_drop_geometry(displaced.sf), 1, one_comm)
  } else {
    numFeatures <- apply(sf::st_drop_geometry(displaced.sf), 1, one_comm)
  }
  
  ml_featureCount.df <- iotools::fdrbind(lapply(numFeatures, `[[`, "ml_featureCount.df"))
  possibleNumFeatures.df <- iotools::fdrbind(lapply(numFeatures, `[[`, "possibleNumFeatures.df"))
  
  list(ml_featureCount.df = ml_featureCount.df, possibleNumFeatures.df = possibleNumFeatures.df)
}