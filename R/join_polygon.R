#' @name ml_polygon
#' @rdname ml_polygon
#' @title  ml_polygon
#'
#' @description  Determines the maximum likely polygon for a given point
#' 
#' @returns returns a dataframe with one row containing the data from most likely health facility joined, and its estimated probability.
#' 
#' @details Replaces the weightedClose function
#'
#'
#' @param bufferPoints.sf Buffer points of a given displaced location containing weights
#' @param polygons.sf Polygons to join the displaced points to. For now, the function assumes there is a column with unique ID for each polygon called "SPAID"
#' @param polygons.id Name of column of polygons.id sf object that contains a unique ID, default is "SPAID".
#' @param weightsCol Name of the column in the bufferPoints.sf object that contains the weights to use, defaults to "layer".
#'
#' @author J.W. Rozelle
#'
#'
#' @export ml_polygon
#' @examples
#'
#' # coming soon!
#' 
#'



ml_polygon <- function(bufferPoints.sf, polygons.sf, polygons.id = "SPAID", weightsCol = "layer") {
  
  #----validation-------
  pb_check_sf(bufferPoints.sf, "bufferPoints.sf")
  pb_check_sf(polygons.sf, "polygons.sf")
  
  pb_check_cols(polygons.sf, polygons.id, "polygons.sf")
  pb_check_cols(bufferPoints.sf, weightsCol, "bufferPoints.sf")
  
  
  # Standardize column names BEFORE intersection so they propagate into the intersection result
  polygons.sf$SPAID <- polygons.sf[[polygons.id]]
  bufferPoints.sf$layer <- bufferPoints.sf[[weightsCol]]
  
  intersection <- st_intersection(
    polygons.sf[, c("SPAID")],
    bufferPoints.sf[, c("layer")]
  )
  
  # normalize intersection weights (not bufferPoints.sf)
  intersection$layer <- pb_normalize_weights(intersection$layer, "intersection weights")
  
  # non_missing <- length(intersection$SPAID)
  
  # group all intersections by the number of possible polygons
  int_result <- intersection %>% 
    st_drop_geometry() %>% 
    group_by(SPAID) %>% 
    summarize(weight = sum(layer, na.rm = TRUE))
  
  # merge the grouped possible polygons and probabilities with the data in the polygons. Keep only the possible polygons
  includedScores <- merge(
    int_result,
    st_drop_geometry(polygons.sf[c("SPAID")]),
    all.x = TRUE, all.y = FALSE,
    by = "SPAID"
  )
  
  # Filter to the maximum likely polygon
  mostLikelyHF.df <- filter(int_result, weight == max(weight, na.rm = T)) |> as.data.frame()
  
  # throw a warning if linked to more than one health facility
  if (nrow(mostLikelyHF.df) > 1) {
    warning(paste0(
      "Warning: linked to more than one polygon. The first of ",
      nrow(mostLikelyHF.df),
      " tied polygons is used."
    ))
    mostLikelyHF.df <- mostLikelyHF.df[1, , drop = FALSE]
  }
  
  
  
  return(mostLikelyHF.df)
}


#' @name pb_polyJoin
#' @rdname pb_polyJoin
#' @title  pb_polyJoin
#'
#' @description  Estimates th
#' 
#' @returns Returns a dataframe of length nrow(displaced.sf), with a DHSID column containing the unique ID of displaced.sf and SPAID containing the unique ID of the most likely polygon that a given point falls in.
#'
#'
#' @param displaced.sf Displaced coordinates in an sf object type.
#' @param polygons.sf Polygons to join the displaced points to. For now, the function assumes there is a column with unique ID for each polygon called "SPAID"
#' @param displaced.id The name of the column that contains a unique ID for the displaced communities.
#' @param densityBuffer A density buffer list object created with the pb_integratedDens or pb_Density functions.
#' @param adminBound (Optional) the administrative boundary that circumscribes the displacement
#' @param adminID The unique ID for the adminBound featurs, defaults to "ID_2".
#' @param n.cores This allows for parallelization using the futures package. This can be slightly unstable, but normally functions well and dramatically speeds compute time when it does. n.cores specifies the number of cores to use. The default is 1, and does not parallelize.
#'
#' @author J.W. Rozelle
#'
#'
#' @export pb_Density
#' @examples
#'
#' # coming soon!
#' 
#'





pb_polyJoin <- function(displaced.sf, polygons.sf, displaced.id = "DHSID", densityBuffer = NULL, adminBound = haiti_adm2, adminID = "ID_2", n.cores = 1) {
  
  # Check to make sure that objects which should be in sf format are in sf format
  #   displaced.sf
  if (!"sf" %in% class(displaced.sf)) {
    stop("The specified displaced.sf is not in sf format. It must be in sf format for this function to work.")
  }
  #   polygons.sf
  if (!"sf" %in% class(polygons.sf)) {
    stop("The specified polygons.sf is not in sf format. It must be in sf format for this function to work.")
  }
  #   adminBound
  if (!"sf" %in% class(displaced.sf)) {
    stop("The specified adminBound is not in sf format. It must be in sf format for this function to work.")
  }
  
  
  # ensure that the adminBound sf object contains the unique ID
  if (!adminID %in% names(adminBound)) {
    stop(paste0(adminID, " is not found in the adminBound sf object."))
  }
  
  if (!is.null(adminBound)) {
    if (sum(duplicated(adminBound[[adminID]])) > 0) {
      stop("adminID must specify a uniquely valid ID for the adminbound object")
    }
  }
  
  # ensure that the displaced.sf sf object contains the unique ID
  if (!displaced.id %in% names(displaced.sf)) {
    stop(paste0(displaced.id, " is not found in the displaced.sf sf object."))
  }
  
  if (sum(duplicated(displaced.sf[[displaced.id]])) > 0) {
    stop("displaced.id must specify a uniquely valid ID for the displaced.sf object")
  }
  
  # rename the adminboundary unique ID to adminID
  adminBound$adminID <- adminBound[[adminID]]
  
  
  if (n.cores > 1) {
    
    tryCatch({
      
      # set up the workers
      cl <- parallel::makeCluster(n.cores)
      plan(cluster, workers = cl)
      
      closestHFs <- future_apply(st_drop_geometry(displaced.sf), 1, function(x) {
        
        rowDHSID <- x[displaced.id]
        
        # extract an sf object for each row of the data frame
        singleComm <- displaced.sf[displaced.sf[[displaced.id]] == rowDHSID, ]
        crs2use <- raster::crs(singleComm)
        
        # If the administrative boundary is defined
        if (!is.null(adminBound)) {
          # get the single admin boundary for a community
          singleAdminBound <- suppressWarnings(sf::st_intersection(singleComm, adminBound))
          singleAdminBound.poly <- dplyr::filter(adminBound, adminID == singleAdminBound$adminID[1])
        }
        if (st_coordinates(singleComm)[2] != 0) {
          
          # rasterize the displacement buffer around the single community
          singleDens.raster <- rasterizeDisplacement(
            densityBuffer, 
            initialCoords = st_coordinates(singleComm), 
            inputCRS = crs2use
          )
          
          
          # If the administrative boundary is specified, trim the probability buffer by the polygon layer
          if (!is.null(adminBound)) {
            singleDens.raster <- pb_trim_probRaster_to_point_admin(
              probRaster = singleDens.raster,
              point_sf   = singleComm,
              adminBound = adminBound,
              adminID    = "adminID"
            )
          }
          
          
          
          # turn these into a point object
          singleDens.sf <- st_rasterAsPoint(singleDens.raster)
          rm(singleDens.raster)
          
          # Get the most probable nearby health facility
          probableHF <- ml_polygon(singleDens.sf, polygons.sf)
          
          # put dhsid, spaid, and weight into a merged, single-row data frame
          results.df <- data.frame(DHSID = rowDHSID)
          results.df <- cbind(results.df, probableHF)
        } else {
          results.df <- NA
        }
        
        return(results.df)
        
      })
    },
    
    # always stop the cluster, even if future_apply fails
    finally = {parallel::stopCluster(cl)}
    )
    
  } else {
    
    closestHFs <- apply(st_drop_geometry(displaced.sf), 1, function(x) {
      
      rowDHSID <- x[displaced.id]
      
      # extract an sf object for each row of the data frame
      singleComm <- displaced.sf[displaced.sf[[displaced.id]] == rowDHSID, ]
      crs2use <- raster::crs(singleComm)
      
      # If the administrative boundary is defined
      if (!is.null(adminBound)) {
        # get the single admin boundary for a community
        singleAdminBound <- suppressWarnings(sf::st_intersection(singleComm, adminBound))
        singleAdminBound.poly <- dplyr::filter(adminBound, adminID == singleAdminBound$adminID[1])
      }
      if (st_coordinates(singleComm)[2] != 0) {
        
        # rasterize the displacement buffer around the single community
        singleDens.raster <- rasterizeDisplacement(
          densityBuffer, 
          initialCoords = st_coordinates(singleComm), 
          inputCRS = crs2use
        )
        
        
        # If the administrative boundary is specified, trim the probability buffer by the polygon layer
        if (!is.null(adminBound)) {
          singleDens.raster <- pb_trim_probRaster_to_point_admin(
            probRaster = singleDens.raster,
            point_sf   = singleComm,
            adminBound = adminBound,
            adminID    = "adminID"
          )
        }
        
        
        
        # turn these into a point object
        singleDens.sf <- st_rasterAsPoint(singleDens.raster)
        rm(singleDens.raster)
        
        # Get the most probable nearby health facility
        probableHF <- ml_polygon(singleDens.sf, polygons.sf)
        
        # put dhsid, spaid, and weight into a merged, single-row data frame
        results.df <- data.frame(DHSID = rowDHSID)
        results.df <- cbind(results.df, probableHF)
      } else {
        results.df <- NA
      }
      
      return(results.df)
      
    })
  }
  
  closestHFs.df <- iotools::fdrbind(closestHFs) |> as.data.frame()
  
  return(closestHFs.df)
}


#' @name pb_weightedPolyJoin
#' @rdname pb_weightedPolyJoin
#' @title  pb_weightedPolyJoin
#'
#' @description For each given displaced community, estimates the raster value based on the probability weighted mean of possible true locations
#' 
#' @returns Returns a dataframe of length nrow(displaced.sf), with a DHSID column containing the unique ID of displaced.sf and SPAID containing the unique ID of the most likely polygon that a given point falls in.
#'
#'
#' @param displaced.sf Displaced coordinates in an sf object type.
#' @param polygons.sf Polygons to join the displaced points to. For now, the function assumes there is a column with unique ID for each polygon called "SPAID"
#' @param metrics Character vector of column names containing the numeric metrics you wish to calculate.
#' @param densityBuffer A density buffer list object created with the pb_integratedDens function.
#' @param adminBound (Optional) the administrative boundary that circumscribes the displacement.
#' @param adminID The unique ID for the adminBound features, defaults to "ID_2".
#' @param n.cores This allows for parallelization using the futures package. This can be slightly unstable, but normally functions well and dramatically speeds compute time when it does. n.cores specifies the number of cores to use. The default is 1, and does not parallelize.
#'
#' @author J.W. Rozelle
#'
#'
#' @export pb_weightedPolyJoin
#' @examples
#'
#' # coming soon!
#' 

pb_weightedPolyJoin <- function(displaced.sf, 
                                polygons.sf, 
                                displaced.id = "DHSID", 
                                metrics = c("sri_score", 
                                            "sri_basicamenities", 
                                            "sri_basicequip", 
                                            "sri_diagcapacity", 
                                            "sri_infprev", 
                                            "sri_med"
                                ), 
                                densityBuffer = NULL, 
                                adminBound = NULL, 
                                adminID = "ID_2",
                                n.cores = 1) {
  
  displaced.sf$DHSID <- displaced.sf[[displaced.id]]
  
  
  if (n.cores > 1) {
    
    tryCatch({
      
      # set up the workers
      cl <- parallel::makeCluster(n.cores)
      plan(cluster, workers = cl)
      
      closestHFs <- future_apply(st_drop_geometry(displaced.sf), 1, function(x) {
        
        rowDHSID <- x[displaced.id]
        
        # extract an sf object for each row of the data frame
        singleComm <- displaced.sf[displaced.sf[[displaced.id]] == rowDHSID, ]
        crs2use <- raster::crs(singleComm)
        
        if (st_coordinates(singleComm)[2] != 0) {
          
          # rasterize the displacement buffer around the single community
          singleDens.raster <- rasterizeDisplacement(
            densityBuffer, 
            initialCoords = st_coordinates(singleComm), 
            inputCRS = crs2use
          )
          
          # If the administrative boundary is defined
          if (!is.null(adminBound)) {
            # get the single admin boundary for a community
            singleAdminBound <- suppressWarnings(sf::st_intersection(singleComm, adminBound))
            singleAdminBound.poly <- dplyr::filter(adminBound, adminID == singleAdminBound$adminID[1])
            
            # Trim the weighted raster
            singleDens.raster <- trimProbBuff(singleDens.raster, adminBound = adminBound)
            rm(singleAdminBound, singleAdminBound.poly)
          }
          
          
          
          # turn these into a point object
          singleDens.sf <- st_rasterAsPoint(singleDens.raster)
          rm(singleDens.raster)
          
          # Get the most probable nearby health facility
          probableMetrics.result <- pb_weighted_metrics(
            bufferPoints.sf = singleDens.sf,
            polygons.sf = polygons.sf,
            poly_id = "SPAID",
            weightVar = "layer",
            metrics = metrics
          )
          probableMetrics.df <- probableMetrics.result$weightedMetrics.df
          probableMetrics.df$DHSID <- rowDHSID
          probableMetrics.df$n_linked <- nrow(probableMetrics.result$hfWeights)
          
          hfWeights.df <- probableMetrics.result$hfWeights
          hfWeights.df$DHSID <- rowDHSID
          
          
          
        } else {
          probableMetrics.df <- NA
          hfWeights.df <- NA
        }
        
        return(list(probableMetrics.df = probableMetrics.df, hfWeights.df = hfWeights.df))
        
      })
      
    }, 
    finally = {parallel::stopCluster(cl)}
    )
    
    
  } else {
    
    closestHFs <- apply(st_drop_geometry(displaced.sf), 1, function(x) {
      
      rowDHSID <- x[displaced.id]
      
      # extract an sf object for each row of the data frame
      singleComm <- displaced.sf[displaced.sf[[displaced.id]] == rowDHSID, ]
      crs2use <- raster::crs(singleComm)
      
      if (st_coordinates(singleComm)[2] != 0) {
        
        # rasterize the displacement buffer around the single community
        singleDens.raster <- rasterizeDisplacement(
          densityBuffer, 
          initialCoords = st_coordinates(singleComm), 
          inputCRS = crs2use
        )
        
        # If the administrative boundary is defined
        if (!is.null(adminBound)) {
          # get the single admin boundary for a community
          singleAdminBound <- suppressWarnings(sf::st_intersection(singleComm, adminBound))
          singleAdminBound.poly <- dplyr::filter(adminBound, adminID == singleAdminBound$adminID[1])
          
          # Trim the weighted raster
          singleDens.raster <- trimProbBuff(singleDens.raster, adminBound = adminBound)
          rm(singleAdminBound, singleAdminBound.poly)
        }
        
        
        
        # turn these into a point object
        singleDens.sf <- st_rasterAsPoint(singleDens.raster)
        rm(singleDens.raster)
        
        # Get the most probable nearby health facility
        probableMetrics.result <- pb_weighted_metrics(
          bufferPoints.sf = singleDens.sf,
          polygons.sf = polygons.sf,
          poly_id = "SPAID",
          weightVar = "layer",
          metrics = metrics
        )
        probableMetrics.df <- probableMetrics.result$weightedMetrics.df
        probableMetrics.df$DHSID <- rowDHSID
        probableMetrics.df$n_linked <- nrow(probableMetrics.result$hfWeights)
        
        hfWeights.df <- probableMetrics.result$hfWeights
        hfWeights.df$DHSID <- rowDHSID
        
        
        
      } else {
        probableMetrics.df <- NA
        hfWeights.df <- NA
      }
      
      return(list(probableMetrics.df = probableMetrics.df, hfWeights.df = hfWeights.df))
      
    })
  }
  
  
  # tidy things up, put into a single data frame
  #   Extract the objects of interest from the list
  #     probableMetrics
  probableMetrics.list <- lapply(closestHFs, function(x) {
    x$probableMetrics.df
  }) 
  #     hfWeights
  hfWeights.list <- lapply(closestHFs, function(x) {
    x$hfWeights.df
  }) 
  
  rm(closestHFs)
  
  #   bind them into a dataframe
  probableMetrics.df <- iotools::fdrbind(probableMetrics.list)
  hfWeights.df <- iotools::fdrbind(hfWeights.list)
  
  rm(probableMetrics.list)
  rm(hfWeights.list)
  
  
  
  return(list(probableMetrics.df = probableMetrics.df, hfWeights.df = hfWeights.df))
}


