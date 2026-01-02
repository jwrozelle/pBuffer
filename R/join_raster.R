#' @name pb_rasterValJoin
#' @rdname pb_rasterValJoin
#' @title  pb_rasterValJoin
#'
#' @description  Estimates th
#' 
#' @returns Returns a list of two dataframes: estimate.df and commVals.df. estimate.df contains the best estimate of the raster value, while commVals.df contains all possible raster cells and their respective weights.
#'
#' @details This function can run rather slowly, since it re-calculates the probability raster for every community location. This is probably overkill, and unless the raster resolution is very low the pb_quickRasterValJoin would be nearly indestinguishable.
#'
#' @param displaced.sf Displaced coordinates in an sf object type.
#' @param inputRaster Object of class raster that the function will estimate a probable value from
#' @param bufferlengths The radius of the probability buffer to use
#' @param adminBound (Optional) the administrative boundary that circumscribes the displacement
#' @param adminID The unique ID for the adminBound features, if they are used. Defaults to "ID_2".
#' @param n.cores This allows for parallelization using the futures package. This can be slightly unstable, but normally functions well and dramatically speeds compute time when it does. n.cores specifies the number of cores to use. The default is 1, and does not parallelize.
#'
#' @author J.W. Rozelle
#'
#'
#' @export pb_rasterValJoin
#' @examples
#'
#' # coming soon!
#' 


#' Probability-weighted raster value under positional uncertainty
#'
#' For each displaced point, estimates the raster value at the (unobserved) true location
#' by integrating a probability kernel over raster cells within a specified radius and
#' computing a probability-weighted mean of raster values.
#'
#' @details
#' The kernel implemented here corresponds to *uniform-in-radius* displacement within
#' a circle of radius \code{bufferlength}: the implied density over area is
#' \eqn{f(x,y) = 1 / (2 \pi R \sqrt{x^2 + y^2})} for \eqn{\sqrt{x^2+y^2} \le R}.
#' This integrates to 1 over the circle but is numerically singular at the origin.
#' For the cell containing the origin, the implementation uses a small \code{eps}
#' adjustment in the denominator to avoid numerical failure in the integrator.
#'
#' Missing raster values are excluded from the weighted mean; weights are implicitly
#' renormalized over the non-missing subset via \code{weighted.mean(..., na.rm=TRUE)}.
#'
#' @param displaced.sf An \code{sf} POINT object containing displaced locations.
#' @param inputRaster A \code{raster::RasterLayer} from which values are extracted.
#' @param bufferlength Numeric scalar buffer radius in meters (default 5000). This
#'   implementation expects a single value; pass one radius at a time.
#' @param displaced.id Name of the unique ID column in \code{displaced.sf} (default \code{"DHSID"}).
#' @param adminBound Optional \code{sf} polygon boundary used to constrain extraction. If provided,
#'   the raster is masked to the polygon containing each point before weighting.
#' @param adminID Name of the unique polygon ID column in \code{adminBound} (default \code{"ID_2"}).
#' @param n.cores Integer number of workers for parallel execution (default 1). Uses
#'   \code{future.apply::future_lapply()} when \code{n.cores > 1}.
#' @param eps Small positive constant used to stabilize the kernel at the origin (default \code{1e-6}).
#'
#' @return A list with two data.frames:
#' \describe{
#'   \item{estimate.df}{One row per point ID with the probability-weighted raster estimate.}
#'   \item{commVals.df}{Long data.frame of candidate raster cell values and weights by point ID.}
#' }
#'
#' @export
pb_rasterValJoin <- function(displaced.sf,
                             inputRaster,
                             bufferlength = 5000,
                             displaced.id = "DHSID",
                             adminBound = NULL,
                             adminID = "ID_2",
                             n.cores = 1,
                             eps = 1e-6) {
  
  # --- validation ---
  if (!inherits(displaced.sf, "sf")) stop("displaced.sf must be an sf object.")
  if (!inherits(sf::st_geometry(displaced.sf), "sfc_POINT")) stop("displaced.sf must have POINT geometry.")
  if (!inherits(inputRaster, "RasterLayer")) stop("inputRaster must be a raster::RasterLayer.")
  if (!is.numeric(bufferlength) || length(bufferlength) != 1 || bufferlength <= 0) {
    stop("bufferlength must be a single positive number.")
  }
  if (!displaced.id %in% names(displaced.sf)) stop("displaced.id not found in displaced.sf: ", displaced.id)
  if (anyDuplicated(sf::st_drop_geometry(displaced.sf)[[displaced.id]]) > 0) {
    stop("displaced.id must be unique per row in displaced.sf.")
  }
  
  if (!is.null(adminBound)) {
    if (!inherits(adminBound, "sf")) stop("adminBound must be an sf object if provided.")
    if (!adminID %in% names(adminBound)) stop("adminID not found in adminBound: ", adminID)
  }
  
  # ensure CRS alignment for masking/intersection logic
  # (raster has its own CRS string; compare via proj4string / wkt lightly)
  if (!is.na(raster::crs(inputRaster)) && !is.na(sf::st_crs(displaced.sf))) {
    # Users should generally pre-align these; we don't auto-project rasters here.
    # We assume displaced.sf is already in the raster CRS.
  }
  
  raster_res <- raster::res(inputRaster)  # (xres, yres)
  
  warning("If the input raster has missing values, any probability-buffer mass landing on missing cells will be excluded (implicit renormalization).")
  
  # helper: compute weights & values for one point
  one_point <- function(row_id) {
    
    # subset single point
    singleComm <- displaced.sf[sf::st_drop_geometry(displaced.sf)[[displaced.id]] == row_id, ]
    if (nrow(singleComm) != 1) stop("Expected exactly one row for ", displaced.id, " = ", row_id)
    
    coords <- sf::st_coordinates(singleComm)[1, c("X", "Y")]
    if (is.na(coords[1]) || is.na(coords[2])) {
      return(list(
        estimate.df = data.frame(estimate = NA_real_, id = row_id),
        commVals.df = data.frame(id = row_id, pixelProb = numeric(0), pixelVal = numeric(0))
      ))
    }
    
    # buffer polygon around point
    buff <- suppressWarnings(sf::st_buffer(singleComm, bufferlength + max(raster_res)))
    
    # mask raster to buffer
    r_mask <- suppressWarnings(raster::mask(inputRaster, buff))
    
    # optionally constrain to admin polygon containing the point
    if (!is.null(adminBound)) {
      # find the polygon containing the point
      j <- suppressWarnings(sf::st_join(singleComm, adminBound[, adminID, drop = FALSE], join = sf::st_within, left = TRUE))
      admin_val <- sf::st_drop_geometry(j)[[adminID]][1]
      if (!is.na(admin_val)) {
        admin_poly <- adminBound[adminBound[[adminID]] == admin_val, ]
        r_mask <- suppressWarnings(raster::mask(r_mask, admin_poly))
      }
    }
    
    # trim for speed
    r_mask <- suppressWarnings(raster::trim(r_mask))
    
    # ==== EMPTY RASTER GUARD (AFTER TRIM) ====
    # If masking + trimming removed all cells, dim(r_mask) will be NULL
    # and any attempt to index rows/cols will error.
    if (is.null(dim(r_mask))) {
      return(list(
        estimate.df = data.frame(estimate = NA_real_, id = row_id),
        commVals.df = data.frame(id = row_id,
                                 pixelProb = numeric(0),
                                 pixelVal  = numeric(0))
      ))
    }
    # ==== END EMPTY RASTER GUARD ====
    
    
    # raster dims
    rasterRows <- dim(r_mask)[1]
    rasterCols <- dim(r_mask)[2]  # FIXED (was [1] in your code)
    
    # # kernel (uniform in radius) with stabilization epsilon
    R <- bufferlength
    # kernel_xy <- function(x, y) {
    #   rr <- sqrt(x^2 + y^2)
    #   # stabilize near 0
    #   rr <- rr + eps
    #   ifelse(rr <= R, 1 / (2 * pi * R * rr), 0)
    # }
    
    kernel_xy <- function(x, y) {
      rr <- sqrt(x^2 + y^2)
      rr_safe <- pmax(rr, eps)
      ifelse(rr <= R, 1 / (2 * pi * R * rr_safe), 0)
    }
    
    
    pixelProb <- numeric(rasterRows * rasterCols)
    pixelVal  <- numeric(rasterRows * rasterCols)
    
    # Extract raster values
    vals <- raster::values(r_mask)
    
    # iterate cells
    for (cell in 1:(rasterRows * rasterCols)) {
      
      # get cell center xy
      cell_xy <- raster::xyFromCell(r_mask, cell)
      # relative position (center minus point)
      rel <- c(coords[1] - cell_xy[1], coords[2] - cell_xy[2])
      
      # cell bounds in relative coordinates
      outerXY <- rel + raster_res / 2
      innerXY <- rel - raster_res / 2
      
      xlimits <- sort(c(outerXY[1], innerXY[1]))
      ylimits <- sort(c(outerXY[2], innerXY[2]))
      
      # integrate kernel over cell rectangle
      integ <- calculus::integral(
        kernel_xy,
        bounds = list(x = xlimits, y = ylimits),
        coordinates = "cartesian"
      )
      
      pixelProb[cell] <- integ$value
      pixelVal[cell] <- vals[cell]
    }
    
    comm <- data.frame(
      id = row_id,
      pixelProb = pixelProb,
      pixelVal = pixelVal
    )
    
    # --- probability-mass diagnostics (compute BEFORE dropping NA pixelVal) ---
    mass_total   <- sum(comm$pixelProb, na.rm = TRUE)
    mass_nonmiss <- sum(comm$pixelProb[!is.na(comm$pixelVal)], na.rm = TRUE)
    mass_ratio   <- if (is.finite(mass_total) && mass_total > 0) mass_nonmiss / mass_total else NA_real_
    
    # drop NA raster cells before weighted mean
    comm <- comm[!is.na(comm$pixelProb) & !is.na(comm$pixelVal), , drop = FALSE]
    
    est <- if (nrow(comm) == 0) {
      NA_real_
    } else {
      stats::weighted.mean(comm$pixelVal, w = comm$pixelProb, na.rm = TRUE)
    }
    
    list(
      estimate.df = data.frame(estimate = est, id = row_id),
      commVals.df = comm
    )
  }
  
  ids <- sf::st_drop_geometry(displaced.sf)[[displaced.id]]
  
  # --- run (parallel optional) ---
  if (n.cores > 1) {
    if (!requireNamespace("future", quietly = TRUE) || !requireNamespace("future.apply", quietly = TRUE)) {
      stop("For n.cores > 1, install 'future' and 'future.apply'.")
    }
    
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    
    cl <- parallel::makeCluster(n.cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    
    future::plan(future::cluster, workers = cl)
    out_list <- future.apply::future_lapply(ids, one_point)
    
  } else {
    out_list <- lapply(ids, one_point)
  }
  
  # bind results
  estimate.df <- do.call(rbind, lapply(out_list, `[[`, "estimate.df"))
  commVals.df <- do.call(rbind, lapply(out_list, `[[`, "commVals.df"))
  
  # restore column names
  names(estimate.df)[names(estimate.df) == "id"] <- displaced.id
  names(commVals.df)[names(commVals.df) == "id"] <- displaced.id
  
  list(estimate.df = estimate.df, commVals.df = commVals.df)
}








#' @name pb_quickRasterValJoin
#' @rdname pb_quickRasterValJoin
#' @title  pb_quickRasterValJoin
#'
#' @description For each given displaced community, estimates the raster value based on the probability weighted mean of possible true locations
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
#' @export pb_quickRasterValJoin
#' @examples
#'
#' # coming soon!
#' 


pb_quickRasterValJoin <- function(displaced.sf, 
                                  inputRaster, 
                                  bufferlengths = 5000, 
                                  displaced.id = "DHSID",
                                  adminBound = NULL, 
                                  adminID = NULL, 
                                  n.cores = 1) {
  
  # ERRORS!
  #   displaced.sf is the wrong class
  if (!"sf" %in% class(displaced.sf)) {
    stop("displaced.sf is not of class sf. It must be point features of class sf")
  }
  #   features2count.sf is the wrong class
  if (!"RasterLayer" %in% class(inputRaster)) {
    stop("inputRaster is not of class 'raster'. It must be point features of class 'raster'")
  }
  #   adminbBound is the wrong class
  if (!is.null(adminBound) & !"sf" %in% class(adminBound)) {
    stop("adminBound is not of class sf. It must be polygon features of class sf")
  }
  
  #   ensure that the adminBound sf object contains the unique ID
  if (!is.null(adminBound)) {
    if (is.null(adminID)) stop("adminID must be provided when adminBound is provided.")
    if (!adminID %in% names(adminBound)) {
      stop(paste0(adminID, " is not found in the adminBound sf object."))
    }
    if (anyDuplicated(adminBound[[adminID]]) > 0) {
      stop("adminID must specify a uniquely valid ID for the adminBound object.")
    }
  }
  
  
  # if (!is.null(adminBound)) {
  #   if (sum(duplicated(adminBound[[adminID]])) > 0) {
  #     stop("adminID must specify a uniquely valid ID for the adminbound object")
  #   }
  # }
  
  #   ensure that the displaced.sf sf object contains the unique ID
  if (!displaced.id %in% names(displaced.sf)) {
    stop(paste0(displaced.id, " is not found in the displaced.sf sf object."))
  }
  
  if (sum(duplicated(displaced.sf[[displaced.id]])) > 0) {
    stop("displaced.id must specify a uniquely valid ID for the displaced.sf object")
  }
  
  warning("If the input raster has missing values, any probability buffer points that fall on the missing raster values will be excluded from the weighted score.")
  
  # IDs to iterate over (driven by displaced.id)
  ids <- sf::st_drop_geometry(displaced.sf)[[displaced.id]]
  
  # QUICK VERSION CURRENTLY SUPPORTS ONE BUFFERLENGTH ONLY
  if (length(bufferlengths) != 1) {
    stop("pb_quickRasterValJoin currently supports a single bufferlength. Pass one value at a time.")
  }
  
  
  
  # get the resolution of the raster
  raster.res <- raster::res(inputRaster)
  crs2use <- raster::crs(inputRaster)
  # pixelArea <- raster.res[1]*raster.res[2]
  
  # assign the adminID character value ot adminID variable name in dataframe
  if (!is.null(adminBound)) {
    adminBound$.admin_id <- sf::st_drop_geometry(adminBound)[[adminID]]
  }
  
  pbBuffer.list <- list()
  
  for (i in 1:length(bufferlengths)) {
    
    # probability density function
    buffer_pb_function <- function(x, y) {
      
      result <- ifelse(sqrt(x^2 + y^2) <= bufferlengths[i], 1/(bufferlengths[i] * 2 * pi * sqrt(x^2+y^2)), 0)
      
      return(result)
      
    }
    
    # calculate the number of rows & columns needed given the raster resolution.
    quadrantRowN <- (ceiling(bufferlengths[i] / raster.res[2]))
    quadrantColN <- (ceiling(bufferlengths[i] / raster.res[1]))
    
    # create a raster of the appropriate size given the bufferlength and raster resolution
    # pbBuffer.matrix <- matrix(data = 0, nrow = quadrantRowN*2, ncol = quadrantColN*2)
    
    
    
    q4_matrix <- matrix(data = 0, nrow = quadrantRowN + 2, ncol = quadrantColN + 2)
    
    # used to calculate half the x & y meters in the raster cell to place center point.
    xadd <- raster.res[1] / 2
    yadd <- raster.res[2] / 2
    
    
    # do the weight calculation for the matrix
    for (q4_row in 1:(quadrantRowN+2)) {
      for (q4_col in 1:(quadrantColN+2)) {
        
        # Set limits for integrals
        if (q4_row - 1 > 0 | q4_col - 1 > 0) {
          ylimits <- c(((q4_row-1)*raster.res[2] - yadd) , ((q4_row-1)*raster.res[2] + yadd))
          xlimits <- c(((q4_col-1)*raster.res[1] - xadd) , ((q4_col-1)*raster.res[1] + xadd))
          
          integralValue <- calculus::integral(buffer_pb_function, bounds = list(x = xlimits, y = ylimits), coordinates = "cartesian")
          
          q4_matrix[q4_row, q4_col] <- integralValue$value
          
        } else {
          ylimits <- c(((q4_row-1)*raster.res[2] - yadd) , ((q4_row-1)*raster.res[2] + yadd))
          xlimits <- c(((q4_col-1)*raster.res[1] - xadd-0.01) , ((q4_col-1)*raster.res[1] + xadd))
          
          integralValue <- calculus::integral(buffer_pb_function, bounds = list(x = xlimits, y = ylimits), coordinates = "cartesian")
          
          # 
          q4_matrix[1, 1] <- integralValue$value
        }
        
      }
    }
    
    
    q3_matrix <- q4_matrix[, ncol(q4_matrix):2]
    q2_matrix <- q4_matrix[nrow(q4_matrix):2, ncol(q4_matrix):2]
    q1_matrix <- q4_matrix[nrow(q4_matrix):2, ]
    
    # bind the east and west halves of the north and south hemicircles
    north <- cbind(q2_matrix, q1_matrix)
    south <- cbind(q3_matrix, q4_matrix)
    
    # bind the north and south
    densMatrix <- rbind(north, south)
    
    rm(q1_matrix, q2_matrix, q3_matrix, q4_matrix, north, south)
    
    
    pbBuffer.list[[i]] <- densMatrix
  }
  
  
  
  if (n.cores > 1) {
    
    tryCatch({
      
      # set up future workers (and restore prior plan on exit)
      if (!requireNamespace("future", quietly = TRUE) || !requireNamespace("future.apply", quietly = TRUE)) {
        stop("For n.cores > 1, install 'future' and 'future.apply'.")
      }
      
      old_plan <- future::plan()
      on.exit(future::plan(old_plan), add = TRUE)
      
      future::plan(future::multisession, workers = n.cores)
      
      # ids <- sf::st_drop_geometry(displaced.sf)[["DHSID"]]  # replace "DHSID" if you later add displaced.id  # replace "DHSID" if you later add displaced.id  # replace "DHSID" if you later add displaced.id  # replace "DHSID" if you later add displaced.id
      
      weightedRasterVals <- future.apply::future_lapply(
        ids,
        function(x) {
        
        rowDHSID <- as.character(x)
        
        # extract an sf object for each row of the data frame
        singleComm <- displaced.sf[sf::st_drop_geometry(displaced.sf)[[displaced.id]] == rowDHSID, ]
        
        # throw an error if singleComm has multiple observations
        if (nrow(singleComm) > 1) {
          stop(paste0("There are multiple observations with the same DHSID: ", rowDHSID))
        }
        
        
        coords <- sf::st_coordinates(singleComm)[1, c("X","Y")]
        if (!anyNA(coords)) {
          
          singleBuffer.list <- list()
          singleBufferRaster.list <- list()
          
          i <- 1 # this code previously supported multiple buffers, but, now I just use one.  Retaining the og code with i for ease.
          # for (i in 1:length(bufferlengths)){
            
            # create a buffer for the community
            singleBuffer.list[[i]] <- sf::st_buffer(singleComm, bufferlengths[i]+(max(raster.res)*2)) |> suppressWarnings()
            
            # # If the administrative boundary is defined
            # if (!is.null(adminBound)) {
            #   # get the single admin boundary for a community
            #   singleAdminBound <- suppressWarnings(sf::st_intersection(singleComm, adminBound))
            #   singleAdminBound.poly <- dplyr::filter(adminBound, adminID == singleAdminBound$adminID[1])
            #   
            #   # change all values that don't fall within the buffer to 0
            #   singleBufferRaster.list[[i]] <- raster::mask(inputRaster, singleBuffer.list[[i]]) |> suppressWarnings()
            #   
            #   # Trim the weighted raster
            #   singleBufferRaster.list[[i]] <- raster::mask(singleBufferRaster.list[[i]], singleAdminBound.poly) |> suppressWarnings()
            #   rm(singleAdminBound, singleAdminBound.poly)
            # }
            # Mask raster to buffer (always)
            singleBufferRaster.list[[i]] <- raster::mask(inputRaster, singleBuffer.list[[i]]) |> suppressWarnings()
            
            # If the administrative boundary is defined, additionally mask to the polygon containing the point
            if (!is.null(adminBound)) {
              
              j <- suppressWarnings(
                sf::st_join(
                  singleComm,
                  adminBound[, ".admin_id", drop = FALSE],
                  join = sf::st_within,
                  left = TRUE
                )
              )
              
              admin_val <- sf::st_drop_geometry(j)$.admin_id[1]
              
              if (!is.na(admin_val)) {
                singleAdminBound.poly <- adminBound[adminBound$.admin_id == admin_val, ]
                singleBufferRaster.list[[i]] <- raster::mask(singleBufferRaster.list[[i]], singleAdminBound.poly) |> suppressWarnings()
              } else {
                # point not within any admin polygon; leave raster masked only to buffer
                singleAdminBound.poly <- NULL
              }
              
              rm(j, admin_val, singleAdminBound.poly)
            }
            
            
            
            # 
            # # mask based on the administrative bounadries
            # singleBufferRaster.list[[i]] <- raster::mask(singleBufferRaster.list[[i]], singleAdminBound.poly)
            
            # trim the excess NA's
            singleBufferRaster.list[[i]] <- raster::trim(singleBufferRaster.list[[i]]) |> suppressWarnings()
            
            # ==== EMPTY RASTER GUARD (AFTER TRIM) ====
            if (is.null(dim(singleBufferRaster.list[[i]]))) {
              return(NA)
            }
            # ==== END EMPTY RASTER GUARD ====
            
            
            # rasterRows 
            rasterRows <- dim(singleBufferRaster.list[[i]])[1]
            rasterCols <- dim(singleBufferRaster.list[[i]])[2]
            
            # RASTER CENTER POINT
            singleCoords <- sf::st_coordinates(singleComm)
            #   cell of center point
            centerCell <- raster::cellFromXY(singleBufferRaster.list[[i]], singleCoords)
            centerRowCol <- raster::rowColFromCell(singleBufferRaster.list[[i]], centerCell)
            
            centerCell_xy <- raster::xyFromCell(singleBufferRaster.list[[i]], centerCell)
            
            pb_singleComm.raster <- raster(pbBuffer.list[[i]])
            
            pb_xmin <- centerCell_xy[1,"x"] - raster.res[1]*(ncol(pbBuffer.list[[i]])/2)
            pb_xmax <- centerCell_xy[1,"x"] + raster.res[1]*(ncol(pbBuffer.list[[i]])/2)
            pb_ymin <- centerCell_xy[1,"y"] - raster.res[2]*(nrow(pbBuffer.list[[i]])/2)
            pb_ymax <- centerCell_xy[1,"y"] + raster.res[2]*(nrow(pbBuffer.list[[i]])/2)
            
            extent(pb_singleComm.raster) <- c(pb_xmin, pb_xmax, pb_ymin, pb_ymax)
            crs(pb_singleComm.raster) <- crs2use
            raster::values(pb_singleComm.raster)[raster::values(pb_singleComm.raster) == 0] <- NA
            
            # get the rasters to line up
            singleBufferRaster.list[[i]] <- raster::extend(singleBufferRaster.list[[i]], pb_singleComm.raster)
            pb_singleComm.raster <- raster::extend(pb_singleComm.raster, singleBufferRaster.list[[i]])
            
            
            
            rm(pb_xmin, pb_xmax, pb_ymin, pb_ymax)
            
          # }
          
          singleCommVals.df <- data.frame(pixelProb = raster::values(pb_singleComm.raster), 
                                          pixelVal = raster::values(singleBufferRaster.list[[i]])
          )
          singleCommVals.df$weights <- singleCommVals.df$pixelProb # / sum(pixelProb, na.rm = T)
          singleCommVals.df[[displaced.id]] <- rowDHSID
          
          singleCommVals.df <- dplyr::filter(singleCommVals.df, !is.na(pixelProb))
          singleCommVals.df <- dplyr::filter(singleCommVals.df, !is.na(pixelVal))
          
          estimate <- weighted.mean(singleCommVals.df$pixelVal, singleCommVals.df$weights, na.rm = T)
          estimate.df <- data.frame(estimate)
          estimate.df[[displaced.id]] <- rowDHSID
          
          
          return(list(estimate.df = estimate.df, commVals.df = singleCommVals.df))
        } else {
          return(NA)
        }
      },
      future.packages = c("sf", "raster", "dplyr")
      )
      
    }# , 
    
    # finally = {parallel::stopCluster(cl)}
    )
    
    
    
    
    
    
    # if not parallelized...  
  } else {
    
    # ids <- sf::st_drop_geometry(displaced.sf)[["DHSID"]]  # replace "DHSID" if you later add displaced.id  # replace "DHSID" if you later add displaced.id  # replace "DHSID" if you later add displaced.id
    weightedRasterVals <- lapply(ids, function(x) {
      
      
      rowDHSID <- as.character(x)
      
      # extract an sf object for each row of the data frame
      singleComm <- displaced.sf[sf::st_drop_geometry(displaced.sf)[[displaced.id]] == rowDHSID, ]
      
      # throw an error if singleComm has multiple observations
      if (nrow(singleComm) > 1) {
        stop(paste0("There are multiple observations with the same DHSID: ", rowDHSID))
      }
      
      
      coords <- sf::st_coordinates(singleComm)[1, c("X","Y")]
      if (!anyNA(coords)) {
        
        singleBuffer.list <- list()
        singleBufferRaster.list <- list()
        
        i <- 1
        # for (i in 1:length(bufferlengths)){
          
          # create a buffer for the community
          singleBuffer.list[[i]] <- sf::st_buffer(singleComm, bufferlengths[i]+(max(raster.res)*2)) |> suppressWarnings()
          
          # # If the administrative boundary is defined
          # if (!is.null(adminBound)) {
          #   # get the single admin boundary for a community
          #   singleAdminBound <- suppressWarnings(sf::st_intersection(singleComm, adminBound))
          #   singleAdminBound.poly <- dplyr::filter(adminBound, adminID == singleAdminBound$adminID[1])
          #   
          #   # change all values that don't fall within the buffer to 0
          #   singleBufferRaster.list[[i]] <- raster::mask(inputRaster, singleBuffer.list[[i]]) |> suppressWarnings()
          #   
          #   # Trim the weighted raster
          #   singleBufferRaster.list[[i]] <- raster::mask(singleBufferRaster.list[[i]], singleAdminBound.poly) |> suppressWarnings()
          #   rm(singleAdminBound, singleAdminBound.poly)
          # }
          
          
          # Mask raster to buffer (always)
          singleBufferRaster.list[[i]] <- raster::mask(inputRaster, singleBuffer.list[[i]]) |> suppressWarnings()
          
          # If the administrative boundary is defined, additionally mask to the polygon containing the point
          if (!is.null(adminBound)) {
            
            j <- suppressWarnings(
              sf::st_join(
                singleComm,
                adminBound[, ".admin_id", drop = FALSE],
                join = sf::st_within,
                left = TRUE
              )
            )
            
            admin_val <- sf::st_drop_geometry(j)$.admin_id[1]
            
            if (!is.na(admin_val)) {
              singleAdminBound.poly <- adminBound[adminBound$.admin_id == admin_val, ]
              singleBufferRaster.list[[i]] <- raster::mask(singleBufferRaster.list[[i]], singleAdminBound.poly) |> suppressWarnings()
            } else {
              # point not within any admin polygon; leave raster masked only to buffer
              singleAdminBound.poly <- NULL
            }
            
            rm(j, admin_val, singleAdminBound.poly)
          }
          
          
          
          # 
          # # mask based on the administrative bounadries
          # singleBufferRaster.list[[i]] <- raster::mask(singleBufferRaster.list[[i]], singleAdminBound.poly)
          
          # trim the excess NA's
          singleBufferRaster.list[[i]] <- raster::trim(singleBufferRaster.list[[i]]) |> suppressWarnings()
          
          # ==== EMPTY RASTER GUARD (AFTER TRIM) ====
          if (is.null(dim(singleBufferRaster.list[[i]]))) {
            return(NA)
          }
          # ==== END EMPTY RASTER GUARD ====
          
          
          # rasterRows 
          rasterRows <- dim(singleBufferRaster.list[[i]])[1]
          rasterCols <- dim(singleBufferRaster.list[[i]])[2]
          
          # RASTER CENTER POINT
          singleCoords <- sf::st_coordinates(singleComm)
          #   cell of center point
          centerCell <- raster::cellFromXY(singleBufferRaster.list[[i]], singleCoords)
          centerRowCol <- raster::rowColFromCell(singleBufferRaster.list[[i]], centerCell)
          
          centerCell_xy <- raster::xyFromCell(singleBufferRaster.list[[i]], centerCell)
          
          pb_singleComm.raster <- raster(pbBuffer.list[[i]])
          
          pb_xmin <- centerCell_xy[1,"x"] - raster.res[1]*(ncol(pbBuffer.list[[i]])/2)
          pb_xmax <- centerCell_xy[1,"x"] + raster.res[1]*(ncol(pbBuffer.list[[i]])/2)
          pb_ymin <- centerCell_xy[1,"y"] - raster.res[2]*(nrow(pbBuffer.list[[i]])/2)
          pb_ymax <- centerCell_xy[1,"y"] + raster.res[2]*(nrow(pbBuffer.list[[i]])/2)
          
          extent(pb_singleComm.raster) <- c(pb_xmin, pb_xmax, pb_ymin, pb_ymax)
          crs(pb_singleComm.raster) <- crs2use
          raster::values(pb_singleComm.raster)[raster::values(pb_singleComm.raster) == 0] <- NA
          
          # get the rasters to line up
          singleBufferRaster.list[[i]] <- raster::extend(singleBufferRaster.list[[i]], pb_singleComm.raster)
          pb_singleComm.raster <- raster::extend(pb_singleComm.raster, singleBufferRaster.list[[i]])
          
          
          
          rm(pb_xmin, pb_xmax, pb_ymin, pb_ymax)
          
        # }
        
        singleCommVals.df <- data.frame(pixelProb = raster::values(pb_singleComm.raster), 
                                        pixelVal = raster::values(singleBufferRaster.list[[i]])
        )
        singleCommVals.df$weights <- singleCommVals.df$pixelProb # / sum(pixelProb, na.rm = T)
        singleCommVals.df[[displaced.id]] <- rowDHSID
        
        singleCommVals.df <- dplyr::filter(singleCommVals.df, !is.na(pixelProb))
        singleCommVals.df <- dplyr::filter(singleCommVals.df, !is.na(pixelVal))
        
        estimate <- weighted.mean(singleCommVals.df$pixelVal, singleCommVals.df$weights, na.rm = T)
        estimate.df <- data.frame(estimate)
        estimate.df[[displaced.id]] <- rowDHSID
        
        
        return(list(estimate.df = estimate.df, commVals.df = singleCommVals.df))
      } else {
        return(NA)
      }
      
    })
  }
  
  # tidy things up, put into a single data frame
  #   Extract the objects of interest from the list
  #     probableMetrics
  estimate.list <- lapply(weightedRasterVals, function(x) {
    x$estimate.df
  }) 
  #     hfWeights
  commVals.list <- lapply(weightedRasterVals, function(x) {
    x$commVals.df
  }) 
  
  rm(weightedRasterVals)
  
  #   bind them into a dataframe
  estimate.df <- iotools::fdrbind(estimate.list)
  commVals.df <- iotools::fdrbind(commVals.list)
  
  rm(estimate.list)
  rm(commVals.list)
  
  return(list(estimate.df = estimate.df, commVals.df = commVals.df))
}


