#' @name pb_Density
#' @rdname pb_Density
#' @title  pb_Density
#'
#' @description  This function integrates the probability density function for raster based on provided resolution, radii, and weights. Generally this function is not useful by itself, but is used as a helper function.
#' 
#' @returns Returns a list with a matrix of the weights from the pdf, the input resolution, and the input boundaries
#'
#'
#' @param  cellMeters Resolution of the probability buffer. Higher resolutions will result in more accuracy, but will take more computational resources to use. Default is 50 meters.
#' @param boundaries A numeric vector of the radius limits. Defaults to 5,000 and 10,000, which is the displacement for rural clusters in the DHS (99% displaced between 0 and 5 km, 1% displaced between 0 and 10 km)
#' @param weights A numeric vector for weights of the boundaries set in the boundaries argument, if using more than one boundary. Defaults to c(0.99, 0.01) based on the rural displacement probability
#'
#' @author J.W. Rozelle
#'
#'
#' @export pb_Density
#' @examples
#'
#' # create the probability density matrix
#' ruralPB <- pb_Density()
#' 
#' ruralPB.matrix <- ruralPB[["weightedCircle"]]
#' 
#'


pb_Density <- function(cellMeters = 50, boundaries = c(5e3, 1e4), weights = c(0.99, 0.01)) {
  
  warning("Warning: This is depreciated. Use pb_integratedDensity instead.")
  
  if (length(boundaries) == 1) {
    weights <- 1
  }
  
  # errors
  if (length(boundaries) != length(weights)) {
    stop("The length of the boundaries and weights vectors must be equal")
  }
  
  if (sum(weights)!= 1) {
    stop("The weights should be a numeric vector that sums to one")
  }
  
  probCircs <- list()
  
  for (i in 1:length(boundaries)) {
    # do the weight calculation
    probCirc <- weightCalc(cellMeters, boundaries[i])
    
    # extract the weighted circle
    probCircs[[i]] <- probCirc$weightedCircle
    
    # replace missing values with 0
    probCircs[[i]][is.na(probCircs[[i]][])] <- 0
    
    # Multiple by the appropriate weight
    probCircs[[i]][] <- probCircs[[i]][]*weights[i]
    
    if (i > 1) {
      # need to add in empty values for the inner boundary to align with the outer boundary
      #   First get the number of rows
      rowColAdd <- nrow(probCircs[[i]]) - nrow(probCircs[[i-1]])
      #   create top block of zeroes and bottom block of zeroes
      topBotBlock <- matrix(0, nrow = rowColAdd/2, ncol = nrow(probCircs[[i-1]]))
      sideBlock <- matrix(0, nrow = nrow(probCircs[[i]]), ncol = rowColAdd/2)
      
      probCircInner <- rbind(topBotBlock, probCircs[[i-1]])
      probCircInner <- rbind(probCircInner, topBotBlock)
      probCircInner <- cbind(probCircInner, sideBlock)
      probCircInner <- cbind(sideBlock, probCircInner)
      
      probCircs[[i]] <- probCircInner + probCircs[[i]]
      
      probCircs[[i]][probCircs[[i]]==0] <- NA
    }
    
  }
  
  return(list(weightedCircle = probCircs[[length(boundaries)]], cellMeters = cellMeters, radiusMeters = boundaries[length(boundaries)]))
  
  
}  




#' @name pb_integratedDensity
#' @rdname pb_integratedDensity
#' @title  pb_integratedDensity
#'
#' @description  This function integrates the probability density function for raster based on provided resolution, radii, and weights. Generally this function is not useful by itself, but is used as a helper function.
#' 
#' @returns Returns a list with a matrix of the weights from the pdf, the input resolution, and the input boundaries
#'
#'
#' @param  cellMeters Resolution of the probability buffer. Higher resolutions will result in more accuracy, but will take more computational resources to use. Default is 50 meters.
#' @param boundaries A numeric vector of the radius limits. Defaults to 5,000 and 10,000, which is the displacement for rural clusters in the DHS (99% displaced between 0 and 5 km, 1% displaced between 0 and 10 km)
#' @param weights A numeric vector for weights of the boundaries set in the boundaries argument, if using more than one boundary. Defaults to c(0.99, 0.01) based on the rural displacement probability
#'
#' @author J.W. Rozelle
#'
#'
#' @export pb_integratedDensity
#' @examples
#'
#' # create the probability density matrix
#' ruralPB <- pb_integratedDensity()
#' 
#' ruralPB.matrix <- ruralPB[["weightedCircle"]]
#' 
#'




pb_integratedDensity <- function(cellMeters = 50, boundaries = c(5e3, 1e4), weights = c(0.99, 0.01)) {
  
  if (length(boundaries) == 1) {
    weights <- 1
  }
  
  # errors
  if (length(boundaries) != length(weights)) {
    stop("The length of the boundaries and weights vectors must be equal")
  }
  
  if (sum(weights)!= 1) {
    stop("The weights should be a numeric vector that sums to one")
  }
  
  probCircs <- list()
  
  for (i in 1:length(boundaries)) {
    
    # xprime <- 0-0.005
    # yprime <- 5-0.005
    # dimx <- 50/1000
    # dimy <- 50/1000
    
    # define matrix to hold values
    # get number of rows / columns
    quadrantRowColN <- (boundaries[i] / cellMeters)
    
    densMatrix <- matrix(data = 0, nrow = quadrantRowColN*2, ncol = quadrantRowColN*2)
    
    # define the pb function for the given boundary
    buffer_pb_function <- function(x, y) {
      
      result <- ifelse(sqrt(x^2 + y^2) <= boundaries[i], 1/(boundaries[i] * 2 * pi * sqrt(x^2+y^2)), 0)
      
      return(result)
      
    }
    
    # do the weight calculation
    for (q4_row in 1:quadrantRowColN) {
      for (q4_col in 1:quadrantRowColN) {
        
        # Set limits for integrals
        if (q4_row - 1 > 0 | q4_col - 1 > 0) {
          ylimits <- c((q4_row-1)*cellMeters, q4_row*cellMeters)
          xlimits <- c((q4_col-1)*cellMeters, q4_col*cellMeters)
        } else {
          ylimits <- c(0, q4_row*cellMeters)
          xlimits <- c(1, q4_col*cellMeters)
        }
        
        # 
        
        integralValue <- calculus::integral(buffer_pb_function, bounds = list(x = xlimits, y = ylimits), coordinates = "cartesian")
        
        densMatrix[c(quadrantRowColN+q4_row, quadrantRowColN-q4_row+1), c(quadrantRowColN+q4_col, quadrantRowColN-q4_col+1)] <- integralValue$value
        
      }
    }
    
    probCircs[[i]] <- densMatrix
    
    # Multiple by the appropriate weight
    probCircs[[i]][] <- probCircs[[i]][]*weights[i]
    
    # stack the probability circle matrices
    if (i > 1) {
      # need to add in empty values for the inner boundary to align with the outer boundary
      #   First get the number of rows
      rowColAdd <- nrow(probCircs[[i]]) - nrow(probCircs[[i-1]])
      #   create top block of zeroes and bottom block of zeroes
      topBotBlock <- matrix(0, nrow = rowColAdd/2, ncol = nrow(probCircs[[i-1]]))
      sideBlock <- matrix(0, nrow = nrow(probCircs[[i]]), ncol = rowColAdd/2)
      
      probCircInner <- rbind(topBotBlock, probCircs[[i-1]])
      probCircInner <- rbind(probCircInner, topBotBlock)
      probCircInner <- cbind(probCircInner, sideBlock)
      probCircInner <- cbind(sideBlock, probCircInner)
      
      probCircs[[i]] <- probCircInner + probCircs[[i]]
      
      probCircs[[i]][probCircs[[i]]==0] <- NA
    }
    
    
  }
  
  return(list(weightedCircle = probCircs[[length(boundaries)]], cellMeters = cellMeters, radiusMeters = boundaries[length(boundaries)]))
  
  
}









