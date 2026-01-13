# ======================================================================
# R/utils_densityKernel.R
#
# Purpose:
#   Internal kernel builders used by pb_densityBuffer().
#   - pb_build_kernel_integral(): numeric integration of mass per cell
#   - pb_build_kernel_mc(): MC cloud binned to grid, then normalized
#
# Outputs:
#   A matrix (rows = y, cols = x) of probability mass per cell. The matrix is
#   origin-centered via symmetric mirroring construction.
#
# Assumptions:
#   - kernel = "uniform_radius" matches r ~ Uniform(0, R) (your earlier code).
#   - kernel = "uniform_area"   matches r = sqrt(U)*R.
#
# Notes:
#   The "uniform_radius" density has a 1/r singularity at r=0, but the cell mass
#   is finite. We avoid division-by-zero by adding a small epsilon at the origin.
# ======================================================================

#' @keywords internal
pb_build_kernel_integral <- function(R, cell_m, kernel = c("uniform_radius", "uniform_area")) {
  kernel <- match.arg(kernel)
  
  if (!requireNamespace("calculus", quietly = TRUE)) {
    stop("Package 'calculus' required for integral kernel. install.packages('calculus').", call. = FALSE)
  }
  
  # Half-cell extents for rectangular integration bounds
  half <- cell_m / 2
  
  # Quadrant size (cells) to cover radius
  q_n <- ceiling(R / cell_m) + 2
  q4 <- matrix(0, nrow = q_n, ncol = q_n)
  
  # Define density in Cartesian coordinates
  # uniform_radius: f(x,y) = 1 / (2*pi*R*r) for r<=R
  # uniform_area:   f(x,y) = 1 / (pi*R^2)   for r<=R
  dens_fun <- function(x, y) {
    rr <- sqrt(x^2 + y^2)
    inside <- rr <= R
    
    if (kernel == "uniform_area") {
      out <- ifelse(inside, 1 / (pi * R^2), 0)
      return(out)
    }
    
    # uniform_radius (singular at rr=0)
    eps <- 1e-9
    out <- ifelse(inside, 1 / (2 * pi * R * pmax(rr, eps)), 0)
    out
  }
  
  for (rr in seq_len(q_n)) {
    for (cc in seq_len(q_n)) {
      
      # Cell center coordinates in meters from origin in Q4:
      # (rr-1, cc-1) because (1,1) is origin cell in quadrant coordinates.
      yc <- (rr - 1) * cell_m
      xc <- (cc - 1) * cell_m
      
      ybounds <- c(yc - half, yc + half)
      xbounds <- c(xc - half, xc + half)
      
      integ <- calculus::integral(
        dens_fun,
        bounds = list(x = sort(xbounds), y = sort(ybounds)),
        coordinates = "cartesian"
      )
      
      q4[rr, cc] <- integ$value
    }
  }
  
  # Mirror Q4 to build full matrix centered at origin
  q3 <- q4[, ncol(q4):2, drop = FALSE]
  q2 <- q4[nrow(q4):2, ncol(q4):2, drop = FALSE]
  q1 <- q4[nrow(q4):2, , drop = FALSE]
  
  north <- cbind(q2, q1)
  south <- cbind(q3, q4)
  M <- rbind(north, south)
  
  # Normalize to sum 1 (numerical integration error)
  s <- sum(M, na.rm = TRUE)
  if (is.finite(s) && s > 0) M <- M / s
  
  M
}

#' @keywords internal
pb_build_kernel_mc <- function(R, cell_m, kernel = c("uniform_radius", "uniform_area"), n_draws = 200000) {
  kernel <- match.arg(kernel)
  
  # Draw offsets under chosen kernel
  theta <- stats::runif(n_draws, 0, 2*pi)
  u     <- stats::runif(n_draws, 0, 1)
  
  rad <- if (kernel == "uniform_area") sqrt(u) * R else u * R
  
  x <- cos(theta) * rad
  y <- sin(theta) * rad
  
  # Bin to grid. We create an origin-centered grid that spans [-R, R] plus padding.
  lim <- R + 2 * cell_m
  breaks_x <- seq(-lim, lim, by = cell_m)
  breaks_y <- seq(-lim, lim, by = cell_m)
  
  # Cut returns factor bins; tabulate gives counts per cell
  bx <- cut(x, breaks = breaks_x, include.lowest = TRUE)
  by <- cut(y, breaks = breaks_y, include.lowest = TRUE)
  
  # Build 2D table (rows=y, cols=x)
  tab <- table(by, bx)
  M <- as.matrix(tab)
  
  # Convert counts to probabilities
  s <- sum(M)
  if (s > 0) M <- M / s
  
  M
}
