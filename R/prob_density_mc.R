#' Monte Carlo approximation to pBuffer displacement density
#'
#' Builds a centered probability kernel on a regular grid via Monte Carlo sampling
#' of displacement offsets (dx, dy) under a radial displacement model.
#'
#' Returns a `densityBuffer` list compatible with pBuffer join functions that
#' expect `list(weightedCircle, cellMeters, radiusMeters)`.
#'
#' @param cellMeters Grid cell size in meters (default 50).
#' @param boundaries Numeric vector of max radii (meters). Example: c(5000, 10000).
#' @param weights Mixture weights for `boundaries`. Must sum to 1.
#' @param n_draws Number of Monte Carlo draws (default 2e6).
#' @param uniform_area If FALSE (default), sample radius uniformly on [0, R]
#'   (matches the density used by pb_integratedDensity). If TRUE, sample uniformly
#'   over area (radius ~ sqrt(U)*R).
#' @param seed Optional integer seed for reproducibility.
#' @param na_as_zero If TRUE, zeros are kept as 0; otherwise zeros set to NA
#'   (default matches existing style: FALSE -> set zeros to NA).
#'
#' @return A list with elements: weightedCircle (matrix), cellMeters, radiusMeters.
#' @export
pb_mcDensity <- function(cellMeters = 50,
                         boundaries = c(5e3, 1e4),
                         weights = c(0.99, 0.01),
                         n_draws = 2e6,
                         uniform_area = FALSE,
                         seed = NULL,
                         na_as_zero = FALSE) {
  
  pb_check_pos_scalar(cellMeters, "cellMeters")
  if (!is.numeric(boundaries) || anyNA(boundaries) || any(boundaries <= 0)) {
    stop("boundaries must be a positive numeric vector.", call. = FALSE)
  }
  
  if (length(boundaries) == 1L) {
    weights <- 1
  }
  
  if (length(boundaries) != length(weights)) {
    stop("boundaries and weights must have equal length.", call. = FALSE)
  }
  pb_check_weights_sum1(weights, tol = 1e-8, arg_name = "weights")
  weights <- pb_normalize_weights(weights, "weights")
  
  n <- as.integer(n_draws)
  if (n < 1) stop("n_draws must be >= 1.", call. = FALSE)
  if (!is.null(seed)) set.seed(seed)
  
  # Choose component per draw
  k <- sample.int(length(boundaries), size = n, replace = TRUE, prob = weights)
  R <- boundaries[k]
  
  # Angle and radius
  theta <- stats::runif(n, 0, 2 * pi)
  u <- stats::runif(n, 0, 1)
  r <- if (uniform_area) sqrt(u) * R else u * R
  
  dx <- cos(theta) * r
  dy <- sin(theta) * r
  
  # Build grid: size mirrors pb_integratedDensity() convention (even dimension)
  radiusMeters <- max(boundaries)
  half_cells <- ceiling(radiusMeters / cellMeters)
  dim <- half_cells * 2L
  extent <- dim * cellMeters
  # grid origin at (-extent/2, -extent/2)
  x0 <- -extent / 2
  y0 <- -extent / 2
  
  # Map offsets to cell indices
  col <- floor((dx - x0) / cellMeters) + 1L
  row <- floor((dy - y0) / cellMeters) + 1L
  
  keep <- col >= 1L & col <= dim & row >= 1L & row <= dim
  col <- col[keep]
  row <- row[keep]
  
  # Accumulate counts (or equal weights; counts normalize to probabilities)
  idx <- (col - 1L) * dim + row
  tab <- tabulate(idx, nbins = dim * dim)
  m <- matrix(tab, nrow = dim, ncol = dim)
  
  # Normalize to probability mass 1
  s <- sum(m)
  if (s == 0) stop("Kernel sum is 0; check parameters.", call. = FALSE)
  m <- m / s
  
  if (!na_as_zero) {
    m[m == 0] <- NA_real_
  }
  
  list(weightedCircle = m, cellMeters = cellMeters, radiusMeters = radiusMeters)
}





#' Monte Carlo DHS protocol density kernel (urban or rural)
#'
#' @param urban Logical. If TRUE: max 2km. If FALSE: 5km + 10km tail.
#' @inheritParams pb_mcDensity
#' @param p_rural_hi Fraction on 10km tail for rural (default 0.01).
#'
#' @export
pb_mcDensity_dhs <- function(urban = FALSE,
                             cellMeters = 50,
                             n_draws = 2e6,
                             p_rural_hi = 0.01,
                             uniform_area = FALSE,
                             seed = NULL,
                             na_as_zero = FALSE) {
  
  if (isTRUE(urban)) {
    boundaries <- 2000
    weights <- 1
  } else {
    boundaries <- c(5000, 10000)
    weights <- c(1 - p_rural_hi, p_rural_hi)
  }
  
  pb_mcDensity(
    cellMeters = cellMeters,
    boundaries = boundaries,
    weights = weights,
    n_draws = n_draws,
    uniform_area = uniform_area,
    seed = seed,
    na_as_zero = na_as_zero
  )
}







