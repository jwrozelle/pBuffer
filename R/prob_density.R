# ======================================================================
# prob_density.R
#
# Purpose:
#   Provide a reusable probability kernel ("density buffer") representing
#   positional uncertainty under DHS-style displacement.
#
# Design:
#   - pb_densityBuffer(): the modern constructor (integral or MC; optional mixture)
#   - pb_Density() and pb_integratedDensity(): backwards-compatible wrappers
#     that return a densityBuffer object with legacy field names used elsewhere
#     in the package.
#
# Notes for future-you:
#   - The "uniform-in-radius" kernel has a 1/r singularity at r=0. This is not
#     a problem for cell-based weights because we integrate over finite-area cells
#     (the mass in the origin cell is finite). We still guard numerically with eps.
#   - Mixture support (e.g., rural 99% 5km + 1% 10km) is supported but optional.
#     When mixture is enabled, kernels are stored for each mixture radius; the
#     *application* of mixture weights happens in pb_apply_densityBuffer().
# ======================================================================

# ---- public constructor ------------------------------------------------

#' Create a reusable DHS-style density buffer object
#'
#' @description
#' Constructs an origin-centered probability kernel on a regular grid. The result
#' is an object consumed by join functions (e.g., \code{pb_rasterValJoin()}), and
#' is intended to be reusable across all future join types (raster/polygon/point).
#'
#' Two construction methods are supported:
#' \itemize{
#'   \item \code{method="integral"}: deterministic numeric integration of cell masses
#'   \item \code{method="mc"}: Monte Carlo approximation binned to a grid
#' }
#'
#' Mixture support is optional (rare): if enabled, kernels are built for each
#' mixture radius and stored in the object; mixture weights are applied later by
#' \code{pb_apply_densityBuffer()} when the kernel is used for a specific point.
#'
#' @param radii_m Numeric vector of radii (meters) to build kernels for.
#' @param cell_m Grid resolution (meters) of the kernel.
#' @param method One of \code{"integral"} or \code{"mc"}.
#' @param kernel One of \code{"uniform_radius"} or \code{"uniform_area"}.
#' @param mixture Logical; if TRUE, store mixture metadata (weights/radii).
#' @param mixture_weights Numeric vector summing to 1 (e.g., c(0.99,0.01)).
#' @param mixture_radii_m Numeric vector of radii for mixture components.
#' @param mixture_applies_to One of \code{"rural_only"} or \code{"global"}.
#' @param n_draws_mc Number of MC draws used to approximate each kernel (method="mc").
#' @param seed Optional seed for reproducibility (method="mc").
#' @param urban_col Name of urban/rural column (default "URBAN_RURA").
#' @param urban_value Value indicating urban (default "U").
#'
#' @return An object of class \code{pb_densityBuffer}.
#' @export
pb_densityBuffer <- function(radii_m,
                             cell_m = 50,
                             method = c("integral", "mc"),
                             kernel = c("uniform_radius", "uniform_area"),
                             mixture = FALSE,
                             mixture_weights = c(0.99, 0.01),
                             mixture_radii_m = c(5000, 10000),
                             mixture_applies_to = c("rural_only", "global"),
                             n_draws_mc = 200000,
                             seed = NULL,
                             urban_col = "URBAN_RURA",
                             urban_value = "U") {
  
  method <- match.arg(method)
  kernel <- match.arg(kernel)
  mixture_applies_to <- match.arg(mixture_applies_to)
  
  radii_m <- as.numeric(radii_m)
  if (any(!is.finite(radii_m)) || any(radii_m <= 0)) {
    stop("radii_m must be positive finite values (meters).", call. = FALSE)
  }
  if (!is.finite(cell_m) || cell_m <= 0) {
    stop("cell_m must be a positive meter resolution.", call. = FALSE)
  }
  
  radii_all <- sort(unique(radii_m))
  
  # If mixture enabled, ensure mixture radii are included in the kernel set
  mix_meta <- NULL
  if (isTRUE(mixture)) {
    if (length(mixture_weights) != length(mixture_radii_m)) {
      stop("mixture_weights and mixture_radii_m must have same length.", call. = FALSE)
    }
    if (any(mixture_weights < 0) || abs(sum(mixture_weights) - 1) > 1e-8) {
      stop("mixture_weights must be nonnegative and sum to 1.", call. = FALSE)
    }
    mixture_radii_m <- as.numeric(mixture_radii_m)
    if (any(!is.finite(mixture_radii_m)) || any(mixture_radii_m <= 0)) {
      stop("mixture_radii_m must be positive finite meters.", call. = FALSE)
    }
    
    radii_all <- sort(unique(c(radii_all, mixture_radii_m)))
    
    mix_meta <- list(
      enabled = TRUE,
      applies_to = mixture_applies_to,
      weights = as.numeric(mixture_weights),
      radii_m = as.numeric(mixture_radii_m)
    )
  }
  
  # Build kernels keyed by radius
  kernels <- lapply(radii_all, function(R) {
    if (method == "integral") {
      pb_build_kernel_integral(R = R, cell_m = cell_m, kernel = kernel)
    } else {
      pb_build_kernel_mc(R = R, cell_m = cell_m, kernel = kernel, n_draws = n_draws_mc, seed = seed)
    }
  })
  names(kernels) <- as.character(radii_all)
  
  obj <- list(
    method = method,
    kernel = kernel,
    cell_m = cell_m,
    radii_m = radii_all,
    kernels = kernels,
    mixture = mix_meta,
    urban_col = urban_col,
    urban_value = urban_value
  )
  class(obj) <- "pb_densityBuffer"
  obj
}

# ---- backwards-compatible wrappers ------------------------------------

#' Backwards-compatible density buffer (legacy API)
#'
#' @description
#' Historical pBuffer function that created a probability circle. This replacement
#' now returns a \code{pb_densityBuffer} object while preserving legacy list fields
#' expected elsewhere in the package.
#'
#' @param cellMeters Grid resolution (meters). Default 50.
#' @param boundaries Radius boundaries (meters). Default c(5000, 10000).
#' @param weights Mixture weights for boundaries. Default c(0.99,0.01).
#'
#' @return A list-like \code{pb_densityBuffer} object with legacy fields:
#'   \code{weightedCircle}, \code{cellMeters}, \code{radiusMeters}.
#' @export
pb_Density <- function(cellMeters = 50,
                       boundaries = c(5000, 10000),
                       weights = c(0.99, 0.01)) {
  
  # Legacy: pb_Density implied mixture. We store mixture metadata and build kernels.
  dens <- pb_densityBuffer(
    radii_m = boundaries,
    cell_m = cellMeters,
    method = "integral",
    kernel = "uniform_radius",
    mixture = TRUE,
    mixture_weights = weights,
    mixture_radii_m = boundaries,
    mixture_applies_to = "global"
  )
  
  # Legacy fields used by older code paths
  dens$weightedCircle <- dens$kernels[[as.character(max(boundaries))]]
  dens$cellMeters <- cellMeters
  dens$radiusMeters <- max(boundaries)
  
  dens
}

#' Backwards-compatible integrated density buffer (legacy API)
#'
#' @description
#' Historically used to create a single-radius probability buffer. This replacement
#' returns a \code{pb_densityBuffer} object but keeps legacy list fields.
#'
#' @param cellMeters Grid resolution (meters). Default 50.
#' @param boundaries Radius boundaries (meters). Default c(5000, 10000).
#' @param weights Mixture weights for boundaries. Default c(0.99,0.01).
#'
#' @return A list-like \code{pb_densityBuffer} object with legacy fields.
#' @export
pb_integratedDensity <- function(cellMeters = 50,
                                 boundaries = c(5000, 10000),
                                 weights = c(0.99, 0.01)) {
  
  # Keep behavior consistent with pb_Density: build kernels for all boundaries,
  # store mixture metadata, and expose the outer-radius matrix as weightedCircle.
  dens <- pb_Density(cellMeters = cellMeters, boundaries = boundaries, weights = weights)
  dens
}
