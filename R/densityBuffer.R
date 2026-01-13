# ======================================================================
# R/densityBuffer.R
#
# Purpose:
#   Define a reusable "density buffer" object representing the DHS-style
#   displacement PDF (or an approximation of it) on a grid.
#
# Key idea:
#   A densityBuffer is origin-centered and contains kernel matrices for one
#   or more radii, built either by:
#     - method = "integral" : deterministic numerical integration over cells
#     - method = "mc"       : Monte Carlo cloud condensed to a probability grid
#
# Why this exists:
#   - Enables a unified workflow across join types: raster, polygon, point.
#   - Separates the expensive kernel construction from repeated joins.
#
# Notes:
#   - The kernel density has a point singularity at r=0 for uniform-in-radius.
#     This is integrable; cell masses are finite. Integral method uses a small
#     epsilon for numerical stability; MC method avoids r=0 almost surely.
# ======================================================================

#' Create a DHS-style density buffer object (integral or MC)
#'
#' @description
#' Builds an origin-centered probability kernel on a regular grid that can be
#' re-used for probability-weighted joins.
#'
#' @param method One of "integral" or "mc".
#' @param radii_m Numeric vector of radii (meters) to build kernels for.
#' @param cell_m Grid cell size (meters) for the kernel grid.
#' @param kernel One of "uniform_radius" or "uniform_area".
#' @param urban_col Name of urban/rural column (default "URBAN_RURA").
#' @param urban_value Value indicating urban (default "U").
#' @param rural_value Value indicating rural (default "R").
#' @param p_rural_hi Tail probability (used only when you want mixture kernels).
#' @param rural_max_hi Tail radius (meters) for rural tail.
#' @param n_draws_mc Number of MC draws for kernel approximation (method="mc").
#' @param seed Optional seed for reproducibility (method="mc").
#'
#' @return An object of class "pb_densityBuffer".
#' @export
pb_densityBuffer <- function(method = c("integral", "mc"),
                             radii_m,
                             cell_m = 100,
                             kernel = c("uniform_radius", "uniform_area"),
                             urban_col = "URBAN_RURA",
                             urban_value = "U",
                             rural_value = "R",
                             # --- mixture support (optional) ---
                             mixture = FALSE,
                             mixture_weights = c(0.99, 0.01),
                             mixture_radii_m = c(5000, 10000),
                             mixture_applies_to = c("rural_only", "global"),
                             # MC options
                             n_draws_mc = 200000,
                             seed = NULL) {
  
  method <- match.arg(method)
  kernel <- match.arg(kernel)
  mixture_applies_to <- match.arg(mixture_applies_to)
  
  # Ensure radii are numeric and positive
  radii_m <- as.numeric(radii_m)
  
  if (isTRUE(mixture)) {
    # Validate mixture vectors
    if (length(mixture_weights) != length(mixture_radii_m)) {
      stop("mixture_weights and mixture_radii_m must have the same length.", call. = FALSE)
    }
    if (any(mixture_weights < 0) || abs(sum(mixture_weights) - 1) > 1e-8) {
      stop("mixture_weights must be nonnegative and sum to 1.", call. = FALSE)
    }
    mixture_radii_m <- as.numeric(mixture_radii_m)
    if (any(!is.finite(mixture_radii_m)) || any(mixture_radii_m <= 0)) {
      stop("mixture_radii_m must be positive finite meters.", call. = FALSE)
    }
    
    # Ensure kernels exist for mixture radii too
    radii_all <- sort(unique(c(radii_m, mixture_radii_m)))
  } else {
    radii_all <- sort(unique(radii_m))
  }
  
  if (any(!is.finite(radii_all)) || any(radii_all <= 0)) {
    stop("radii_m must be positive finite meters.", call. = FALSE)
  }
  if (!is.finite(cell_m) || cell_m <= 0) {
    stop("cell_m must be positive.", call. = FALSE)
  }
  
  if (!is.null(seed)) set.seed(seed)
  
  kernels <- lapply(radii_all, function(R) {
    if (method == "integral") {
      pb_build_kernel_integral(R = R, cell_m = cell_m, kernel = kernel)
    } else {
      pb_build_kernel_mc(R = R, cell_m = cell_m, kernel = kernel, n_draws = n_draws_mc)
    }
  })
  names(kernels) <- as.character(radii_all)
  
  mix_obj <- NULL
  if (isTRUE(mixture)) {
    mix_obj <- list(
      enabled = TRUE,
      applies_to = mixture_applies_to,
      components = lapply(seq_along(mixture_weights), function(i) {
        list(
          name = paste0("comp", i),
          weight = mixture_weights[i],
          radius_m = mixture_radii_m[i]
        )
      })
    )
  }
  
  obj <- list(
    method = method,
    kernel = kernel,
    radii_m = radii_all,
    cell_m = cell_m,
    urban_col = urban_col,
    urban_value = urban_value,
    rural_value = rural_value,
    mixture = mix_obj,
    kernels = kernels
  )
  
  class(obj) <- "pb_densityBuffer"
  obj
}

