# ======================================================================
# join_raster.R
#
# Purpose:
#   Probability-weighted raster joins under DHS-style positional uncertainty.
#
# Key changes in this replacement:
#   - pb_rasterValJoin() now accepts a pb_densityBuffer object (from pb_densityBuffer(),
#     or legacy pb_Density/pb_integratedDensity outputs which now return that class).
#   - Mixture support is handled internally by pb_apply_densityBuffer() and therefore
#     works automatically for all future join types.
#   - pb_quickRasterValJoin() is retained for backwards compatibility as an alias-style
#     wrapper; it uses the same implementation path.
#
# Notes:
#   - "quick" vs "non-quick" historically differed in how weights were computed.
#     With densityBuffer caching, the distinction largely disappears, but we keep
#     both functions to preserve your public API.
# ======================================================================

#' Probability-weighted raster value under positional uncertainty
#'
#' @description
#' For each displaced point, estimates the raster value at the unobserved true location
#' by applying a probability kernel (density buffer) centered at the displaced point and
#' taking a probability-weighted mean of raster cell values.
#'
#' @details
#' The kernel is provided by \code{densityBuffer}, created via \code{pb_densityBuffer()}.
#' The densityBuffer may be integral-based or MC-based, and may optionally include mixture
#' metadata (e.g., rural 99% 5km + 1% 10km). Mixture logic and optional admin trimming are
#' applied internally and transparently.
#'
#' Missing raster values are excluded from the weighted mean; weights are implicitly
#' renormalized after trimming by \code{pb_apply_densityBuffer()}.
#'
#' @param displaced.sf An \code{sf} POINT object with displaced locations.
#' @param inputRaster A \code{raster::RasterLayer} to extract values from.
#' @param densityBuffer A \code{pb_densityBuffer} object (or legacy object returned by
#'   \code{pb_Density()} / \code{pb_integratedDensity()}).
#' @param displaced.id Name of unique ID column in \code{displaced.sf}.
#' @param adminBound Optional \code{sf} polygons used to constrain weights to the polygon
#'   containing each point.
#' @param adminID Admin polygon ID column name in \code{adminBound} (default "ID_2").
#' @param radius_m Optional scalar radius override (single-kernel mode). If NULL,
#'   uses the largest radius available in densityBuffer.
#' @param n.cores Number of parallel workers (default 1). Uses \code{future.apply} if > 1.
#'
#' @return A list with:
#'   \itemize{
#'     \item \code{estimate.df}: data.frame with \code{displaced.id} and \code{estimate}
#'     \item \code{commVals.df}: data.frame placeholder (kept for API compatibility)
#'   }
#'
#' @export
pb_rasterValJoin <- function(displaced.sf,
                             inputRaster,
                             densityBuffer,
                             displaced.id = "DHSID",
                             adminBound = NULL,
                             adminID = "ID_2",
                             radius_m = NULL,
                             n.cores = 1) {
  
  pb_check_sf(displaced.sf, "displaced.sf")
  pb_check_cols(displaced.sf, displaced.id, "displaced.sf")
  pb_check_unique_id(sf::st_drop_geometry(displaced.sf), displaced.id, "displaced.sf")
  
  if (!inherits(inputRaster, "RasterLayer")) {
    stop("inputRaster must be a raster::RasterLayer", call. = FALSE)
  }
  
  # Accept legacy density buffer objects (now class pb_densityBuffer from our wrappers)
  if (!inherits(densityBuffer, "pb_densityBuffer")) {
    stop("densityBuffer must be a pb_densityBuffer object created by pb_densityBuffer()/pb_Density()/pb_integratedDensity().",
         call. = FALSE)
  }
  
  if (!is.null(adminBound)) {
    pb_check_sf(adminBound, "adminBound")
    pb_check_cols(adminBound, adminID, "adminBound")
    pb_check_unique_id(sf::st_drop_geometry(adminBound), adminID, "adminBound")
    
    if (sf::st_crs(adminBound) != sf::st_crs(displaced.sf)) {
      stop("adminBound and displaced.sf must share CRS.", call. = FALSE)
    }
  }
  
  # Per-point estimator (thin): build weights raster, then weighted mean of raster values
  one_point <- function(i) {
    
    pt <- displaced.sf[i, , drop = FALSE]
    
    w_r <- pb_apply_densityBuffer(
      point_sf = pt,
      densityBuffer = densityBuffer,
      templateRaster = inputRaster,
      radius_m = radius_m,
      adminBound = adminBound,
      adminID = adminID
    )
    
    if (is.null(w_r)) return(NA_real_)
    
    # Align the input raster to the weights raster locally and compute weighted mean
    r_loc <- suppressWarnings(raster::crop(inputRaster, w_r))
    r_loc <- suppressWarnings(raster::extend(r_loc, w_r))
    w_r   <- suppressWarnings(raster::extend(w_r, r_loc))
    
    v <- raster::values(r_loc)
    w <- raster::values(w_r)
    
    ok <- !is.na(v) & !is.na(w)
    if (!any(ok)) return(NA_real_)
    
    stats::weighted.mean(v[ok], w = w[ok])
  }
  
  # Parallel option (future.apply)
  if (n.cores > 1) {
    if (!requireNamespace("future", quietly = TRUE) || !requireNamespace("future.apply", quietly = TRUE)) {
      stop("For n.cores > 1, install 'future' and 'future.apply'.", call. = FALSE)
    }
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    future::plan(future::multisession, workers = n.cores)
    
    est <- future.apply::future_sapply(seq_len(nrow(displaced.sf)), one_point, future.seed = TRUE)
  } else {
    est <- vapply(seq_len(nrow(displaced.sf)), one_point, numeric(1))
  }
  
  estimate.df <- data.frame(
    stringsAsFactors = FALSE,
    check.names = FALSE,
    displaced_id_value = as.character(sf::st_drop_geometry(displaced.sf)[[displaced.id]]),
    estimate = as.numeric(est)
  )
  names(estimate.df)[names(estimate.df) == "displaced_id_value"] <- displaced.id
  
  # API compatibility placeholder: you can later populate with diagnostics if desired
  commVals.df <- data.frame()
  
  list(estimate.df = estimate.df, commVals.df = commVals.df)
}

#' Quick probability-weighted raster join (API-compatible wrapper)
#'
#' @description
#' Historically a faster variant of \code{pb_rasterValJoin()}. With densityBuffer
#' caching, \code{pb_rasterValJoin()} is already the "quick" path. This function
#' is retained as a wrapper for backwards compatibility.
#'
#' @inheritParams pb_rasterValJoin
#' @export
pb_quickRasterValJoin <- function(displaced.sf,
                                  inputRaster,
                                  densityBuffer,
                                  displaced.id = "DHSID",
                                  adminBound = NULL,
                                  adminID = "ID_2",
                                  radius_m = NULL,
                                  n.cores = 1) {
  
  pb_rasterValJoin(
    displaced.sf = displaced.sf,
    inputRaster = inputRaster,
    densityBuffer = densityBuffer,
    displaced.id = displaced.id,
    adminBound = adminBound,
    adminID = adminID,
    radius_m = radius_m,
    n.cores = n.cores
  )
}
