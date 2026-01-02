# ---- internal validation helpers (not exported) ----

# Check object is sf
pb_check_sf <- function(x, arg_name = "object") {
  if (!inherits(x, "sf")) {
    stop(arg_name, " must be an sf object.", call. = FALSE)
  }
  invisible(TRUE)
}

# Check geometry type (e.g., POINT)
pb_check_geom_type <- function(x, geom = "POINT", arg_name = "object") {
  pb_check_sf(x, arg_name)
  geom_class <- class(sf::st_geometry(x))[1]
  expected <- paste0("sfc_", geom)
  if (geom_class != expected) {
    stop(arg_name, " must have ", geom, " geometry.", call. = FALSE)
  }
  invisible(TRUE)
}

# Check CRS exists
pb_check_crs <- function(x, arg_name = "object") {
  pb_check_sf(x, arg_name)
  if (is.null(sf::st_crs(x))) {
    stop(arg_name, " has no CRS. Set a CRS before proceeding.", call. = FALSE)
  }
  invisible(TRUE)
}

# Check CRS is projected (units in meters)
pb_check_projected <- function(x, arg_name = "object") {
  pb_check_crs(x, arg_name)
  if (sf::st_is_longlat(x)) {
    stop(
      arg_name,
      " is in a geographic (lon/lat) CRS. ",
      "Project to a CRS with meter units before using this function.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

# Check required columns exist
pb_check_cols <- function(x, cols, arg_name = "object") {
  missing <- setdiff(cols, names(x))
  if (length(missing) > 0) {
    stop(
      arg_name,
      " is missing required column(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

# Normalize a numeric weight vector safely
pb_normalize_weights <- function(w, arg_name = "weights") {
  if (!is.numeric(w)) {
    stop(arg_name, " must be numeric.", call. = FALSE)
  }
  s <- sum(w, na.rm = TRUE)
  if (is.na(s) || s <= 0) {
    stop(arg_name, " sum to 0 or NA; cannot normalize.", call. = FALSE)
  }
  w / s
}

# Validate ID column uniqueness
pb_check_unique_id <- function(x, id_col, arg_name = "object") {
  pb_check_cols(x, id_col, arg_name)
  if (anyDuplicated(x[[id_col]]) > 0) {
    stop(
      arg_name,
      " column '", id_col, "' must be unique.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}
