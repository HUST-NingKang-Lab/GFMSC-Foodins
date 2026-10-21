# Shared input checks for the Kefir analysis scripts.
require_columns <- function(x, cols) {
  missing <- setdiff(cols, names(x))
  if(length(missing)) stop("Missing columns: ", paste(missing, collapse=", "))
}
check_id <- function(x, col) {
  require_columns(x, col)
  if(anyNA(x[[col]]) || any(!nzchar(x[[col]])) || anyDuplicated(x[[col]]))
    stop(col, " must be nonempty and unique")
}
check_numeric <- function(x, cols) {
  if(!length(cols)) stop("No numeric measurement columns")
  for(n in cols) if(!is.numeric(x[[n]]) || any(!is.finite(x[[n]])) || any(x[[n]] < 0))
    stop("Expected finite, nonnegative numeric values in: ", n)
}
read_mapping <- function(path) {
  x <- unique(readr::read_csv(path, show_col_types=FALSE)[, c("MAG_id", "SGB_id")])
  check_id(x, "MAG_id")
  if(anyNA(x$SGB_id) || any(!nzchar(x$SGB_id))) stop("Missing SGB mapping")
  x
}
read_summary <- function(path) {
  x <- readr::read_csv(path, show_col_types=FALSE)
  check_id(x, "SGB_id")
  require_columns(x, c("species", "mean_abundance", "prevalence", "n_amp_mag", "amp_mag_fraction", "sgb_category"))
  check_numeric(x, c("mean_abundance", "prevalence", "n_amp_mag", "amp_mag_fraction"))
  if(any(x$prevalence>1 | x$amp_mag_fraction>1)) stop("Fractions must be in [0,1]")
  x$species <- sub("^s__", "", x$species)
  x
}
read_counts <- function(path) {
  x <- readr::read_csv(path, show_col_types=FALSE)
  check_id(x, "MAG_id")
  check_numeric(x, setdiff(names(x), "MAG_id"))
  x
}
read_samples <- function(path, x) {
  samples <- readLines(path, warn=FALSE)
  if(length(samples)<3 || anyDuplicated(samples) || any(!nzchar(samples)))
    stop("Provide at least 3 unique sample IDs, one per line, no header")
  require_columns(x, samples)
  check_numeric(x, samples)
  samples
}
