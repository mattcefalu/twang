#' Check vector for NA or NAN values.
#'
#' `check_missing` rasises and error if the data contains.
#' NA or NAN values.
#'
#' @param x numeric
#'   The data set to check for NA or NAN values.
#'
#' @export
check_missing <- function(x) {
  
  # grab the name of the variable
  name <- deparse(substitute(x))
  
  if (any(sapply(x, is.nan))) {stop(paste('There are NaN values in', name, sep = ' '))}
  if (any(sapply(x, is.na)))  {stop(paste('There are NA values in', name, sep = ' '))}
}

#' Calculate a weighted mean.
#'
#' `weighted_mean` calculates a weighted mean, given a vector.
#'
#' @param x numeric
#'   The the data set
#' @param weights numeric
#'   The weights
#' @param multiplier
#'   An additional vector to multiply
#'   Default : `NULL`
#' @param na.rm
#'   Whether to remove NA values.
#'   Default: `TRUE`
#'
#' @return numeric
#'   The weighted mean of the data.
#'
#' @export
weighted_mean <- function(x, weights, multiplier = NULL, na.rm = TRUE) {
  
  # if multiplier is NULL, then just multiply by 1
  if (is.null(multiplier)){
    multiplier <- rep(1, length(x))
  }
  return(sum(x * multiplier * weights, na.rm = na.rm) / sum(multiplier * weights, na.rm = na.rm))
}
