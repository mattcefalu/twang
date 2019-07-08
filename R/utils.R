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

#' Check whether the variables in d1 are a subset
#' of the variables in d2. The data sets must be of
#' equal length.
is_subset <- function(d1, d2) {
  stopifnot(dim(d1)[1] == dim(d2)[1])
  dim(setdiff(d1, d2))[2] == 0
}

#' Check whether the Xs in a `ps` object are a 
#' subset of the Xs in a data frame or matrix,
#' and vice versa
#'
#' @export
check_subset_equal <- function(y_vars, x_vars, raise_error = TRUE) {
  
  if (class(y_vars) == 'ps') {
    y_names <- y_vars$gbm.obj$var.names
    y_vars <- y_vars$data[y_names]
  }

  a_subset_m <- is_subset(y_vars, x_vars)
  m_subset_a <- is_subset(x_vars, y_vars)
  a_equal_m <- (a_subset_m == m_subset_a) && (a_subset_m == TRUE) && (m_subset_a == TRUE)

  if (raise_error && !m_subset_a) {
    stop('There are covariates in A that do not appear in M.')
  }
  return(a_equal_m)
  }
