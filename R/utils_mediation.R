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

#' Check vector for NA or NAN values.
#'
#' `check_missing` rasises and error if the data contains.
#' NA or NAN values.
#'
#' @param x numeric
#'   The data set to check for NA or NAN values.
check_missing <- function(x) {
  
  # grab the name of the variable
  name <- deparse(substitute(x))
  
  if (any(sapply(x, is.nan))) {stop(paste('There are NaN values in', name, sep = ' '))}
  if (any(sapply(x, is.na)))  {stop(paste('There are NA values in', name, sep = ' '))}
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

#' A little helper function that raises an error if the
#' user provides weights, and they are not equal to the
#' number of stopping methods
#'
#' Call this in the `weighted_mediation()` function.
#'
#' @param weights numeric
#'   The weights
#' @param stopping_methods
#'   The stopping method or methods.
#' @param weights_name
#'   The name of the weights we are checking.
#'   Default : 'total_effects_wts'
check_equal_wts_stopping <- function(weights,
                                     stopping_methods,
                                     weights_name = 'total_effects_wts') {
  n_cols_weights <- ncol(weights)
  n_stopping_methods <- length(c(stopping_methods))
  if (!(n_cols_weights == n_stopping_methods)) {
    stop(paste("The number of columns in the", weights_name,
               "must equal the number of stopping methods",
               n_cols_weights, "!=", n_stopping_methods, sep = " "))
  }
}

#' Calculate the actual effects
#' 
#' @param w_11 The Y(1, M(1)) weights
#' @param w_00 The Y(0, M(0)) weights
#' @param w_10 The Y(1, M(0)) weights
#' @param w_01 The Y(0, M(1)) weights
#' @param y_outcome The Y variable
calculate_effects <- function(w_11, w_00, w_10, w_01, y_outcome,sampw=NULL) {
  
  if (is.null(sampw)) {
    sampW <- rep(1, length(y_outcome))
   } 
   else { sampW <- sampw }

  # Calculate the weighted means for each of the conditions
  E_11 <- weighted_mean(y_outcome, w_11*sampW)
  E_00 <- weighted_mean(y_outcome, w_00*sampW)
  E_10 <- weighted_mean(y_outcome, w_10*sampW)
  E_01 <- weighted_mean(y_outcome, w_01*sampW)
  
  # Calculate total, natural direct, and natural indirect effects
  # These are, as follows:
  # - TE = E(Y(1, M(1))) - E(Y(0, M(0)))
  # - NDE(0) = E(Y(1, M(0))) - E(Y(0, M(0)))
  # - NIE(1) = E(Y(1, M(1))) - E(Y(1, M(0)))
  # - NDE(1) = E(Y(1, M(1))) - E(Y(0, M(1)))
  # - NIE(0) = E(Y(0, M(1))) - E(Y(0, M(0)))
  total_effect      <- E_11 - E_00
  natural_direct0   <- E_10 - E_00   # holding the mediator constant at 0
  natural_indirect1 <- E_11 - E_10   # holding the exposure constant at 1
  natural_direct1   <- E_11 - E_01   # holding the mediator constant at 1
  natural_indirect0 <- E_01 - E_00   # holding the exposure constant at 0
  
  res <- data.frame(estimate=c(total_effect,
                               natural_direct0,
                               natural_indirect1,
                               natural_direct1,
                               natural_indirect0,
                               E_00, E_11, E_10, E_01),
                    row.names=c('TE', 'NDE_0', 'NIE_1','NDE_1', 'NIE_0', 
                                'E[Y(0), M(0)]', 'E[Y(1), M(1)]',
                                'E[Y(1), M(0)]', 'E[Y(0), M(1)]'))
  
  return(res)

}
