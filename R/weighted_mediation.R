#' Weighted mediation analysis.
#'
#' Estimate causal mediation mechanism of a treatment
#' using propensity score weighting.
#'
#' For users comfortable with [ps], any options prefaced with
#' `ps_` are passed directly to the `ps()` function.
#'
#' @param a_treatment integer
#'   The treatement variable, a, which must be
#'   dichotomous.
#' @param m_mediator numeric
#'   The mediator variable, m, which can be continuous,
#'   dichotomous, or categorical.
#' @param x_covariates numeric
#'   The covariates, x, for the mediation model, Model M.
#' @param y_outcome numeric, optional
#'   The outcome variable, y. If this is not provided, then
#'   no effects will be calculated and a warning will be raised.
#'   (Default : `NULL`).
#' @param ax_ps : ps object, optional
#'   If `ax_ps` is provided, then the Model A weights will be
#'   extracted from this object.
#'   (Default : `NULL`).
#' @param ax_weights : numeric, optional
#'   If `ax_weights` is provided, then these will be used as
#'   the Model A weights.
#'   (Default : `NULL`).
#' @param estimate_ax logical, optional
#'   If `estimate_ax=TRUE`, then Model A will be estimated,
#'   and the user does not have to provide Model A weights (`ax_weights`)
#'   or a ps object (`ax_ps`) directly. If `ax_covariates` is provided,
#'   then the these covariates will be used to estimate the model p(A = a | X).
#'   Note that if `estimate_ax=TRUE`, Model A will be estimated regardless
#'   of whether `ax_weights` or `ax_ps` were provided.
#'   (Default : `NULL`).
#' @param ax_covariates numeric, optional
#'   The covariates to use when calculating Model A, if `estimate_ax=TRUE`.
#'   If `ax_covariates` is not provided and `estimate_ax=TRUE`, then
#'   `x_covariates` will be used instead, and a warning will be raised.
#'   (Default : `NULL`).
#' @param ps_n.trees integer, optional
#'   Number of gbm iterations passed on to [gbm] or [xgboost] in the [ps]
#'   function.
#'   (Default : `10000`).
#' @param ps_interaction.depth integer, optional
#'   The tree depth passed on to [gbm] or [xgboost] in the [ps]
#'   function.
#'   (Default : `3`).
#' @param ps_shrinkage numerc, optional
#'   The learning rate passed on to [gbm] or [xgboost] in the [ps]
#'   function.
#'   (Default : `0.005`).
#' @param ps_bag.fraction numerc, optional
#'   The fraction of observations randomly selected in each iterations,
#'   passed on to [gbm] or [xgboost] in the [ps] function.
#'   (Default : `1.0`).
#' @param ps_perm.test.iters integer, optional
#'   A non-negative integer giving the number of iterations
#'   of the permutation test for KS, passed to [ps].
#'   (Default : `0`).
#' @param ps_verbose logical, optional 
#'   If `TRUE`, lots of information will be printed to monitor the
#'   the progress of the [ps] fitting. 
#'   (Default: `FALSE`).
#' @param ps_stop.method integer, optional
#'   A method or methods of measuring and summarizing balance
#'   across pretreament variables. Current options are `ks.mean`,
#'   `ks.max`, `es.mean`, and `es.max`
#'   (Default : `c("ks.mean", "ks.max")`).
#' @param ps_version character, optional
#'  Which version of [ps] to use.
#'     * `"legacy"` Uses the prior implementation of the `ps`` function.
#'     * All other options use `ps.fast` implementation.
#'   (Default : "legacy")
#' @param ps_booster character, optional
#'   `gbm` or `xgboost`
#'   (Default : 'gbm').
#' @param ps_sampw numeric, optional
#'   Optional sampling weights
#'   (Default : `NULL`).
#' @param ps_tree_method character, optional
#'   Only used with `xgboost`. See [xgboost] for further details. 
#'   (Default : "hist")
#' @param ps_n.keep  integer, optional
#'   A numeric variable indicating the algorithm should only
#'   consider every `n.keep`-th iteration of the propensity score model and
#'   optimize balance over this set instead of all iterations. Only used
#'   with `xgboost`.
#'   (Default : `1``).
#' @param ps_n.grid integer, optional
#'   A numeric variable that sets the grid size for an initial
#'   search of the region most likely to minimize the `stop.method`. A
#'   value of `n.grid=50` uses a 50 point grid from `1:n.trees`. It
#'   finds the minimum, say at grid point 35. It then looks for the actual
#'   minimum between grid points 34 and 36.Only used with `xgboost`.
#'   (Default: `25``).
#' @param ... list, optional
#'   Additional arguments passed to [ps].
#' @return mediation object
#'   The `mediation` object includes the following:
#'   - `model_a` The model A `ps()` results.
#'   - `model_m` The model M `ps()` results.
#'   - `model_a_weights` The model A weights.
#'   - `model_m_weights` The model M weights.
#'   - `natural_direct_weights` The natural direct weights.
#'   - `data` The data set used to compute model M.
#'   - `stopping_methods` The stopping methods passed to `stop.method`.
#'   - `datestamp` The date when the analysis was run.
#'   - For each `stop.method`, a list with the following:
#'     * `overall_effect` The overall effect.
#'     * `natural_direct_effect` The natural direct effect.
#'     * `natural_indirect_effect` The natural indirect effect.
#'     * `expected_treatment0_mediator0` E(Y(0, M(0)))
#'     * `expected_treatment1_mediator1` E(Y(1, M(1)))
#'     * `expected_treatment1_mediator0` E(Y(0, M(1)))
#'
#' @seealso [ps]
#' @keywords models, multivariate
#'
#' @export
weighted_mediation <- function(a_treatment,
                               m_mediator,
                               x_covariates,
                               y_outcome = NULL,
                               ax_ps = NULL,
                               ax_weights = NULL,
                               estimate_ax = FALSE,
                               ax_covariates = NULL,
                               ps_n.trees = 10000,
                               ps_interaction.depth = 3,
                               ps_shrinkage = 0.005,
                               ps_bag.fraction = 1.0,
                               ps_perm.test.iters = 0,
                               ps_verbose = FALSE,
                               ps_stop.method = c("ks.mean", "ks.max"),
                               ps_version = "fast",
                               ps_booster = "gbm",
                               ps_sampw = NULL,
                               ps_tree_method = "hist",
                               ps_n.keep = 1,
                               ps_n.grid = 25,
                               ...) {

  # Get the list of initial arguments for `ps()`
  ps_args <- list(formula = A ~ .,
                  n.trees = ps_n.trees,
                  interaction.depth = ps_interaction.depth,
                  shrinkage = ps_shrinkage,
                  bag.fraction = ps_bag.fraction,
                  perm.test.iters = ps_perm.test.iters,
                  verbose = ps_verbose,
                  stop.method = ps_stop.method, 
                  version = ps_version,
                  sampw = ps_sampw)

  # If we are not using the legacy version of `ps()`,
  # then there are some additional arguments to pass
  if (ps_version != 'legacy'){
    ps_args <- c(ps_args, list(booster = ps_booster,
                               tree_method = ps_tree_method,
                               n.keep = ps_n.keep,
                               n.grid = ps_n.grid))
  }
  ps_args = c(ps_args, list(...))

  # Check for nulls in treatment and mediator, which must exist
  check_missing(a_treatment)
  check_missing(m_mediator)

  # Check to make sure that either `a_weights` was provided, or `estimate_ax`
  # was set to `TRUE`; otherwise, we raise an error
  if (is.null(ax_ps) && is.null(ax_weights) && estimate_ax == FALSE) {
    stop(paste("If no weights are provided via `ax_weights` and",
               "no `ps` object is provided via `ax_ps,",
               "then users must explicitly set `estimate_ax`",
               "to `TRUE`.", sep = " "))
  }

  # If `ax_ps` was provided, we check to make sure that it is actually a `ps` object
  if (!is.null(ax_ps) && !(class(ax_ps) == 'ps')) {
    stop(paste("The `ax_ps` argument, if provided, must be a",
               "`ps` object, not",
               class(ax_ps),
               sep = " "))
  }

  # If `estimate_ax == TRUE`, we check to see if `ax_covariates`
  # was provided; if not, we raise a warning and set the value equal to `x_covariates`
  if ((estimate_ax == TRUE) && is.null(ax_covariates)) {
    warning(paste("The `estimate_ax` argument was set to `TRUE`,",
                  "but no `ax_covariates` were provided; using `x_covariates`",
                  "to estimate Model A instead.\n", sep = " "))

    ax_covariates <- x_covariates
  }

  # Check whether `y_outcome` exists. If it doesn't, we print a
  # warning about outcomes being calculated; if it does, we check to make
  # sure no NA or NaN values exist.
  if (is.null(y_outcome)) {
    warning(paste("The `y_outcome` parameter is NULL. Therefore, only",
                  "weights will be returned; no effects will be calculated.\n",
                  sep = " "))
  } else {
    check_missing(y_outcome)
  }

  # Check the number of unique values in `a_treatment`, which must equal 2 (for now)
  unique_treatment_categories <- length(unique(a_treatment))
  if (unique_treatment_categories != 2) {
    stop(paste("The `a_treatment` argument must be",
               "dichotomous, but this data has",
               unique_treatment_categories,
               "unique values.", sep = " "))
  }

  # Default `data_a` and `model_a_res` to `NULL`, since
  # the user may have only provided weights directly for Model A 
  data_a <- NULL
  model_a_res <- NULL

  # Get the Model A weights, possibly by estimating the model
  # if `estimate_ax == TRUE`
  if (estimate_ax == TRUE) {
    data_a <- data.frame("A" = a_treatment, "X" = ax_covariates)
    model_a_res <- do.call(ps, c(list(data = data_a, estimand = "ATE"), ps_args))
  } else if (!is.null(ax_ps)) {
    model_a_res <- ax_ps
  } else if (!is.null(ax_weights)) {
    # If the number of columns in `ax_weights` does not match the
    # number of stopping methods, then we raise an error
    if (ncol(ax_weights) != length(ps_stop.method)) {
      stop(paste("If you are using `ax_weights`, the number of columns must",
                 "equal the length of `ps_stop.method`:",
                 ncol(ax_weights), "!=", length(ps_stop.method), sep = " "))
    }
    # Otherwise, we make sure that the names in `ax_weights` are
    # the same as the stopping methods + estimand type
    names(ax_weights) <- paste(ps_stop.method, "ATE", sep = ".")
  }

  # Create model for (1 - A) = X + M
  data_m <- data.frame("A" = (1 - a_treatment), "X" = x_covariates, "M" = m_mediator)
  model_m_res <- do.call(ps, c(list(data = data_m, estimand = "ATT"), ps_args))

  # Swap the treatment and control groups back to their original values
  data_m[['A']] <- (1 - data_m[['A']])

  # Get the weights for Model A, and the natural indirect weights
  model_m_weights <- model_m_res$w
  model_a_weights <- if (!is.null(ax_weights)) ax_weights else model_a_res$w
  natural_direct_weights <- model_m_weights * model_a_weights

  # Collect the results into a list; if we don't have a model for
  # Model A, then this will just be NULL
  results <- list(model_a = model_a_res,
                  model_m = model_m_res,
                  model_a_weights = model_a_weights,
                  model_m_weights = model_m_weights,
                  natural_direct_weights = natural_direct_weights,
                  data = data_m,
                  stopping_methods = ps_stop.method,
                  datestamp = date())
  class(results) <- "mediation"
  class(results) <- "weighted"  

  # If there is no `y_outcome`, we just return the results at this point;
  # otherwise, we move on to calculate the actual effects
  if (is.null(y_outcome)) {
    return(results)
  }

  # If we have the `y_outcome`, we calculate the
  # effects for each stopping method
  for (stopping_method in c(ps_stop.method)) {

    # Get the name of the stopping method column
    stopping_method_a <- paste(stopping_method, "ATE", sep = ".")
    stopping_method_nde <- paste(stopping_method, "ATT", sep = ".")

    # Grab the weights from the `ps` object
    weights_a <- results$model_a_weights[, stopping_method_a]
    weights_nde <- results$natural_direct_weights[, stopping_method_nde]

    # Calculate the weighted means for each of the conditions
    treatment1_mediator1 <- weighted_mean(y_outcome, weights_a, a_treatment == 1)
    treatment0_mediator0 <- weighted_mean(y_outcome, weights_a, a_treatment == 0)
    treatment1_mediator0 <- weighted_mean(y_outcome, weights_nde, a_treatment == 1)

    # Calculate overall, natural direct, and natural indirect effects
    overall_effect <- treatment1_mediator1 - treatment0_mediator0
    natural_direct <- treatment1_mediator1 - treatment1_mediator0
    natural_indirect <- treatment1_mediator0 - treatment0_mediator0

    # Collect the results for this stopping method, and add
    # them back into the original results object
    effects_name = paste(stopping_method, "effects", sep = "_")
    results[[effects_name]] <- list("datestamp" = date(),
                                    "overall_effect" = overall_effect,
                                    "natural_direct_effect" = natural_direct,
                                    "natural_indirect_effect" = natural_indirect,
                                    "expected_treatment0_mediator0" = treatment0_mediator0,
                                    "expected_treatment1_mediator1" = treatment1_mediator1,
                                    "expected_treatment1_mediator0" = treatment1_mediator0)

  }
  return(results)
}