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
#'   no effects will be calculated and a warning will be raised. Default : `NULL`.
#' @param ax_ps : ps object, optional
#'   If `ax_ps` is provided, then the Model A weights will be
#'   extracted from this object. Default : `NULL`.
#' @param ax_weights : numeric, optional
#'   If `ax_weights` is provided, then these will be used as
#'   the Model A weights. Default : `NULL`.
#' @param estimate_ax logical, optional
#'   If `estimate_ax=TRUE`, then Model A will be estimated,
#'   and the user does not have to provide Model A weights (`ax_weights`)
#'   or a ps object (`ax_ps`) directly. If `ax_covariates` is provided,
#'   then the these covariates will be used to estimate the model p(A = a | X).
#'   Note that if `estimate_ax=TRUE`, Model A will be estimated regardless
#'   of whether `ax_weights` or `ax_ps` were provided. Default : `NULL`.
#' @param ax_covariates numeric, optional
#'   The covariates to use when calculating Model A, if `estimate_ax=TRUE`.
#'   If `ax_covariates` is not provided and `estimate_ax=TRUE`, then
#'   `x_covariates` will be used instead, and a warning will be raised. Default : `NULL`.
#' @param ps_n.trees integer, optional
#'   Number of gbm iterations passed on to [gbm]. Default: 10000.
#' @param ps_interaction.depth integer, optional
#'   A positive integer denoting the tree depth used in
#'   gradient boosting. Default: 3.
#' @param ps_shrinkage numerc, optional
#'   A numeric value between 0 and 1 denoting the learning rate.
#'   See [gbm] for more details. Default: 0.01.
#' @param ps_bag.fraction numerc, optional
#'   A numeric value between 0 and 1 denoting the fraction of
#'   the observations randomly selected in each iteration of the gradient
#'   boosting algorithm to propose the next tree. See [gbm] for
#'   more details. Default: 1.0.
#' @param ps_n.minobsinnode An integer specifying the minimum number of observations 
#'   in the terminal nodes of the trees used in the gradient boosting.  See [gbm] for
#'   more details. Default: 10.
#' @param ps_perm.test.iters integer, optional
#'   A non-negative integer giving the number of iterations
#'   of the permutation test for the KS statistic. If `perm.test.iters=0`
#'   then the function returns an analytic approximation to the p-value. Setting
#'   `perm.test.iters=200` will yield precision to within 3\% if the true
#'   p-value is 0.05. Use `perm.test.iters=500` to be within 2\%. Default: 0.
#' @param ps_verbose logical, optional 
#'   If `TRUE`, lots of information will be printed to monitor the
#'   the progress of the fitting. Default: `FALSE`.
#' @param ps_stop.method integer, optional
#'   A method or methods of measuring and summarizing balance across pretreatment
#'   variables. Current options are `ks.mean`, `ks.max`, `es.mean`, and `es.max`. `ks` refers to the
#'   Kolmogorov-Smirnov statistic and es refers to standardized effect size. These are summarized
#'   across the pretreatment variables by either the maximum (`.max`) or the mean (`.mean`). 
#'   Default: `c("ks.mean", "ks.max")`.
#' @param ps_version character, optional
#'  "gbm", "xgboost", or "legacy", indicating which version of the twang package to use.
#'   * `"gbm"` uses gradient boosting from the [gbm] package.
#'   * `"xgboost"` uses gradient boosting from the [xgboost] package.
#'   * `"legacy"` uses the prior implementation of the `ps` function.
#' @param ps_ks.exact `NULL` or a logical indicating whether the
#'   Kolmogorov-Smirnov p-value should be based on an approximation of exact
#'   distribution from an unweighted two-sample Kolmogorov-Smirnov test. If
#'   `NULL`, the approximation based on the exact distribution is computed
#'   if the product of the effective sample sizes is less than 10,000.
#'   Otherwise, an approximation based on the asymptotic distribution is used.
#'   **Warning:** setting `ks.exact = TRUE` will add substantial
#'   computation time for larger sample sizes. Default: `NULL`.
#' @param ps_sampw numeric, optional
#'   Optional sampling weights Default : `NULL`.
#' @param ps_n.keep  integer, optional
#'   A numeric variable indicating the algorithm should only
#'   consider every `n.keep`-th iteration of the propensity score model and
#'   optimize balance over this set instead of all iterations. Only used
#'   with `xgboost`. Default : 1.
#' @param ps_n.grid integer, optional
#'   A numeric variable that sets the grid size for an initial
#'   search of the region most likely to minimize the `stop.method`. A
#'   value of `n.grid=50` uses a 50 point grid from `1:n.trees`. It
#'   finds the minimum, say at grid point 35. It then looks for the actual
#'   minimum between grid points 34 and 36.Only used with `xgboost`. Default: 25.
#' @return mediation object
#'   The `mediation` object includes the following:
#'   - `model_a_res` The model A `ps()` results.
#'   - `model_m_res` The model M `ps()` results.
#'   - `model_a_wts` The model A weights.
#'   - `model_m_ets` The model M weights.
#'   - `natural_direct_wts` The natural direct weights.
#'   - `natural_indirect_wts` The natural indirect weights.
#'   - `total_effect_wts` The total effects weights
#'   - `data` The data set used to compute models
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
                               total_effect_wts = NULL,
                               total_effect_ps = NULL,
                               estimate_total_effect_wts = FALSE,
                               ax_covariates = NULL,
                               ps_n.trees = 10000,
                               ps_interaction.depth = 3,
                               ps_shrinkage = 0.01,
                               ps_bag.fraction = 1.0,
                               ps_n.minobsinnode = 10,
                               ps_perm.test.iters = 0,
                               ps_verbose = FALSE,
                               ps_stop.method = c("ks.mean", "ks.max"),
                               ps_version = "gbm",
                               ps_ks.exact = NULL,
                               ps_sampw = NULL,
                               ps_n.keep = 1,
                               ps_n.grid = 25) {
  
  # Get the list of initial arguments for `ps()`
  ps_args <- list(formula = A ~ .,
                  n.trees = ps_n.trees,
                  interaction.depth = ps_interaction.depth,
                  shrinkage = ps_shrinkage,
                  bag.fraction = ps_bag.fraction,
                  n.minobsinnode = ps_n.minobsinnode,
                  perm.test.iters = ps_perm.test.iters,
                  verbose = ps_verbose,
                  stop.method = ps_stop.method, 
                  version = ps_version,
                  sampw = ps_sampw,
                  ks.exact = ps_ks.exact,
                  n.keep = ps_n.keep,
                  n.grid = ps_n.grid)
  
  # Check for nulls in treatment and mediator, which must exist
  check_missing(a_treatment)
  check_missing(m_mediator)
  
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

  # Case 1, the user provides the total effect weights, W_{A=0|X_A}, directly
  if (!is.null(total_effect_wts)) {

    # The total effects weights length must equal the number of stopping methods
    stopifnot(ncol(total_effect_wts) == length(c(ps_stop.method)))
    
    # We still need to estimate W_{A=0|X_M} because we don't know if X_A == X_M
    data_a = data.frame("A" = a_treatment, "X" = x_covariates)
    model_a_res <- do.call(ps, c(list(data = data_a, estimand = "ATE"), ps_args))
    model_a_wts <- calculate_weights(model_a_res$ps, a_treatment, use_opposite = TRUE)
    
  } else {
    
    # Case 2, the user provides the `ps` object directly, which contains
    # the total effect weights, W_{A=0|X_A}
    if (!is.null(total_effect_ps)) {
      
      # The total effects weights length must equal the number of stopping methods
      stopifnot(ncol(total_effect_ps$w) == length(c(ps_stop.method)))
      
      # We want to check if X_A is a subset of X_M and vice versa.
      # This will raise an error if there are elements in X_A that are not in X_M
      a_equals_m <- check_subset_equal(total_effect_ps, x_covariates)
      if (a_equals_m) {
        # The total effect weights W_{A=0|X_A} and W_{A=0|X_M} are equivalent
        total_effect_wts <- calculate_weights(total_effect_ps$ps, a_treatment)
        model_a_res <- total_effect_ps
        model_a_wts <- calculate_weights(total_effect_ps$ps, a_treatment, use_opposite = TRUE)
        
        # If X_A is a subset of X_M, then we need to estimate P(A=0|X_M)
      } else {
        # The total effect weights are 1 / (1 - total_effect_ps)
        total_effect_wts <- calculate_weights(total_effect_ps$ps, a_treatment)
        
        # We still need to estimate W_{A=0|X_M}
        data_a = data.frame("A" = a_treatment, "X" = x_covariates)
        model_a_res <- do.call(ps, c(list(data = data_a, estimand = "ATE"), ps_args))
        model_a_wts <- calculate_weights(model_a_res$ps, a_treatment, use_opposite = TRUE)

      }

      # Case 3, we need to estimate W_{A=0|X_A} and W_{A=0|X_M} ourselves,
      # but it is possible that they are the same thing
    } else {

      # If estimate_ax is not TRUE, then we raise an error
      if (!isTRUE(estimate_total_effect_wts)) {
        stop(paste('You have not provided `total_effect_wts` or `total_effect_ps`;',
                   'Therefore, you must explicitly set `estimate_total_effect_wts=TRUE`.',
                   sep = ' '))
      }

      # We want to check if X_A is a subset of X_M and vice versa.
      # By default, we assume this is true; we only check if `ax_covariates`
      # we provided. Note that this will raise an error if if there are
      # elements in X_A that are not in X_M
      a_equals_m <- TRUE
      if (!is.null(ax_covariates)) {
        a_equals_m <- check_subset_equal(ax_covariates, x_covariates)
      } else {
        # Since `ax_covariates` were not provided, we use `x_covariates`.
        # Implicitly, this means that X_A will equal X_M.
        ax_covariates <- x_covariates
      }

      # We need to estimate W_{A=0|X_M} no matter what
      data_a = data.frame("A" = a_treatment, "X" = x_covariates)
      model_a_res <- do.call(ps, c(list(data = data_a, estimand = "ATE"), ps_args))
      model_a_wts <- calculate_weights(model_a_res$ps, a_treatment, use_opposite = TRUE)
      
      if (!a_equals_m) {
        # We need to estimate W_{A=0|X_A} because X_A != X_M
        data_ax = data.frame("A" = a_treatment, "X" = ax_covariates)
        model_ax_res <- do.call(ps, c(list(data = data_ax, estimand = "ATE"), ps_args))
        total_effect_wts <- calculate_weights(model_ax_res$ps, a_treatment)
      
      } else {
        # We don't need to estimate W_{A=0|X_A} because X_A = X_M
        total_effect_wts <- calculate_weights(model_a_res$ps, a_treatment)
      }
    }
  }

  # Esimate models for (1 - A) = X + M and A = X + M; together, 
  # these will give us the model M weights of interest: 
  #  - P(A = 0 | M, X_M) /  P(A = 1 | M, X_M)
  #  - P(A = 1 | M, X_M) /  P(A = 0 | M, X_M)
  data_m0 <- data.frame("A" = (1 - a_treatment), "X" = x_covariates, "M" = m_mediator)
  data_m1 <- data.frame("A" = a_treatment, "X" = x_covariates, "M" = m_mediator)
  model_m0_res <- do.call(ps, c(list(data = data_m0, estimand = "ATT"), ps_args))
  model_m1_res <- do.call(ps, c(list(data = data_m1, estimand = "ATT"), ps_args)) 

  # Grab the indexes for control group  
  ctrl <- which(a_treatment == 0)

  # Calculate the model M weights
  model_m_wts <- (model_m0_res$ps / model_m1_res$ps)
  model_m_wts[ctrl,] <- (model_m1_res$ps[ctrl,] / model_m0_res$ps[ctrl,])

  # Calculate the NDE weights
  natural_direct_wts <- model_m_wts * model_a_wts

  # Calculate the NIE weights (which are just TE weights - NDE weights)
  natural_indirect_wts <- total_effect_wts - natural_direct_wts

  # Let's remove the data from these `ps` objects, so they're not so bulky
  model_m0_res[['data']] <- NULL
  model_m1_res[['data']] <- NULL
  model_a_res[['data']] <- NULL

  # Collect the results into a list
  results <- list(natural_indirect_wts = natural_indirect_wts,
                  natural_direct_wts = natural_direct_wts,
                  total_effect_wts = total_effect_wts,
                  model_m_wts = model_m_wts,
                  model_m0 = model_m0_res,
                  model_m1 = model_m1_res,
                  model_a_wts = model_a_wts,
                  model_a = model_a_res,
                  stopping_methods = ps_stop.method,
                  data = data_m1,
                  datestamp = date())
  class(results) <- "mediation"

  # If there is no `y_outcome`, we just return the results at this point;
  # otherwise, we move on to calculate the actual effects
  if (is.null(y_outcome)) {
    return(results)
  }

  # If we have the `y_outcome`, we calculate the
  # effects for each stopping method
  stop_methods <- c(ps_stop.method)
  for (idx in 1:length(stop_methods)) {

    stop_method <- stop_methods[idx]

    # Grab the weights from the `ps` object
    te_wts <- results$total_effect_wts[, idx]
    nd_wts <- results$natural_direct_wts[, idx]

    # Calculate the weighted means for each of the conditions
    treatment1_mediator1 <- weighted_mean(y_outcome, te_wts, a_treatment == 1)
    treatment0_mediator0 <- weighted_mean(y_outcome, te_wts, a_treatment == 0)
    treatment1_mediator0 <- weighted_mean(y_outcome, nd_wts, a_treatment == 1)
    treatment0_mediator1 <- weighted_mean(y_outcome, nd_wts, a_treatment == 0)

    # Calculate overall, natural direct, and natural indirect effects
    total_effect <- treatment1_mediator1 - treatment0_mediator0
    natural_direct <- treatment1_mediator0 - treatment0_mediator0
    natural_indirect <- treatment1_mediator1 - treatment1_mediator0

    # Collect the results for this stopping method, and add
    # them back into the original results object
    effects_name = paste(stop_method, "effects", sep = "_")
    results[[effects_name]] <- list(total_effect = total_effect,
                                    natural_direct_effect = natural_direct,
                                    natural_indirect_effect = natural_indirect,
                                    expected_treatment0_mediator0 = treatment0_mediator0,
                                    expected_treatment1_mediator1 = treatment1_mediator1,
                                    expected_treatment1_mediator0 = treatment1_mediator0,
                                    expected_treatment0_mediator1 = treatment0_mediator1)
    
  }
  return(results)
}