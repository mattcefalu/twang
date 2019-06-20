#' Weighted mediation analysis.
#'
#' Estimate causal mediation mechanism of a treatment under
#' using propensity score weighting.
#'
#' For users comfortable with [ps], any options prefaced with
#' `ps_` are passed directly to the `ps()` function.
#'
#' @param a_treatment integer
#'   The treatement variable, a, which must be
#'   dichotomous.
#' @param m_mediator numeric
#'   The mediator variable, m.
#' @param x_covariates numeric
#'   The covariates, x
#'   (Default : `NULL`).
#' @param x_covariates_mediator numeric
#'   The covariates, x, if different for the mediator model
#'   (Default : `NULL`).
#' @param y_outcome numeric
#'   The outcome variable, y. If this is not provided, then
#'   no outcomes will be calculated
#'   (Default : `NULL`).
#' @param ps_n.trees integer
#'   Number of gbm iterations passed on to `gbm`
#'   (Default : `10000`).
#' @param ps_interaction.depth integer
#'   Interaction.depth passed on to `gbm`
#'   (Default : `3`).
#' @param ps_shrinkage numerc
#'   Shrinkage passed on to `gbm`
#'   (Default : `0.005`).
#' @param ps_bag.fraction numerc
#'   Bag fraction passed on to `gbm`
#'   (Default : `1.0`).
#' @param ps_perm.test.iters integer
#'   A non-negative integer giving the number of iterations
#'   of the permutation test, passed to `gbm`
#'   (Default : `0`).
#' @param ps_iterlim integer
#'   Maximum number of iterations for the direct optimization
#'   (Default : `1000`).
#' @param ps_stop.method integer
#'   A method or methods of measuring and summarizing balance
#'   across pretreament variables. Current options are `ks.mean`,
#'   `ks.max`, `es.mean`, and `es.max`
#'   (Default : `c("ks.mean", "ks.max")`).
#' @param ps_sampw numeric
#'   Optional sampling weights
#'   (Default : `NULL`).
#' @param ps_obj ps object
#'   Optional `ps` object. If this is provided
#'   it will be used instead of fitting a model
#'   with `x_covariates`.
#'   (Default : `NULL`).
#' @param ps_weights data frame
#'   Optional data.frame object. If this is provided
#'   it will be used instead of fitting a model
#'   with `x_covariates`.
#'   (Default : `NULL`).
#' @param ... list
#'   Additional arguments passed to [twang::ps].
#' @return mediation object
#'   The `mediation` object includes the following:
#'   - `model_a` The model A `ps()` results.
#'   - `model_m` The model M `ps()` results.
#'   - `model_a_weights` The model A weights.
#'   - `model_m_weights` The model M weights.
#'   - `natural_indirect_weights` The natural indirect weights.
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
                               x_covariates = NULL,
                               x_covariates_mediator = NULL,
                               y_outcome = NULL,
                               ps_n.trees = 10000,
                               ps_interaction.depth = 3,
                               ps_shrinkage = 0.005,
                               ps_bag.fraction = 1.0,
                               ps_perm.test.iters = 0,
                               ps_iterlim = 1000,
                               ps_stop.method = c("ks.mean", "ks.max"),
                               ps_sampw = NULL,
                               ps_obj = NULL,
                               ps_weights = NULL,
                               ...) {
  
  # NOTE: For now, we force `verbose = FALSE`
  # and `distribution = 'bernoulli'`; otherwise,
  # users may pass arguments to the `ps()` function
  
  # get the list of initial arguments for `ps()`
  # TODO : Eventually, we will want to make separate sets of
  # ps arguments for Model A and Model M
  ps_args <- list(formula = A ~ .,
                  # distribution = "bernoulli",
                  verbose = FALSE,
                  n.trees = ps_n.trees,
                  interaction.depth = ps_interaction.depth,
                  shrinkage = ps_shrinkage,
                  bag.fraction = ps_bag.fraction,
                  perm.test.iters = ps_perm.test.iters,
                  # iterlim = ps_iterlim,
                  version = 'legacy',
                  stop.method = ps_stop.method,
                  sampw = ps_sampw)
  ps_args = c(ps_args, list(...))
  
  # 1. Run preliminary checks on data --------------------------
  # ------------------------------------------------------------
  
  # first, we check for nulls in treatment and mediator, which must exist.
  check_missing(a_treatment)
  check_missing(m_mediator)
  
  # second, we check to see if `x_covariates` exists. If it doesn't,
  # we make sure `ps_obj` was provided and raise an error otherwise.
  # If it does, we check to make sure no NA or NaN values exist.
  if (is.null(x_covariates)) {
    if (is.null(ps_obj) && is.null(ps_weights)) {
      stop("You must provide either `x_covariates`, `ps_obj`, or `ps_weights`.")
    } else {
      # if the `ps_obj` exists, make sure it's actually a `ps` object
      if (!(class(ps_obj) == 'ps')){
        stop(paste("The `ps_obj` must be a `ps` object, not",
                   class(ps_obj),
                   sep = " "))
      }
    }
  }
  
  # third, we chekc to see if `x_covariates_mediator`  exists. If it doesn't,
  # we see if `x_covariates` exists (and assign that); otherwise, we raise an
  # error. If it does, we check to make sure no NA or NaN values exist.
  if (is.null(x_covariates_mediator)) {
    if (!is.null(x_covariates)) {
      x_covariates_mediator <- x_covariates
    } else {
      stop("You must provide either `x_covariates` or `x_covariates_mediator`.")
    }
  }
  
  # fourth, we check whether `y_outcome` exists. If it doesn't we print a
  # warning about outcomes being calculated; if is does, we check to make
  # sure no NA or NaN values exist.
  if (is.null(y_outcome)) {
    warning(paste("The `y_outcome` parameter is NULL. Therefore, only",
                  "weights will be returned; no outcome will be calculated.",
                  sep=" "))
  } else {
    check_missing(y_outcome)
  }
  
  # fifth, we check the number of unique values in `a_treatment`;
  # if it is not equal to 2, we raise an error, since it must be dichotomous.
  unique_treatment_categories <- length(unique(a_treatment))
  if (unique_treatment_categories != 2) {
    error_message <- paste("The `a.treatment` must be ",
                           "dichotomous, but this data has",
                           unique_treatment_categories,
                           "unique values.", sep=" ")
    stop(error_message)
  }
  
  # 2. Collect data, and fit the models ------------------------
  # ------------------------------------------------------------
  
  # create model for A = X; if `ps_obj` was passed, use that instead
  # and if `ps_weights` were passed, use those instead
  data_a <- NULL
  model_a_res <- NULL
  if (is.null(ps_obj) && is.null(ps_weights)) {
    data_a <- data.frame("A" = a_treatment, "X" = x_covariates)
    model_a_res <- do.call(ps, c(list(data = data_a, estimand = "ATE"), ps_args))
  } else if (!is.null(ps_obj)) {
    model_a_res <- ps_obj
  } else if (!is.null(ps_weights)) {
    if (ncol(ps_weights) != length(ps_stop.method)) {
      stop(paste("If you are using `ps_weights`, the number of columns must",
                 "equal the length of `ps_stop.method`:",
                 ncol(ps_weights), "!=", length(ps_stop.method), sep=" "))
      # we make sure that the names in the `ps_weights` object are
      # the same as the stopping methods + estimand type, as expected
      names(ps_weights) <- paste(ps_stop.method, "ATE", sep=".")
    }
  }
  
  # create model for (1 - A) = X + M
  data_m <- data.frame("A" = (1 - a_treatment), "X" = x_covariates_mediator, "M" = m_mediator)
  model_m_res <- do.call(ps, c(list(data = data_m, estimand = "ATT"), ps_args))
  
  # get the weights for Model A, and the natural indirect weights
  model_m_weights <- model_m_res$w
  model_a_weights <- if (!is.null(ps_weights)) ps_weights else model_a_res$w
  natural_indirect_weights <- model_m_weights * (1 / (1 - (1 / model_a_weights)))
  
  # we set the control group weights to `NA` for now
  model_m_weights[which(a_treatment == 0), ] <- NA
  natural_indirect_weights[which(a_treatment == 0), ] <- NA
  
  # 3. Collect results -----------------------------------------
  # ------------------------------------------------------------
  
  # collect the results into a list; if we don't have a model for
  # model A, then we just make this NULL; we assign the entire
  # list the class `ps_mediation
  results <- list(model_a = model_a_res,
                  model_m = model_m_res,
                  model_a_weights = model_a_weights,
                  model_m_weights = model_m_weights,
                  natural_indirect_weights = natural_indirect_weights,
                  data = data_m,
                  stopping_methods = ps_stop.method,
                  datestamp = date())
  class(results) <- "mediation"
  
  # if there is no `y_outcome`, we just return the results at this point;
  # otherwise, we move on to calculate the actual effects
  if (is.null(y_outcome)) {
    return(results)
  }
  
  # 4. Calculate effects ---------------------------------------
  # ------------------------------------------------------------
  
  # if we have the `y_outcome`, we calculate the effects
  # for each stopping method
  for (stopping_method in c(ps_stop.method)) {
    
    # get the name of the stopping method column
    stopping_method_mlda <- paste(stopping_method, "ATE", sep=".")
    stopping_method_ntrl <- paste(stopping_method, "ATT", sep=".")
    
    # grab the weights from the PS object
    mlda_weights <- results$model_a_weights[, stopping_method_mlda]
    ntrl_weights <- results$natural_indirect_weights[, stopping_method_ntrl]
    
    # calculate the weighted means for each of the following conditions:
    #   - E(Y(1, M(1))); E(Y(0, M(0))); E(Y(1, M(0)))
    treatment1_mediator1 <- weighted_mean(y_outcome, mlda_weights, a_treatment == 1)
    treatment0_mediator0 <- weighted_mean(y_outcome, mlda_weights, a_treatment == 0)
    treatment1_mediator0 <- weighted_mean(y_outcome, ntrl_weights, a_treatment == 1)
    
    # calculate overall, natural direct, and natural indirect effects
    overall_effect <- treatment1_mediator1 - treatment0_mediator0
    natural_direct <- treatment1_mediator0 - treatment0_mediator0
    natural_indirect <- treatment1_mediator1 - treatment1_mediator0
    
    # collect the results for this stopping method, and add
    # them back into the original results object
    effects_name = paste(stopping_method, "effects", sep="_")
    results[[effects_name]] <- list("overall_effect" = overall_effect,
                                    "natural_direct_effect" = natural_direct,
                                    "natural_indirect_effect" = natural_indirect,
                                    "expected_treatment0_mediator0" = treatment0_mediator0,
                                    "expected_treatment1_mediator1" = treatment1_mediator1,
                                    "expected_treatment1_mediator0" = treatment1_mediator0)
    
  }
  return(results)
}