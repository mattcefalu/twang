weighted_mediation2 <- function(a_treatment,
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
                                ps_version = "legacy",
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
                  "to estimate Model A instead.\n", sep=" "))
    
    ax_covariates <- x_covariates
  }
  
  # Check whether `y_outcome` exists. If it doesn't, we print a
  # warning about outcomes being calculated; if it does, we check to make
  # sure no NA or NaN values exist.
  if (is.null(y_outcome)) {
    warning(paste("The `y_outcome` parameter is NULL. Therefore, only",
                  "weights will be returned; no effects will be calculated.\n",
                  sep=" "))
  } else {
    check_missing(y_outcome)
  }
  
  # Check the number of unique values in `a_treatment`, which must equal 2 (for now)
  unique_treatment_categories <- length(unique(a_treatment))
  if (unique_treatment_categories != 2) {
    stop(paste("The `a_treatment` argument must be ",
               "dichotomous, but this data has",
               unique_treatment_categories,
               "unique values.", sep=" "))
  }
  
  # Default `data_a` and `model_a_res` to `NULL`, since
  # the user may have only provided weights directly for Model A 
  data_a <- NULL
  model_a_res <- NULL
  
  # Get the Model A weights, possibly by estimating the model
  # if `estimate_ax == TRUE`
  if (estimate_ax == TRUE) {
    data_a <- data.frame("A" = (1 - a_treatment), "X" = ax_covariates)
    model_a_res <- do.call(ps, c(list(data = data_a, estimand = "ATE"), ps_args))
  } else if (!is.null(ax_ps)) {
    model_a_res <- ax_ps
  } else if (!is.null(ax_weights)) {
    # If the number of columns in `ax_weights` does not match the
    # number of stopping methods, then we raise an error
    if (ncol(ax_weights) != length(ps_stop.method)) {
      stop(paste("If you are using `ax_weights`, the number of columns must",
                 "equal the length of `ps_stop.method`:",
                 ncol(ax_weights), "!=", length(ps_stop.method), sep=" "))
    }
    # Otherwise, we make sure that the names in `ax_weights` are
    # the same as the stopping methods + estimand type
    names(ax_weights) <- paste(ps_stop.method, "ATE", sep=".")
  }

  # Create model for (1 - A) = X + M
  data_m1 <- data.frame("A" = a_treatment, "X" = x_covariates, "M" = m_mediator)
  data_m2 <- data.frame("A" = (1 - a_treatment), "X" = x_covariates, "M" = m_mediator)

  model_m1_res <- do.call(ps, c(list(data = data_m1, estimand = "ATT"), ps_args))
  model_m2_res <- do.call(ps, c(list(data = data_m2, estimand = "ATT"), ps_args))

  # get the weights for Model A, and the natural indirect weights
  model_m1_weights <- model_m1_res$w
  model_m2_weights <- model_m2_res$w
  model_a_weights <- if (!is.null(ax_weights)) ax_weights else model_a_res$w

  # the below gives us E(Y(a, M(a'))) weighting, using the equation
  # [p(A = a' | M = m, X = x) /  p(A = a | M = m, X = x)] *  [1 / p(A = a' | X = x)]

  # we extract the natural indirect weights for A == 1;
  # we multply by `1 / ((aw - 1) / aw)`, which is just Model A `ps`,
  # since we may *only* have A weights from the user
  natural_indirect_weights <- ((model_m2_res$ps /
                                model_m1_res$ps) *
                               (1 / ((model_a_weights - 1) / model_a_weights)))

  # we extrat the natural indirect weights for A == 0
  # we multply by `1 + (1 / (aw - 1))`, which is just Model A `ps`,
  # since we may *only* have A weights from the user
  mask <- which(data_m1[['A']] == 0)
  natural_indirect_weights[mask, ] <- (((1 - model_m2_res$ps[mask, ]) /
                                        (1 - model_m1_res$ps[mask, ])) *
                                       (1 + (1 / (model_a_weights[mask, ] - 1))))
  
  # collect the results into a list; if we don't have a model for
  # model A, then we just make this NULL; we assign the entire
  # list the class `ps_mediation
  results <- list(model_a_prime = model_a_res,
                  model_m = model_m1_res,
                  model_m_prime = model_m2_res,
                  model_a_prime_weights = model_a_weights,
                  model_m_weights = model_m1_weights,
                  model_m_prime_weights = model_m2_weights,
                  natural_indirect_weights = natural_indirect_weights,
                  data = data_m1,
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
    stopping_method_a <- paste(stopping_method, "ATE", sep=".")
    stopping_method_nie <- paste(stopping_method, "ATT", sep=".")
    
    # Grab the weights from the `ps` object
    weights_a <- results$model_a_prime_weights[, stopping_method_a]
    weights_nie <- results$natural_indirect_weights[, stopping_method_nie]
    
    # Calculate the weighted means for each of the conditions
    treatment1_mediator1 <- weighted_mean(y_outcome, weights_a, a_treatment == 1)
    treatment0_mediator0 <- weighted_mean(y_outcome, weights_a, a_treatment == 0)
    treatment1_mediator0 <- weighted_mean(y_outcome, weights_nie, a_treatment == 1)
    
    # Calculate overall, natural direct, and natural indirect effects
    overall_effect <- treatment1_mediator1 - treatment0_mediator0
    natural_direct <- treatment1_mediator0 - treatment0_mediator0
    natural_indirect <- treatment1_mediator1 - treatment1_mediator0
    
    # Collect the results for this stopping method, and add
    # them back into the original results object
    effects_name = paste(stopping_method, "effects", sep = "_")
    results[[effects_name]] <- list("overall_effect" = overall_effect,
                                    "natural_direct_effect" = natural_direct,
                                    "natural_indirect_effect" = natural_indirect,
                                    "expected_treatment0_mediator0" = treatment0_mediator0,
                                    "expected_treatment1_mediator1" = treatment1_mediator1,
                                    "expected_treatment1_mediator0" = treatment1_mediator0)
    
  }
  return(results)
}
