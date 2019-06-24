#' Outcome mediation with [gbm] or [xgboost].
#'
#' @param a_treatment integer
#'   The treatement variable, a, which must be
#'   dichotomous.
#' @param m_mediator numeric
#'   The mediator variable, m.
#' @param x_covariates numeric
#'   The covariates, x
#' @param y_outcome numeric
#'   The outcome variable, y. If this is not provided, then
#'   no outcomes will be calculated
#'   (Default : `NULL`).
#' @param boost_n_trees integer
#'   Number of iterations passed on to [gbm] or [xgboost]
#'   (Default : `10000`).
#' @param boost_n_folds integer
#'   Number of cross-validation folds, passed only to [xgboost]
#'   (Default : `5`).
#' @param boost_maxdepth integer
#'   The maximum interaction depth passed on to [gbm] or [xgboost]
#'   (Default : `3`).
#' @param boost_eta numerc
#'   Shrinkage passed on to [gbm] or [xgboost]
#'   (Default : `0.005`).
#' @param boost_bag_fraction numeric
#'   Bag fraction passed on to [gbm]
#'   (Default : `1.0`).
#' @param boost_train_fraction numeric
#'   Bag training fraction passed only to [gbm]
#'   (Default : `0.9`).
#' @param boost_subsample numeric
#'   The subsample fraction, passed to [xgboost]
#'   (Default : `0.5`)
#' @param boost_early_stopping_rounds integer
#'   The number of rounds of no change after which to
#'   to stop early, passsed to [xgboost]
#'   (Default : `5`)
#' @param booster character
#'   Either `gbm` or `xgboost`.
#'   (Default : `'gbm'`).
#' @param verbose If `TRUE`, lots of information will be printed to monitor
#'   the progress of the fitting.
#'   (Default: `TRUE`.)
#' @param nsims integer
#'   The number of simulations to perform.
#'   (Default : `100`)
#'
#' @export
outcome_mediation <- function(a_treatment,
                              m_mediator,
                              x_covariates,
                              y_outcome,
                              boost_n_trees = 10000,
                              boost_n_folds = 5,
                              boost_maxdepth = 3,
                              boost_eta = 0.005,
                              boost_bag_fraction = 1.0,
                              boost_train_fraction = 0.9,
                              boost_subsample = 0.5,
                              boost_early_stopping_rounds = 5,
                              booster = "gbm",
                              verbose = FALSE,
                              nsims = 100,
                              ...) {

  # get the number of observations
  n_obs = length(y_outcome)

  # make sure booster is either `gbm` or `xgboost`
  if (!(booster %in% c("gbm", "xgboost"))) {
    stop("The `booster` argument must be either 'gbm' or 'xgboost'.")
  }

  # This block is only if we're using `gbm`
  if (booster == "gbm") {

    # figure out the distribution for M
    # TODO : This will need to be updated at some point
    dist_m <- if (length(unique(m_mediator)) == 2) 'bernoulli'  else 'gaussian'
    dist_y <- if (length(unique(y_outcome)) == 2) 'bernoulli'  else 'gaussian'

    # get the list of initial arguments for `gbm()`
    params <- list(verbose = verbose,
                   n.trees = boost_n_trees,
                   cv.folds = boost_n_folds,
                   interaction.depth = boost_maxdepth,
                   shrinkage = boost_eta,
                   train.fraction = boost_train_fraction)
    
    # fit the model P(M|A, X)
    model_m_res <- do.call(gbm::gbm, c(list(formula = M ~ .,
                                            data = data.frame("M" = m_mediator,
                                                              "A" = a_treatment,
                                                              "X" = x_covariates),
                                            distribution = dist_m), params))
    
    # fit the model P(Y|M, A, X)
    model_y_res <- do.call(gbm::gbm, c(list(formula = Y ~ .,
                                            data = data.frame("Y" = y_outcome,
                                                              "M" = m_mediator,
                                                              "A" = a_treatment,
                                                              "X" = x_covariates),
                                            distribution = dist_y), params))
    
    # get the predictions for M1 and M0
    M1 <- predict(model_m_res,
                  newdata = data.frame('X' = x_covariates,
                                       'A' = 1),
                  type = 'response',
                  n.trees = model_m_res$n.trees)
    M0 <- predict(model_m_res,
                  newdata = data.frame('X' = x_covariates,
                                       'A' = 0),
                  type = 'response',
                  n.trees = model_m_res$n.trees)

  # This block is only if we're using `xgboost`
  } else {

    # figure out whether we're using linear or logistic objective
    # functions for both the mediator and the outcome variables
    # TODO : This will need to be updated at some point
    m_objective <- if (length(unique(m_mediator)) == 2) 'reg:logistic'  else 'reg:linear'
    y_objective <- if (length(unique(y_outcome)) == 2) 'reg:logistic'  else 'reg:linear'

    # get the list of initial arguments for `xgboost()`
    m_params <- list(objective = m_objective,
                     eta = boost_eta,
                     maxdepth = boost_maxdepth,
                     subsample = boost_subsample)
    y_params <- list(objective = y_objective,
                     eta = boost_eta,
                     maxdepth = boost_maxdepth,
                     subsample = boost_subsample)

    # create the the xgboost DMatrix for both M and Y
    m_xgb_data <- xgboost::xgb.DMatrix(as.matrix(cbind('A' = a_treatment, x_covariates)),
                                       label = m_mediator)
    y_xgb_data <- xgboost::xgb.DMatrix(as.matrix(cbind('A' = a_treatment, 'M' = m_mediator, x_covariates)),
                                       label = y_outcome)

    # perform xgboost cross-validation for both M and Y
    m_cv <- xgboost::xgb.cv(params = m_params,
                            data = m_xgb_data,
                            nfold = boost_n_folds,
                            nrounds = boost_n_trees,
                            early_stopping_rounds = boost_early_stopping_rounds,
                            # print_every_n = 100,
                            verbose = verbose)
    
    y_cv <- xgboost::xgb.cv(params = y_params,
                            data = y_xgb_data,
                            nfold = boost_n_folds,
                            nrounds = boost_n_trees,
                            early_stopping_rounds = boost_early_stopping_rounds,
                            # print_every_n = 100,
                            verbose = verbose)

    # fit the final xgboost models, using the number if iterations from cross-validation
    model_m_res <- xgboost::xgboost(data = m_xgb_data,
                                    params = m_params,
                                    nrounds = m_cv$niter,
                                    verbose = FALSE)
    model_y_res <- xgboost::xgboost(data = y_xgb_data,
                                    params = y_params,
                                    nrounds = y_cv$niter,
                                    verbose = FALSE)

    # get the predictions for M1 and M0
    cf_A1 <- xgboost::xgb.DMatrix(as.matrix(cbind('A' = 1, x_covariates)))
    cf_A0 <- xgboost::xgb.DMatrix(as.matrix(cbind('A' = 0, x_covariates)))
    
    M1 <- predict(model_m_res, newdata = cf_A1, type = 'response')
    M0 <- predict(model_m_res, newdata = cf_A0, type = 'response')

  }

  # do simulations
  res <- matrix(NA, nrow = nsims, ncol = 4)
  for(i in 1:nsims){
    
    #if (dist_m == 'bernoulli') {
    m1 <- rbinom(n = n_obs, 1, prob = M1)
    m0 <- rbinom(n = n_obs, 1, prob = M0)
      
    #  # TODO : This will need to be changed
    #} else if (dist_m == 'gaussian') {
    #  m1 <- rnorm(n = n_obs, mean = M1)
    #  m0 <- rnorm(n = n_obs, mean = M0)
    # } else {
    #   stop('Something is wrong...')
    # }

    if (booster == "gbm") {
      a1m1 <- list(newdata = data.frame('X' = x_covariates,
                                        'A' = 1,
                                        'M' = m1),
                   n.trees = model_y_res$n.trees)
      a0m0 <- list(newdata = data.frame('X' = x_covariates,
                                        'A' = 0,
                                        'M' = m0),
                   n.trees = model_y_res$n.trees)
      a1m0 <- list(newdata = data.frame('X' = x_covariates,
                                        'A' = 1,
                                        'M' = m0),
                   n.trees = model_y_res$n.trees)
      a0m1 <- list(newdata = data.frame('X' = x_covariates,
                                        'A' = 0,
                                        'M' = m1),
                   n.trees = model_y_res$n.trees)
    } else {
      a1m1 <- list(xgboost::xgb.DMatrix(as.matrix(cbind('A' = 1, 'M' = m1, x_covariates))))
      a1m0 <- list(xgboost::xgb.DMatrix(as.matrix(cbind('A' = 1, 'M' = m0, x_covariates))))
      a0m1 <- list(xgboost::xgb.DMatrix(as.matrix(cbind('A' = 0, 'M' = m1, x_covariates))))
      a0m0 <- list(xgboost::xgb.DMatrix(as.matrix(cbind('A' = 0, 'M' = m0, x_covariates))))
    }
    
    treatment1_mediator1 <- do.call(predict, c(list(model_y_res), a1m1))
    treatment0_mediator0 <- do.call(predict, c(list(model_y_res), a0m0))
    treatment0_mediator1 <- do.call(predict, c(list(model_y_res), a0m1))
    treatment1_mediator0 <- do.call(predict, c(list(model_y_res), a1m0))
    
    res[i, ] <- c(mean(treatment1_mediator1),
                  mean(treatment0_mediator0),
                  mean(treatment0_mediator1),
                  mean(treatment1_mediator0))
  }
  
  # collect results
  res_means <- apply(res, 2, mean)

  treatment1_mediator1 <- res_means[1]
  treatment0_mediator0 <- res_means[2]
  treatment0_mediator1 <- res_means[3]
  treatment1_mediator0 <- res_means[4]

  results <- list('overall_effect' = treatment1_mediator1 - treatment0_mediator0,
                  'natural_direct_effect' = treatment1_mediator0 - treatment0_mediator0,
                  'natural_indirect_effect' = treatment1_mediator1 - treatment1_mediator0,
                  'treatment0_mediator0' = treatment0_mediator0,
                  'treatment1_mediator1' = treatment1_mediator1,
                  'treatment0_mediator1' = treatment0_mediator1,
                  'treatment1_mediator0' = treatment1_mediator0)
  class(results) <- "mediation"
  class(results) <- "outcome"
  return(results)
}