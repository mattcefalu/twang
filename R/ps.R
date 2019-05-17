#' Gradient boosted propensity score estimation
#'
#' `ps` calculates propensity scores using gradient boosted logistic
#' regression and diagnoses the resulting propensity scores using a variety of
#' methods
#' 
#' For user more comfortable with the options of [xgboost],
#' the options for `ps` controlling the behavior of the gradient boosting
#' algorithm can be specified using the [xgboost] naming
#' scheme. This includes `nrounds`, `max_depth`, `eta`, and
#' `subsample`. In addition, the list of parameters passed to
#' [xgboost] can be specified with `params`.
#' 
#' Note that unlike earlier versions of `twang`, the plotting functions are
#' no longer included in the `ps` function. See  [plot] for
#' details of the plots.
#' 
#' @param formula An object of class [formula]: a symbolic
#'   description of the propensity score model to be fit with the treatment
#'   indicator on the left side of the formula and the potential confounding
#'   variables on the right side.
#' @param data A dataset that includes the treatment indicator as well as the
#'   potential confounding variables.
#' @param n.trees Number of gbm iterations passed on to [gbm].
#' @param interaction.depth A positive integer denoting the tree depth used in
#'   gradient boosting. Only one of `interaction.depth` and
#'   `max_depth` can be specified. Default: 3.
#' @param shrinkage A numeric value between 0 and 1 denoting the learning rate.
#'   See [gbm] for more details. Only one of `shrinkage`
#'   and `eta` can be specified. Default: 0.01.
#' @param bag.fraction A numeric value between 0 and 1 denoting the fraction of
#'   the observations randomly selected in each iteration of the gradient
#'   boosting algorithm to propose the next tree. See [gbm] for
#'   more details. Only one of `bag.fraction` and `subsample` can be
#'   specified. Default: 1.0.
#' @param params ...
#' @param perm.test.iters A non-negative integer giving the number of iterations
#'   of the permutation test for the KS statistic. If `perm.test.iters=0`
#'   then the function returns an analytic approximation to the p-value. Setting
#'   `perm.test.iters=200` will yield precision to within 3\% if the true
#'   p-value is 0.05. Use `perm.test.iters=500` to be within 2\%. Default: 0.
#' @param print.level The amount of detail to print to the screen. Default: 2.
#' @param verbose If `TRUE`, lots of information will be printed to monitor the
#'   the progress of the fitting. Default: `TRUE`.
#' @param estimand "ATE" (average treatment effect) or "ATT" (average treatment
#'   effect on the treated) : the causal effect of interest. ATE estimates the
#'   change in the outcome if the treatment were applied to the entire
#'   population versus if the control were applied to the entire population. ATT
#'   estimates the analogous effect, averaging only over the treated population.
#'   Default: "ATE".
#' @param stop.method A method or methods of measuring and summarizing balance across pretreatment
#'   variables. Current options are `ks.mean`, `ks.max`, `es.mean`, and `es.max`. `ks` refers to the
#'   Kolmogorov-Smirnov statistic and es refers to standardized effect size. These are summarized
#'   across the pretreatment variables by either the maximum (`.max`) or the mean (`.mean`). 
#'   Default: `c("ks.mean", "es.mean")`.
#' @param sampw Optional sampling weights.
#' @param multinom `TRUE` if used for multinomial propensity scores. Default: `FALSE`.
#' @param version Legacy?
#' @param ks.exact `NULL` or a logical indicating whether the
#'   Kolmogorov-Smirnov p-value should be based on an approximation of exact
#'   distribution from an unweighted two-sample Kolmogorov-Smirnov test. If
#'   `NULL`, the approximation based on the exact distribution is computed
#'   if the product of the effective sample sizes is less than 10,000.
#'   Otherwise, an approximation based on the asymptotic distribution is used.
#'   **Warning:** setting `ks.exact = TRUE` will add substantial
#'   computation time for larger sample sizes. Default: `NULL`.
#' @param booster `gbm` or `xgboost`. Default: `gbm`.
#' @param tree_method `xgboost` param. Default: "hist".
#' @param n.keep A numeric variable indicating the algorithm should only
#'   consider every `n.keep`-th iteration of the propensity score model and
#'   optimize balance over this set instead of all iterations. Default: 1.
#' @param n.grid A numeric variable that sets the grid size for an initial
#'   search of the region most likely to minimize the `stop.method`. A
#'   value of `n.grid=50` uses a 50 point grid from `1:n.trees`. It
#'   finds the minimum, say at grid point 35. It then looks for the actual
#'   minimum between grid points 34 and 36. Default: 25.
#'
#' @return Returns an object of class `ps`, a list containing 
#'   * `gbm.obj` The returned [gbm] object.
#'   * `treat` The treatment variable.
#'   * `desc` A list containing balance tables for each method selected in
#'     `stop.methods`. Includes a component for the unweighted
#'     analysis names \dQuote{unw}. Each `desc` component includes
#'     a list with the following components
#'     - `ess` The effective sample size of the control group.
#'     - `n.treat` The number of subjects in the treatment group.
#'     - `n.ctrl` The number of subjects in the control group.
#'     - `max.es` The largest effect size across the covariates.
#'     - `mean.es` The mean absolute effect size.
#'     - `max.ks` The largest KS statistic across the covariates.
#'     - `mean.ks` The average KS statistic across the covariates.
#'     - `bal.tab` a (potentially large) table summarizing the quality of the 
#'       weights for equalizing the distribution of features across 
#'       the two groups. This table is best extracted using the
#'       [bal.table] method. See the help for [bal.table] for details
#'       on the table's contents.
#'     - `n.trees` The estimated optimal number of [gbm]
#'       iterations to optimize the loss function for the associated 
#'        `stop.methods`.
#'     - `ps` a data frame containing the estimated propensity scores. Each
#'       column is associated with one of the methods selected in `stop.methods`.
#'     - `w` a data frame containing the propensity score weights. Each
#'       column is associated with one of the methods selected in `stop.methods`.
#'       If sampling weights are given then these are incorporated into these weights.
#'     - `estimand` The estimand of interest (ATT or ATE).
#'  * `datestamp` Records the date of the analysis.
#'  * `parameters` Saves the `ps` call.
#'  * `alerts` Text containing any warnings accumulated during the estimation.
#'  * `iters` A sequence of iterations used in the GBM fits used by `plot` function.
#'  * `balance` The balance measures for the pretreatment covariates, with a column for each
#'    `stop.method`.
#'  * `n.trees` Maximum number of trees considered in GBM fit.
#'  * `data` Data as specified in the `data` argument.
#'
#' @seealso [gbm]
#' @keywords models, multivariate
#'
#' @references Dan McCaffrey, G. Ridgeway, Andrew Morral (2004). "Propensity
#'   Score Estimation with Boosted Regression for Evaluating Adolescent
#'   Substance Abuse Treatment", *Psychological Methods* 9(4):403-425.
#'
#' @export
ps<-function(formula = formula(data),
             data,                         # data
             # boosting options -- gbm version first, then xgboost name
             n.trees = 10000, nrounds = n.trees,
             interaction.depth = 3, max_depth = interaction.depth,
             shrinkage = 0.01, eta = shrinkage , 
             bag.fraction = 1.0, subsample = bag.fraction,
             params = NULL,
             perm.test.iters = 0,
             print.level = 2,       
             verbose = TRUE,
             estimand = "ATE", 
             stop.method = c("ks.mean", "es.mean"), 
             sampw = NULL, 
             multinom = FALSE, 
             version = "fast",
             ks.exact = NULL,
             booster = "gbm",
             tree_method = "hist",
             #save.propensities = FALSE,
             #file = NULL,
             n.keep = 1,
             n.grid = 25,
             #n.grid.ks = 25,
             #n.grid.es = NULL,
             ...){
   
   ## throw some errors if the user specifies two versions of the same option
   if (!missing(n.trees) & !missing(nrounds)) stop("Only one of n.trees and nrounds can be specified.")
   if (!missing(interaction.depth) & !missing(max_depth)) stop("Only one of interaction.depth and max_depth can be specified.")
   if (!missing(shrinkage) & !missing(eta)) stop("Only one of shrinkage and eta can be specified.")
   if (!missing(bag.fraction) & !missing(subsample)) stop("Only one of shrinkage and eta can be specified.")
   
   # throw error if user specifies params with other options
   if (!missing(interaction.depth) | !missing(max_depth) | !missing(shrinkage) | !missing(eta) | !missing(bag.fraction) | !missing(subsample) ){
      if (!is.null(params)) stop("params cannot be specified with any of interaction.depth, max_depth, shrinkage, eta, bag.fraction, or subsample.")
   }
   
   if (version=="legacy"){
      ## throw some errors if the user specifies an option not allowed in legacy version of ps
      if (!is.null(ks.exact))     stop("Option ks.exact is not allowed with version='legacy'")
      if (!is.null(params))     stop("Option params is not allowed with version='legacy'")
      if (!missing(booster))      stop("Option booster is not allowed with version='legacy'")
      if (!missing(tree_method))  stop("Option tree_method is not allowed with version='legacy'")
      if (!missing(n.keep))       stop("Option n.keep is not allowed with version='legacy'")
      if (!missing(n.grid))       stop("Option n.grid is not allowed with version='legacy'")
      #if (!missing(n.grid.ks))    stop("Option n.grid.ks is not allowed with version='legacy'")
      #if (!missing(n.grid.es))    stop("Option n.grid.es is not allowed with version='legacy'")
      
      return(ps.old(formula = formula,
                     data=data,                         # data
                     n.trees=nrounds,                 # gbm options
                     interaction.depth=max_depth,
                     shrinkage=eta,
                     bag.fraction = subsample,
                     perm.test.iters=perm.test.iters,
                     print.level=print.level,       
                     verbose=verbose,
                     estimand=estimand, 
                     stop.method = stop.method, 
                     sampw = sampw, 
                     multinom = multinom, 
                     ...))
  }else{
    # throw error if user specifies params with booster=="gbm"
    if ( booster=="gbm" & !is.null(params) ) stop("params cannot be specified when booster='gbm'.")
    return(ps.fast(formula = formula,
                  data=data,                         # data
                  n.trees=nrounds,                 # gbm options
                  interaction.depth=max_depth,
                  shrinkage=eta,
                  bag.fraction = subsample,
                  params=params,
                  perm.test.iters=perm.test.iters,
                  print.level=print.level,       
                  verbose=verbose,
                  estimand=estimand, 
                  stop.method = stop.method, 
                  sampw = sampw, 
                  multinom = multinom,
                  ks.exact=ks.exact,
                  booster=booster,
                  tree_method=tree_method,
                  save.propensities=save.propensities,
                  file=file,
                  n.keep = n.keep,
                  n.grid = n.grid,
                 # n.grid.ks = n.grid.ks,
                  #n.grid.es = n.grid.es,
                  ...))
  }
  
}
  