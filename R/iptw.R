#' Inverse probability of treatment weighting for marginal structural models.
#'
#' `iptw` uses [gbm] to estimate propensity scores for sequential treatments.
#'
#' This function uses generalized boosted models to estimate inverse probability of
#' treatment weights for sequential treatments.
#'
#' @param formula Either a single formula (long format) or a list with formulas.
#' @param data The dataset, includes treatment assignment as well as covariates.
#' @param timeInvariant An optional formula (with no left-hand variable) specifying time-invariant
#'   chararacteristics.
#' @param cumulative If `TRUE`, the time t model includes time-varying characteristics from times 1
#'   through t-1. Default: `TRUE`.
#' @param timeIndicators For long format fits, a vector of times for each observation.
#' @param ID For long format fits, a vector of numeric identifiers for unique analytic units.
#' @param priorTreatment For long format fits, includes treatment levels from previous times if `TRUE`. This
#'   argument is ignored for wide format fits. Default: `TRUE`.
#' @param n.trees Number of gbm iterations passed on to [gbm].
#' @param nrounds Equivalent to `n.trees`.
#' @param interaction.depth A positive integer denoting the tree depth used in
#'   gradient boosting. Only one of `interaction.depth` and
#'   `max_depth` can be specified. Default: 3.
#' @param max_depth Equivalent to `interaction.depth`.
#' @param shrinkage A numeric value between 0 and 1 denoting the learning rate.
#'   See [gbm] for more details. Only one of `shrinkage`
#'   and `eta` can be specified. Default: 0.01.
#' @param eta Equivalent to `shrinkage`.
#' @param bag.fraction A numeric value between 0 and 1 denoting the fraction of
#'   the observations randomly selected in each iteration of the gradient
#'   boosting algorithm to propose the next tree. See [gbm] for
#'   more details. Only one of `bag.fraction` and `subsample` can be
#'   specified. Default: 1.0.
#' @param subsample Equivalent to `bag.fraction`.
#' @param params ...
#' @param perm.test.iters A non-negative integer giving the number of iterations
#'   of the permutation test for the KS statistic. If `perm.test.iters=0`
#'   then the function returns an analytic approximation to the p-value. Setting
#'   `perm.test.iters=200` will yield precision to within 3\% if the true
#'   p-value is 0.05. Use `perm.test.iters=500` to be within 2\%. Default: 0.
#' @param print.level The amount of detail to print to the screen. Default: 2.
#' @param verbose If `TRUE`, lots of information will be printed to monitor the
#'   the progress of the fitting. Default: `TRUE`.
#' @param stop.method A method or methods of measuring and summarizing balance across pretreatment
#'   variables. Current options are `ks.mean`, `ks.max`, `es.mean`, and `es.max`. `ks` refers to the
#'   Kolmogorov-Smirnov statistic and es refers to standardized effect size. These are summarized
#'   across the pretreatment variables by either the maximum (`.max`) or the mean (`.mean`). 
#'   Default: `c("es.max")`.
#' @param sampw Optional sampling weights.
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
#' @param save.propensities Whether to save the propensity scores. Default: `FALSE`.
#' @param file Optional RDS file path where the `ps` object will be saved.
#' @param n.keep A numeric variable indicating the algorithm should only
#'   consider every `n.keep`-th iteration of the propensity score model and
#'   optimize balance over this set instead of all iterations. Default: 1.
#' @param n.grid A numeric variable that sets the grid size for an initial
#'   search of the region most likely to minimize the `stop.method`. A
#'   value of `n.grid=50` uses a 50 point grid from `1:n.trees`. It
#'   finds the minimum, say at grid point 35. It then looks for the actual
#'   minimum between grid points 34 and 36. Default: 25.
#' @param ... Additional arguments that are passed to ps function.
#'
#' @return Returns an object of class `iptw`, a list containing
#'   * `psList` A list of `ps` objects with length equal to the number of time periods.
#'   * `estimand` The specified estimand.
#'   * `stop.methods` The stopping rules used to optimize `iptw` balance.
#'   * `nFits` The number of `ps` objects (i.e., the number of distinct time points).
#'   * `uniqueTimes` The unique times in the specified model.
#'
#' @seealso [ps]
#'
#' @export
iptw <- function(formula, 
                 data,  
                 # iptw options
                 timeInvariant = NULL,
                 cumulative = TRUE, 
                 timeIndicators = NULL, 
                 ID = NULL, 
                 priorTreatment = TRUE,
                 # boosting options -- gbm version first, then xgboost name
                 n.trees = 10000, nrounds = n.trees,
                 interaction.depth = 3, max_depth = interaction.depth,
                 shrinkage = 0.01, eta = shrinkage , 
                 bag.fraction = 1.0, subsample = bag.fraction,
                 params = NULL,
                 perm.test.iters = 0,
                 print.level = 2,       
                 verbose = TRUE,
                 stop.method = c("es.max"), 
                 sampw = NULL, 
                 version = "fast",
                 ks.exact = NULL,
                 booster = "gbm",
                 tree_method = "hist",
                 save.propensities = FALSE,
                 file = NULL,
                 n.keep = 1,
                 n.grid = 25,
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
      
      return(iptw.old(formula = formula,
                      data=data,                         # data
                      # iptw options
                      timeInvariant = timeInvariant,
                      cumulative = cumulative, 
                      timeIndicators = timeIndicators, 
                      ID = ID, 
                      priorTreatment = priorTreatment,
                      # gbm options
                      n.trees=nrounds,                 
                      interaction.depth=max_depth,
                      shrinkage=eta,
                      bag.fraction = subsample,
                      perm.test.iters=perm.test.iters,
                      print.level=print.level,       
                      verbose=verbose,
                      stop.method = stop.method, 
                      sampw = sampw, 
                      ...))
   }else{
      # throw error if user specifies params with booster=="gbm"
      if ( booster=="gbm" & !is.null(params) ) stop("params cannot be specified when booster='gbm'.")
      return(iptw.fast(formula = formula,
                       data=data,                         # data
                       # iptw options
                       timeInvariant = timeInvariant,
                       cumulative = cumulative, 
                       timeIndicators = timeIndicators, 
                       ID = ID, 
                       priorTreatment = priorTreatment,
                       # boosting options
                       n.trees=nrounds,                 
                       interaction.depth=max_depth,
                       shrinkage=eta,
                       bag.fraction = subsample,
                       params=params,
                       perm.test.iters=perm.test.iters,
                       print.level=print.level,       
                       verbose=verbose,
                       stop.method = stop.method, 
                       sampw = sampw, 
                       ks.exact=ks.exact,
                       booster=booster,
                       tree_method=tree_method,
                       save.propensities=save.propensities,
                       file=file,
                       n.keep = n.keep,
                       n.grid = n.grid,
                       ...))
   }
	
}