#' Inverse probability of treatment weighting for marginal structural models.
#'
#' `iptw` calculates propensity scores for sequential treatments using gradient boosted logistic
#' regression and diagnoses the resulting propensity scores using a variety of
#' methods
#' 
#' For user more comfortable with the options of [xgboost],
#' the options for `iptw` controlling the behavior of the gradient boosting
#' algorithm can be specified using the [xgboost] naming
#' scheme. This includes `nrounds`, `max_depth`, `eta`, and
#' `subsample`. In addition, the list of parameters passed to
#' [xgboost] can be specified with `params`.
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
#' @param interaction.depth A positive integer denoting the tree depth used in
#'   gradient boosting. Default: 3.
#' @param shrinkage A numeric value between 0 and 1 denoting the learning rate.
#'   See [gbm] for more details. Default: 0.01.
#' @param bag.fraction A numeric value between 0 and 1 denoting the fraction of
#'   the observations randomly selected in each iteration of the gradient
#'   boosting algorithm to propose the next tree. See [gbm] for
#'   more details. Default: 1.0.
#'  @param n.minobsinnode An integer specifying the minimum number of observations 
#'   in the terminal nodes of the trees used in the gradient boosting.  See [gbm] for
#'   more details. Default: 10.
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
#' @param version "gbm", "xgboost", or "legacy", indicating which version of the twang package to use.
#'   * `"gbm"` Uses gradient boosting from the `gbm` package.
#'   * `"xgboost"` Uses gradient boosting from the `xgboost` package.
#'   * `"legacy"` Uses the prior implementation of the `ps`` function.
#' @param ks.exact `NULL` or a logical indicating whether the
#'   Kolmogorov-Smirnov p-value should be based on an approximation of exact
#'   distribution from an unweighted two-sample Kolmogorov-Smirnov test. If
#'   `NULL`, the approximation based on the exact distribution is computed
#'   if the product of the effective sample sizes is less than 10,000.
#'   Otherwise, an approximation based on the asymptotic distribution is used.
#'   **Warning:** setting `ks.exact = TRUE` will add substantial
#'   computation time for larger sample sizes. Default: `NULL`.
#' @param n.keep A numeric variable indicating the algorithm should only
#'   consider every `n.keep`-th iteration of the propensity score model and
#'   optimize balance over this set instead of all iterations. Default: 1.
#' @param n.grid A numeric variable that sets the grid size for an initial
#'   search of the region most likely to minimize the `stop.method`. A
#'   value of `n.grid=50` uses a 50 point grid from `1:n.trees`. It
#'   finds the minimum, say at grid point 35. It then looks for the actual
#'   minimum between grid points 34 and 36. If specified with `n.keep>1`, `n.grid` 
#'   corresponds to a grid of points on the kept iterations as defined by ```n.keep```. Default: 25.
#' @param ... Additional arguments that are passed to ps function.
#'
#' @return Returns an object of class `iptw`, a list containing
#'   * `psList` A list of `ps` objects with length equal to the number of time periods.
#'   * `estimand` The specified estimand.
#'   * `stop.methods` The stopping rules used to optimize `iptw` balance.
#'   * `nFits` The number of `ps` objects (i.e., the number of distinct time points).
#'   * `uniqueTimes` The unique times in the specified model.
#'
#' @seealso [ps], [mnps], [gbm], [xgboost], [plot.iptw], [bal.table]
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
                 n.trees = 10000,
                 interaction.depth = 3,
                 shrinkage = 0.01,
                 bag.fraction = 1.0,
                 n.minobsinnode =10,
                 perm.test.iters = 0,
                 print.level = 2,       
                 verbose = TRUE,
                 stop.method = c("es.max"), 
                 sampw = NULL, 
                 version = "gbm",
                 ks.exact = NULL,
                 n.keep = 1,
                 n.grid = 25,
                 ...){

   # collect named arguments from dots
   args         <- list(...)
   args_named   <- names(args)

   # parse hidden options and xgboost parameters
   params       <- args$params
   max_depth    <- if (!is.null(args$max_depth)) args$max_depth else interaction.depth
   subsample    <- if (!is.null(args$subsample)) args$subsample else bag.fraction
   nrounds      <- if (!is.null(args$nrounds)) args$nrounds else n.trees
   eta          <- if (!is.null(args$eta)) args$eta else shrinkage
   min_child_weight <- if (!is.null(args$min_child_weight)) args$min_child_weight else n.minobsinnode
   

   # throw some errors if the user specifies two versions of the same option
   if (!missing(n.trees) & ('nrounds' %in% args_named))             stop("Only one of n.trees and nrounds can be specified.")
   if (!missing(interaction.depth) & ('max_depth' %in% args_named)) stop("Only one of interaction.depth and max_depth can be specified.")
   if (!missing(shrinkage) & ('eta' %in% args_named))               stop("Only one of shrinkage and eta can be specified.")
   if (!missing(bag.fraction) & ('subsample' %in% args_named))      stop("Only one of shrinkage and eta can be specified.")
   if (!missing(n.minobsinnode) & ('min_child_weight' %in% args_named))      stop("Only one of n.minobsinnode and min_child_weight can be specified.")
   
   # throw error if user specifies params with other options
   if (!missing(interaction.depth) | ('max_depth' %in% args_named) | !missing(shrinkage) | ('eta' %in% args_named) | !missing(bag.fraction) | ('subsample' %in% args_named) ){
      if (!is.null(params)) stop("params cannot be specified with any of interaction.depth, max_depth, shrinkage, eta, bag.fraction, or subsample.")
   }

   if (version=="legacy"){
      ## throw some errors if the user specifies an option not allowed in legacy version of ps
      if (!is.null(ks.exact))     stop("Option ks.exact is not allowed with version='legacy'")
      if (!is.null(params))       stop("Option params is not allowed with version='legacy'")
      if (!missing(tree_method))  stop("Option tree_method is not allowed with version='legacy'")
      if (!missing(n.keep))       stop("Option n.keep is not allowed with version='legacy'")
      if (!missing(n.grid))       stop("Option n.grid is not allowed with version='legacy'")
      if (!is.null(args$tree_method))  stop("Option tree_method is not allowed with version='legacy'")
      if (!missing(n.minobsinnode) | !is.null(args$min_child_weight)) stop("Options n.minobsinnode or min_child_weight are not allowed with version='legacy'")
      
      return(iptw.old(formula = formula,
                      data = data,
                      # iptw options
                      timeInvariant = timeInvariant,
                      cumulative = cumulative, 
                      timeIndicators = timeIndicators, 
                      ID = ID, 
                      priorTreatment = priorTreatment,
                      # gbm options
                      n.trees = nrounds,                 
                      interaction.depth = max_depth,
                      shrinkage = eta,
                      bag.fraction = subsample,
                      perm.test.iters = perm.test.iters,
                      print.level = print.level,       
                      verbose = verbose,
                      stop.method = stop.method, 
                      sampw = sampw, 
                      ...))
   }else{
      # xgboost tree method  
      tree_method  <- if (!is.null(args$tree_method)) args$tree_method else "hist"
      
      # throw error if user specifies params with version=="gbm"
      if ( version=="gbm" & !is.null(params) ) stop("params cannot be specified when version='gbm'.")
      return(iptw.fast(formula = formula,
                       data = data,
                       # iptw options
                       timeInvariant = timeInvariant,
                       cumulative = cumulative, 
                       timeIndicators = timeIndicators, 
                       ID = ID, 
                       priorTreatment = priorTreatment,
                       # boosting options
                       n.trees = nrounds,                 
                       interaction.depth = max_depth,
                       shrinkage = eta,
                       bag.fraction = subsample,
                       n.minobsinnode = min_child_weight,
                       params = params,
                       perm.test.iters = perm.test.iters,
                       print.level = print.level,       
                       verbose = verbose,
                       stop.method = stop.method, 
                       sampw = sampw, 
                       ks.exact = ks.exact,
                       version = version,
                       tree_method = tree_method,
                       n.keep = n.keep,
                       n.grid = n.grid))
   }
	
}