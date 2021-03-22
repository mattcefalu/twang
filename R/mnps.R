#' Propensity score estimation for multiple treatments
#'
#' \code{mnps} calculates propensity scores for more than two treatment groups using gradient boosted
#' logistic regression, and diagnoses the resulting propensity scores using a variety of methods.
#' 
#' For user more comfortable with the options of \code{\link{xgboost}},
#' the options for \code{mnps} controlling the behavior of the gradient boosting
#' algorithm can be specified using the \code{\link{xgboost}} naming
#' scheme. This includes \code{nrounds}, \code{max_depth}, \code{eta}, and
#' \code{subsample}. In addition, the list of parameters passed to
#' \code{\link{xgboost}} can be specified with \code{params}.
#' 
#' Note that unlike earlier versions of \code{twang}, the plotting functions are
#' no longer included in the \code{mnps} function. See \code{\link[twang:plot.mnps]{plot}} for
#' details of the plots.
#' 
#' @param formula A formula for the propensity score model with the treatment
#'   indicator on the left side of the formula and the potential
#'   confounding variables on the right side.
#' @param data The dataset, includes treatment assignment as well as covariates.
#' @param n.trees Number of gbm iterations passed on to \code{\link{gbm}}. Default: 10000.
#' @param interaction.depth A positive integer denoting the tree depth used in
#'   gradient boosting. Default: 3.
#' @param shrinkage A numeric value between 0 and 1 denoting the learning rate.
#'   See \code{\link{gbm}} for more details. Default: 0.01.
#' @param bag.fraction A numeric value between 0 and 1 denoting the fraction of
#'   the observations randomly selected in each iteration of the gradient
#'   boosting algorithm to propose the next tree. See \code{\link{gbm}} for
#'   more details. Default: 1.0.
#' @param n.minobsinnode An integer specifying the minimum number of observations 
#'   in the terminal nodes of the trees used in the gradient boosting. See \code{\link{gbm}} for
#'   more details. Default: 10.
#' @param perm.test.iters A non-negative integer giving the number of iterations
#'   of the permutation test for the KS statistic. If \code{perm.test.iters=0}
#'   then the function returns an analytic approximation to the p-value. Setting
#'   \code{perm.test.iters=200} will yield precision to within 3\% if the true
#'   p-value is 0.05. Use \code{perm.test.iters=500} to be within 2\%. Default: 0.
#' @param print.level The amount of detail to print to the screen. Default: 2.
#' @param verbose If \code{TRUE}, lots of information will be printed to monitor the
#'   the progress of the fitting. Default: \code{TRUE}.
#' @param estimand \code{"ATE"} (average treatment effect) or \code{"ATT"} (average treatment
#'   effect on the treated) : the causal effect of interest. ATE estimates the
#'   change in the outcome if the treatment were applied to the entire
#'   population versus if the control were applied to the entire population. ATT
#'   estimates the analogous effect, averaging only over the treated population.
#'   Default: \code{"ATE"}.
#' @param stop.method A method or methods of measuring and summarizing balance across pretreatment
#'   variables. Current options are \code{ks.mean}, \code{ks.max}, \code{es.mean}, and \code{es.max}. \code{ks} refers to the
#'   Kolmogorov-Smirnov statistic and \code{es} refers to standardized effect size. These are summarized
#'   across the pretreatment variables by either the maximum (\code{.max}) or the mean (\code{.mean}). 
#'   Default: \code{c("es.mean")}.
#' @param sampw Optional sampling weights.
#' @param version \code{"gbm"}, \code{"xgboost"}, or \code{"legacy"}, indicating which version of the twang package to use.
#'   \itemize{
#'     \item{\code{"gbm"}}{ uses gradient boosting from the \code{\link{gbm}} package.}
#'     \item{\code{"xgboost"}}{ uses gradient boosting from the \code{\link{xgboost}} package.}
#'     \item{\code{"legacy"}}{ uses the prior implementation of the \code{\link{ps}} function.}
#'   }
#' @param ks.exact \code{NULL} or a logical indicating whether the
#'   Kolmogorov-Smirnov p-value should be based on an approximation of exact
#'   distribution from an unweighted two-sample Kolmogorov-Smirnov test. If
#'   \code{NULL}, the approximation based on the exact distribution is computed
#'   if the product of the effective sample sizes is less than 10,000.
#'   Otherwise, an approximation based on the asymptotic distribution is used.
#'   **Warning:** setting \code{ks.exact = TRUE} will add substantial
#'   computation time for larger sample sizes. Default: \code{NULL}.
#' @param n.keep A numeric variable indicating the algorithm should only
#'   consider every \code{n.keep}-th iteration of the propensity score model and
#'   optimize balance over this set instead of all iterations. Default: 1.
#' @param n.grid A numeric variable that sets the grid size for an initial
#'   search of the region most likely to minimize the \code{stop.method}. A
#'   value of \code{n.grid=50} uses a 50 point grid from \code{1:n.trees}. It
#'   finds the minimum, say at grid point 35. It then looks for the actual
#'   minimum between grid points 34 and 36. If specified with \code{n.keep>1}, \code{n.grid} 
#'   corresponds to a grid of points on the kept iterations as defined by \code{n.keep}. Default: 25.
#' @param treatATT If the estimand is specified to be ATT, this argument is 
#'   used to specify which treatment condition is considered 'the treated'.
#'   It must be one of the levels of the treatment variable. It is ignored for
#'   ATE analyses.
#' @param ... Additional arguments that are passed to \code{\link{ps}} function.
#'
#' @return Returns an object of class \code{mnps}, which consists of the following.  
#'   \itemize{
#'     \item{\code{psList}}{ A list of \code{\link{ps}} objects with length equal to the number of time periods.}
#'     \item{\code{nFits}}{ The number of \code{\link{ps}} objects (i.e., the number of distinct time points).}
#'     \item{\code{estimand}}{ The specified estimand.}
#'     \item{\code{treatATT}}{ For ATT fits, the treatment category that is considered "the treated".}
#'     \item{\code{treatLev}}{ The levels of the treatment variable.}
#'     \item{\code{levExceptTreatAtt}}{ The levels of the treatment variable, excluding the \code{treatATT} level.}
#'     \item{\code{data}}{ The data used to fit the model.}
#'     \item{\code{treatVar}}{ The vector of treatment indicators.}
#'     \item{\code{stopMethods}}{ The stopping rules specified in the call to \code{mnps}.}
#'     \item{\code{sampw}}{ Sampling weights provided to \code{mnps}, if any.}
#'   }
#'
#' @author Lane Burgette `<burgette@rand.org>`, Beth Ann Griffin `<bethg@rand.org>`,
#'   Dan Mc- Caffrey `<danielm@rand.org>`
#'
#' @references Dan McCaffrey, G. Ridgeway, Andrew Morral (2004). "Propensity
#'   Score Estimation with Boosted Regression for Evaluating Adolescent
#'   Substance Abuse Treatment", *Psychological Methods* 9(4):403-425.
#'
#' @seealso \code{\link{ps}}, \code{\link{gbm}}, \code{\link{xgboost}}, \code{\link[twang:plot.mnps]{plot}}, \code{\link{bal.table}}
#'
#' @export
mnps<-function(formula,
               data,
               # boosting options
               n.trees=10000,
               interaction.depth=3,
               shrinkage=0.01,
               bag.fraction = 1.0,
               n.minobsinnode =10,
               perm.test.iters=0,
               print.level=2,       
               verbose=TRUE,
               estimand="ATE", 
               stop.method = c("es.max"), 
               sampw = NULL, version="gbm",
               ks.exact=NULL,
               n.keep = 1,
               n.grid = 25,
               treatATT = NULL,
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
   
   # check that data is not a tibble or data.table
   if( class(data)[1] %in% c("tbl_df","tbl","data.table")){
      stop("The twang package currently does not support data.table or tibble. Please convert your data object to a data.frame.")
   }else{
      if (class(data)[1] != "data.frame"){
         warning("Data classes other than data.frame may cause errors." , call.=FALSE)
      }
   }
   
   if (version=="legacy"){
      ## throw some errors if the user specifies an option not allowed in legacy version of ps
      if (!is.null(ks.exact))     stop("Option ks.exact is not allowed with version='legacy'")
      if (!is.null(params))       stop("Option params is not allowed with version='legacy'")
      if (!is.null(args$tree_method))  stop("Option tree_method is not allowed with version='legacy'")
      if (!missing(n.keep))       stop("Option n.keep is not allowed with version='legacy'")
      if (!missing(n.grid))       stop("Option n.grid is not allowed with version='legacy'")
      if (!is.null(args$tree_method))  stop("Option tree_method is not allowed with version='legacy'")
      if (!missing(n.minobsinnode) | !is.null(args$min_child_weight)) stop("Options n.minobsinnode or min_child_weight are not allowed with version='legacy'")
      
      return(mnps.old(formula = formula,
                    data=data,                         # data
                    n.trees=nrounds,                   # gbm options
                    interaction.depth=max_depth,
                    shrinkage=eta,
                    bag.fraction = subsample,
                    perm.test.iters=perm.test.iters,
                    print.level=print.level,       
                    verbose=verbose,
                    estimand=estimand, 
                    stop.method = stop.method, 
                    sampw = sampw, 
                    treatATT = treatATT, 
                    ...))
   }else{
      # xgboost tree method  
      tree_method  <- if (!is.null(args$tree_method)) args$tree_method else "hist"
      
      # throw error if user specifies params with version=="gbm"
      if ( version=="gbm" & !is.null(params) ) stop("params cannot be specified when version='gbm'.")
      return(mnps.fast(formula = formula,
                     data=data,                         # data
                     n.trees=nrounds,                 # gbm options
                     interaction.depth=max_depth,
                     shrinkage=eta,
                     bag.fraction = subsample,
                     n.minobsinnode = min_child_weight,
                     params=params,
                     perm.test.iters=perm.test.iters,
                     print.level=print.level,       
                     verbose=verbose,
                     estimand=estimand, 
                     stop.method = stop.method, 
                     sampw = sampw, 
                     ks.exact=ks.exact,
                     version=version,
                     tree_method=tree_method,
                     n.keep = n.keep,
                     n.grid = n.grid,
                     treatATT = treatATT))
   }
   
}
