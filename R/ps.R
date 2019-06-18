#' Gradient boosted propensity score estimation
#'
#' \code{ps} calculates propensity scores using gradient boosted logistic
#' regression and diagnoses the resulting propensity scores using a variety of
#' methods
#' 
#' For user more comfortable with the options of \code{\link[xgboost]{xgboost}},
#' the options for \code{ps} controlling the behavior of the gradient boosting
#' algorithm can be specified using the \code{\link[xgboost]{xgboost}} naming
#' scheme. This includes \code{nrounds}, \code{max_depth}, \code{eta}, and
#' \code{subsample}. In addition, the list of parameters passed to
#' \code{\link[xgboost]{xgboost}} can be specified with \code{params}.
#' 
#' Note that unlike earlier versions of \code{twang}, the plotting functions are
#' no longer included in the \code{ps} function. See  \code{\link{plot}} for
#' details of the plots.
#' 
#' @param formula an object of class "\code{\link{formula}}": a symbolic
#'   description of the propensity score model to be fit with the treatment
#'   indicator on the left side of the formula and the potential confounding
#'   variables on the right side.
#' @param data a dataset that includes the treatment indicator as well as the
#'   potential confounding variables.
#' @param estimand "ATE" (average treatment effect) or "ATT" (average treatment
#'   effect on the treated) : the causal effect of interest. ATE estimates the
#'   change in the outcome if the treatment were applied to the entire
#'   population versus if the control were applied to the entire population. ATT
#'   estimates the analogous effect, averaging only over the treated population.
#'   Default: "ATE"
#' @param stop.method a character vector of the method or methods to be used to
#'   optimize the balance of pretreatment variables. Current options are
#'   \code{ks.mean}, \code{ks.max}, \code{es.mean}, and \code{es.max}. \code{ks}
#'   refers to the Kolmogorov-Smirnov statistic and \code{es} refers to
#'   standardized effect size. These are summarized across the pretreatment
#'   variables by either the maximum (\code{.max}) or the mean (\code{.mean}).
#' @param n.trees a positive integer denoting the maximum number of iterations
#'   to be used in the gradient boosting algorithm. Only one of \code{n.trees}
#'   and \code{nrounds} can be specified. Default: 10000
#' @param interaction.depth a positive integer denoting the tree depth used in
#'   gradient boosting. Only one of \code{interaction.depth} and
#'   \code{max_depth} can be specified. Default: 3
#' @param shrinkage a numeric value between 0 and 1 denoting the learning rate.
#'   See \code{\link[gbm]{gbm}} for more details. Only one of \code{shrinkage}
#'   and \code{eta} can be specified. Default: 0.01
#' @param bag.fraction a numeric value between 0 and 1 denoting the fraction of
#'   the observations randomly selected in each iteration of the gradient
#'   boosting algorithm to propose the next tree. See \code{\link[gbm]{gbm}} for
#'   more details. Only one of \code{bag.fraction} and \code{subsample} can be
#'   specified. Default: 1
#' @param perm.test.iters a non-negative integer giving the number of iterations
#'   of the permutation test for the KS statistic. If \code{perm.test.iters=0}
#'   then the function returns an analytic approximation to the p-value. Setting
#'   \code{perm.test.iters=200} will yield precision to within 3\% if the true
#'   p-value is 0.05. Use \code{perm.test.iters=500} to be within 2\%. Default:
#'   0
#' @param print.level the amount of detail to print to the screen
#' @param verbose if TRUE, lots of information will be printed to monitor the
#'   the progress of the fitting
#' @param sampw Optional sampling weights.
#' @param version Legacy?
#' @param ks.exact \code{NULL} or a logical indicating whether the
#'   Kolmogorov-Smirnov p-value should be based on an approximation of exact
#'   distribution from an unweighted two-sample Kolmogorov-Smirnov test. If
#'   \code{NULL}, the approximation based on the exact distribution is computed
#'   if the product of the effective sample sizes is less than 10,000.
#'   Otherwise, an approximation based on the asymptotic distribution is used.
#'   \strong{Warning:} setting \code{ks.exact = TRUE} will add substantial
#'   computation time for larger sample sizes.
#' @param booster gbm or xgboost
#' @param tree_method xgboost param
#' @param n.keep a numeric variable indicating the algorithm should only
#'   consider every \code{n.keep}-th iteration of the propensity score model and
#'   optimize balance over this set instead of all iterations.
#' @param n.grid a numeric variable that sets the grid size for an initial
#'   search of the region most likely to minimize the \code{stop.method}. A
#'   value of \code{n.grid=50} uses a 50 point grid from \code{1:n.trees}. It
#'   finds the minimum, say at grid point 35. It then looks for the actual
#'   minimum between grid points 34 and 36.
#' @param moderator a string containing the name of the moderator variable, which
#'   must correspond to a factor variable in \code{data} that appears on the right-hand
#'   side of \code{formula}. If non-\code{NULL}, balance will be assessed not only
#'   overall, but also within subgroups of the moderator.
#'   
#' @references Dan McCaffrey, G. Ridgeway, Andrew Morral (2004). "Propensity
#'   Score Estimation with Boosted Regression for Evaluating Adolescent
#'   Substance Abuse Treatment", \emph{Psychological Methods} 9(4):403-425.

ps<-function(formula = formula(data),
             data,                         # data
             # boosting options -- gbm version first, then xgboost name
             n.trees=10000, nrounds=n.trees,
             interaction.depth=3, max_depth=interaction.depth,
             shrinkage=0.01, eta=shrinkage , 
             bag.fraction = 1.0, subsample=bag.fraction,
             params=NULL,
             perm.test.iters=0,
             print.level=2,       
             verbose=TRUE,
             estimand="ATE", 
             stop.method = c("ks.mean", "es.mean"), 
             sampw = NULL, 
             multinom = FALSE, 
             version="fast",
             ks.exact=NULL,
             booster="gbm",
             tree_method="hist",
             #save.propensities=FALSE,
             #file=NULL,
             n.keep = 1,
             n.grid = 25,
             #n.grid.ks = 25,
             #n.grid.es = NULL,
             moderator = NULL,
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
      if (!is.null(moderator))     stop("Option moderator is not allowed with version='legacy'")
      
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
                 moderator = moderator,
                  ...))
  }
  
}
  