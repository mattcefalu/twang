mnps<-function(formula,
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
             stop.method = c("es.max"), 
             sampw = NULL, version="fast",
             ks.exact=NULL,
             booster="xgboost",
             tree_method="hist",
             save.propensities=FALSE,
             file=NULL,
             n.keep = 1,
             n.grid = NULL,
             n.grid.ks = 25,
             n.grid.es = NULL,
             treatATT = NULL,
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
      if (!missing(n.grid.ks))    stop("Option n.grid.ks is not allowed with version='legacy'")
      if (!missing(n.grid.es))    stop("Option n.grid.es is not allowed with version='legacy'")
      
      return(mnps.old(formula = formula,
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
                    treatATT = treatATT, 
                    ...))
   }else{
      # throw error if user specifies params with booster=="gbm"
      if ( booster=="gbm" & !is.null(params) ) stop("params cannot be specified when booster='gbm'.")
      return(mnps.fast(formula = formula,
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
                     ks.exact=ks.exact,
                     booster=booster,
                     tree_method=tree_method,
                     save.propensities=save.propensities,
                     file=file,
                     n.keep = n.keep,
                     n.grid = n.grid,
                     n.grid.ks = n.grid.ks,
                     n.grid.es = n.grid.es,
                     treatATT = treatATT, 
                     ...))
   }
   
}
