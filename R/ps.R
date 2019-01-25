ps<-function(formula = formula(data),
             data,                         # data
             n.trees=10000,                 # gbm options
             interaction.depth=3,
             shrinkage=0.01,
             bag.fraction = 1.0,
             perm.test.iters=0,
             print.level=2,       
             verbose=TRUE,
             estimand="ATE", 
             stop.method = c("ks.mean", "es.mean"), 
             sampw = NULL, multinom = FALSE, version="fast",
             ks.exact=NULL,
             booster="xgboost",
             tree_method="hist",
             save.propensities=FALSE,
             file=NULL,
             n.keep = 1,
             n.grid = NULL,
             n.grid.ks = 25,
             n.grid.es = NULL,
             ...){
  
  if (version=="legacy"){
    if ( !is.null(ks.exact) ){
      warning("Option ks.exact is not allowed with version='legacy'")
    }
    return(ps.old(formula = formula,
                  data=data,                         # data
                  n.trees=n.trees,                 # gbm options
                  interaction.depth=interaction.depth,
                  shrinkage=shrinkage,
                  bag.fraction = bag.fraction,
                  perm.test.iters=perm.test.iters,
                  print.level=print.level,       
                  verbose=verbose,
                  estimand=estimand, 
                  stop.method = stop.method, 
                  sampw = sampw, 
                  multinom = multinom, 
                  ...))
  }else{
    return(ps.fast(formula = formula,
                  data=data,                         # data
                  n.trees=n.trees,                 # gbm options
                  interaction.depth=interaction.depth,
                  shrinkage=shrinkage,
                  bag.fraction = bag.fraction,
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
                  n.grid.ks = n.grid.ks,
                  n.grid.es = n.grid.es,
                  ...))
  }
  
}
  