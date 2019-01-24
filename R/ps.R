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
                  ...))
  }
  
}
  