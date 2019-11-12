cb.print.evaluation <- function (period = 10) {
   callback <- function(env = parent.frame()) {
      if (is.null(env$time)){
         env$time = Sys.time()
      }
      i <- env$iteration 
      if (i == env$begin_iteration){
         cat('Iteration', i  ,"\n")
      }
      if ( (i>env$begin_iteration & i<(10+env$begin_iteration)) || (i%%10 == 0 & i<(100+env$begin_iteration)) || (i%%100 == 0)  || i == env$end_iteration) {
         cat('Iteration', i , ": Total boosting time =" , round(difftime(Sys.time(),env$time  , units="mins"),2)   ,"minutes\n")
      }
   }
   attr(callback, 'call') <- match.call()
   attr(callback, 'name') <- 'cb.print.evaluation'
   callback
}
