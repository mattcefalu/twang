cb.evaluation.log <- function(n.keep=1) {
   
   callback <- function(env = parent.frame(), finalize = FALSE) {
      
      if ( (env$iteration - env$begin_iteration + 1)%%n.keep == 0){
         if (is.data.table(env$evaluation_log)){
            env$evaluation_log[ , (paste0("pscore",env$iteration)) := env$bst_evaluation ]
         }else{
            env$evaluation_log = data.table( env$bst_evaluation)
            names(env$evaluation_log) = paste0("pscore",env$iteration)
         }
      }
      
      # if (!is.null(n.save)){
      #    if ( (env$iteration - env$begin_iteration +1) == n.save ){
      #       fwrite(transpose(env$evaluation_log) , file=file, sep = "," , row.names = F , col.names=F , quote=FALSE ) 
      #       env$evaluation_log[ , ( names(env$evaluation_log) ) := NULL ]
      #    }else{
      #       if ( (env$iteration - env$begin_iteration +1)%%n.save == 0 ){
      #          fwrite(transpose(env$evaluation_log) , file=file, sep = "," , row.names = F , col.names=F , quote=FALSE , append=T) 
      #          env$evaluation_log[ , ( names(env$evaluation_log) ) := NULL ]
      #       }
      #    }
      # }
      
      if (finalize)
         return()
   }
   attr(callback, 'call') <- match.call()
   attr(callback, 'name') <- 'cb.evaluation.log'
   callback
}
