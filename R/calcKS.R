calcKS <- function(data, w , treat , multinomATE=FALSE, sw=NULL, moderator.var = NULL){
  
  if (!is.null(moderator.var)) {
    moderator = data[,moderator.var]
    groups <- c("all", levels(moderator))
    data <- data[,!colnames(data) %in% moderator.var, drop = FALSE]
  }
  else {
    groups <- ""
  }
  
  W1 = sweep(w , 1 , treat , "*") #apply(w,2,function(x) x*treat )
  if (multinomATE){
    W0 = as.matrix(sw)
  }else{
    W0 = sweep(w , 1 , 1-treat , "*") #apply(w,2,function(x) x*(1-treat))
  }
  
  ks.effect.list <- lapply(groups, function(g) {
    
    if (g %in% c("", "all")) {
      in. <- rep(TRUE, nrow(data))
      data_g <- data
      if (g == "") g_name <- ""
      else g_name <- paste0(g, "_")
    }
    else {
      in. <- moderator == g
      #Drop moderator variable from subgroup balance
      data_g <- data[in.,!colnames(data) %in% paste(moderator.var, levels(moderator), sep = ":"), drop = FALSE]
      g_name <- paste0(g, "_")
    }
    
    W1_g <- W1[in.,,drop = FALSE]
    W0_g <- W0[in.,,drop = FALSE]
    sw_g <- sw[in.]
    
    ks.effect = NULL
    for (j in 1:ncol(data_g)){
      
      # order of the data, dropping NA. Only numeric variables have NAs.
      index = order(data_g[,j], na.last = NA)
      
      # consider last occurance of unique values
      dups = !duplicated(data_g[index,j], fromLast=T)
      
      if (multinomATE){
        # handles the fact that W0 is a vector when multinomATE is true
        dW = abs(sweep(apply(W1_g[index,,drop=FALSE] , 2 , function(y) cumsum(y)/sum(y) ) , 1 ,  cumsum(sw_g[index])/sum(sw_g[index]) , "-" ) )
      }else{
        dW = abs(apply(W1_g[index,,drop=FALSE] , 2 , function(y) cumsum(y)/sum(y) ) - apply(W0_g[index,,drop=FALSE] , 2 , function(y) cumsum(y)/sum(y) ))
      }
      
      if ( all(data_g[,j] %in% c(0,1))){
        # factors and NA are coded 0/1, and dW is sorted, so the ks statistic corresponds to last data point with 0
        ks.effect = cbind(ks.effect , dW[sum(data_g[,j]==0),])
      }else{  
        # numeric ks is just the max
        ks.effect = cbind(ks.effect , apply(dW[dups,,drop=FALSE] , 2 , max))
      }
    }
    
    colnames(ks.effect) = paste0(g_name, colnames(data_g))
    
    return(ks.effect) 
  })
  
  out <- do.call("cbind", ks.effect.list)

  return(out)
   
}

