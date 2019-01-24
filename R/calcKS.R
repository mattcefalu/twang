calcKS <- function(data, w , treat , multinomATE=FALSE, sw=NULL ){
  
  W1 = sweep(w , 1 , treat , "*") #apply(w,2,function(x) x*treat )
  if (multinomATE){
    W0 = as.matrix(sw)
  }else{
    W0 = sweep(w , 1 , 1-treat , "*") #apply(w,2,function(x) x*(1-treat))
  }
  
  ks.effect = NULL
  for (j in 1:ncol(data)){
    # order of the data, dropping NA. Only numeric variables have NAs.
    index = order( data[,j] , na.last = NA)
    
    if (multinomATE){
      # handles the fact that W0 is a vector when multinomATE is true
      dW = abs( sweep(apply(W1[index,,drop=FALSE] , 2 , function(y) cumsum(y)/sum(y) ) , 1 ,  cumsum(sw[index])/sum(sw[index]) , "-" ) )
    }else{
      dW = abs(apply(W1[index,,drop=FALSE] , 2 , function(y) cumsum(y)/sum(y) ) - apply(W0[index,,drop=FALSE] , 2 , function(y) cumsum(y)/sum(y) ))
    }
    
    if ( all(data[,j] %in% c(0,1))){
      # factors and NA are coded 0/1, and dW is sorted, so the ks statistic corresponds to last data point with 0
      ks.effect = cbind(ks.effect , dW[sum(data[,j]==0),])
    }else{  
      # numeric ks is just the max
      ks.effect = cbind(ks.effect , apply(dW , 2 , max))
    }
  }
  colnames(ks.effect) = colnames(data)
    
  return(ks.effect)  
}

