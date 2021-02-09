### calculate weighted balance statistic
bal.stat.fast <- function(data,vars=NULL,treat.var,w.all, sampw, 
                     get.means=TRUE,
                     get.ks=TRUE,
                     na.action="level", estimand, multinom, fillNAs = FALSE,
                     ks.exact=NULL)
{
  if(is.null(vars)) vars<-names(data)[names(data)!=treat.var]
  
  # make Row.names null to avoid cran check warning about binding 
  Row.names = NULL
  
  ##### Calculate stats for numeric variables
  bal.data = gen.bal.data(data=data , var.names=vars )
  factor.vars = bal.data$factor.vars
  numeric.vars = bal.data$numeric.vars
  
  treat = data[,treat.var,drop=TRUE]
  
  ret = calcES(data=bal.data$bal.data,treat=treat,w=as.matrix(w.all),numeric.vars=numeric.vars,estimand=estimand,multinom=multinom,sw=sampw,get.means=TRUE)
  ret.ks = t(calcKS(data=bal.data$bal.data,w=as.matrix(w.all),treat=treat , multinomATE=(estimand=="ATE" & multinom) , sw=sampw ))
  colnames(ret.ks) = "ks"
  
  ks.p = NULL
  ret.test = NULL

  for (var in factor.vars){
    x = data[,var]
    if((sum(is.na(x)) > 0) && (na.action %in% c("level","lowest"))){
      x <- factor(x, levels = c(levels(x), "<NA>"))
      x[is.na(x)] <- "<NA>"
    }
    test = suppressWarnings(surveytable(~x+t,data=data.frame(x=x,t=treat,w=w.all)))
  
    stat <- c(test$statistic, rep(NA, nrow(test$p) - 1))
    pval <- c(test$p.value, rep(NA, nrow(test$p) -1))
    test = cbind(stat,pval)
    rownames(test) =  paste(var , levels(x) , sep=":")
    ret.test = rbind(ret.test , test)
  }
  
  # only want to run t-test for those that arent a factor
  ## surveyglm always uses total sample, not number of non-NA
  N = nrow(data)
  X = cbind(1 , treat)
  WX = sweep( X , 1 , w.all , "*")
  Di = crossprod( WX , X)
  ess1 <- sapply(split(w.all,treat), function(y){sum(y)^2/sum(y^2)})		
  for (var in setdiff(colnames(bal.data$bal.data) , rownames(ret.test)) ){
      x = bal.data$bal.data[,var]
      index = is.na(x)
      x[index] = 0
      D= solve(Di - crossprod( WX[index,] , X[index,]) ) 
      beta = D %*% crossprod( WX , x ) 
      resid = x - X%*%beta
      resid[index] = 0
      a = sweep(WX,1,resid,"*") #
      a = beta[2]/ sqrt(crossprod(a%*%t(D))[2,2]*N/(N-1) ) 
      ret.test = rbind(ret.test ,matrix( c(a , 2*pt(-abs(a) ,df=N-2 )) , nrow=1 , dimnames = list(var , c("stat","p"))) )
      
      # ks stat pvalue
      if (any(index)){
        ess <- sapply(split(w.all[!index],treat[!index]), function(y){sum(y)^2/sum(y^2)})	 
      }else{
        ess <- ess1
      }
      if ( ifelse( is.null(ks.exact) , prod(ess)<10000 , ks.exact) ){
        pval <- 1- .Call("pSmirnov2x", as.double(ret.ks[var,"ks"]), as.integer(ess[2]), as.integer(ess[1]))
      }else{
         ## copied from ks.test
         pkstwo <- function(x, tol = 1e-06) {
            if (is.numeric(x)) 
               x <- as.double(x)
            else stop("argument 'x' must be numeric")
            p <- rep(0, length(x))
            p[is.na(x)] <- NA
            IND <- which(!is.na(x) & (x > 0))
            if (length(IND)) 
               p[IND] <- .Call("pKS2", p = x[IND], tol)
            p
         }
         pval <- 1 - pkstwo(sqrt(prod(ess)/sum(ess))*ret.ks[var,"ks"])
      }
      ks.p = rbind(ks.p , data.frame( ks.pval=pval ,row.names = var))
  }
  
  res = transform(merge(ret , ret.test , by="row.names"),row.names=Row.names , Row.names=NULL)
  res = transform(merge(res , ret.ks , by="row.names"),row.names=Row.names , Row.names=NULL)
  res = transform(merge(x=res , y=ks.p , by="row.names" , all.x=T),row.names=Row.names , Row.names=NULL)
  
  ## factors are given their chi2 p-value -- they are missing at this point
  if (length(factor.vars)>0){
    res$ks.pval[grepl(paste0("^(",paste0(factor.vars,collapse = "|"),"):") , rownames(res))] = res$pval[grepl(paste0("^(",paste0(factor.vars,collapse = "|"),"):") , rownames(res))] 
  }

  ##### Calculate stats for factor variables
  ## maintain order of variables
  var.order = NULL
  for (var in vars){
    if (is.factor(data[,var])){
      var.order = c(var.order , paste(var, levels(data[,var]) , sep=":" ) )
    }else{
      var.order = c(var.order , var)
    }
    if (any(is.na(data[,var]))){
      var.order = c(var.order , paste0(var, ":<NA>") )
    }
  }
  
  res <- list(results=data.frame(res[var.order,]))
  return(res)
}

