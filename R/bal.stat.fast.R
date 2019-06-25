### calculate weighted balance statistic
bal.stat.fast <- function(data,vars=NULL,treat.var,w.all, sampw, 
                          get.means=TRUE,
                          get.ks=TRUE,
                          na.action="level", estimand, multinom, fillNAs = FALSE,
                          ks.exact=NULL,
                          moderator.var = NULL)
{
  if(is.null(vars)) vars <- names(data)[names(data)!=treat.var]
  
  # is.fac   <- sapply(data[,vars,drop=FALSE],is.factor)
  # fac      <- vars[is.fac]
  # not.fac  <- vars[!is.fac]
  # 
  # ret <- vector("list",length(vars))
  # names(ret) <- vars
  
  #   multinom <- FALSE
  
  #sampW <- sampw
  
  ##### Calculate stats for numeric variables
  #ret[vars] <- ps.summary.new3(data=data[,vars,drop=FALSE], 
  # t=data, w=w.all, sampw = sampW,
  # get.means=get.means, get.ks=get.ks,
  # na.action=na.action,
  # estimand=estimand, multinom = multinom, fillNAs = fillNAs)
  bal.data = gen.bal.data(data=data, var.names=vars, mod.name = moderator.var)
  factor.vars = bal.data$factor.vars
  numeric.vars = bal.data$numeric.vars
  moderator.var = bal.data$moderator.var
  moderator.col <- which(colnames(bal.data$bal.data) == moderator.var)
  treat = data[,treat.var]
  
  if (!is.null(moderator.var)) {
    moderator = data[,moderator.var]
    groups <- c("all", levels(moderator))
    # data <- data[,!colnames(data) %in% moderator.var, drop = FALSE]
  }
  else {
    groups <- ""
  }
  
  res.list <- lapply(groups, function(g) {
    if (g %in% c("", "all")) {
      in. <- rep(TRUE, nrow(data))
      data_g <- data
      bal.data_g <- bal.data$bal.data[in.,!colnames(bal.data$bal.data) %in% moderator.var, drop = FALSE]
      vars_g <- vars
      factor.vars_g <- factor.vars
      numeric.vars_g <- numeric.vars
      if (g == "") g_name <- ""
      else g_name <- paste0(g, "_")
    }
    else {
      in. <- moderator == g
      #Drop moderator variable from subgroup balance
      data_g <- data[in.,!colnames(data) %in% paste(moderator.var, levels(moderator), sep = ":"), drop = FALSE]
      bal.data_g <- bal.data$bal.data[in.,!colnames(bal.data$bal.data) %in% c(moderator.var, paste(moderator.var, levels(moderator), sep = ":")), drop = FALSE]
      vars_g <- vars[vars != moderator.var]
      factor.vars_g <- factor.vars[factor.vars != moderator.var]
      numeric.vars_g <- numeric.vars
      g_name <- paste0(g, "_")
    }
    
    treat_g <- treat[in.]
    w.all_g <- w.all[in.]

    ret <- calcES(data=bal.data_g, treat=treat_g, w=as.matrix(w.all_g), numeric.vars=numeric.vars, estimand=estimand, multinom=multinom, sw=sampw[in.], get.means=TRUE, moderator.var = NULL)
    ret.ks <- t(calcKS(data=bal.data_g, w=as.matrix(w.all_g), treat=treat_g, multinomATE=(estimand=="ATE" && multinom), sw=sampw[in.], moderator.var = NULL))
    colnames(ret.ks) = "ks"
    
    
    # ks.p = NULL
    # ret.test = NULL

    ret.test.factor.list <- lapply(factor.vars_g, function(var){
      x = data_g[in.,var]
      if((sum(is.na(x)) > 0) && (na.action %in% c("level","lowest"))){
        x <- factor(x, levels = c(levels(x), "<NA>"))
        x[is.na(x)] <- "<NA>"
      }
      
      test = suppressWarnings(surveytable(~x+t,data=data.frame(x=x,t=treat_g,w=w.all_g)))
      
      stat <- c(test$statistic, rep(NA, nrow(test$p) - 1))
      pval <- c(test$p.value, rep(NA, nrow(test$p) -1))
      test <- cbind(stat,pval)
      colnames(test) <- c("stat","pval")
      rownames(test) <-  paste(var , levels(x) , sep=":")
      return(test)
    })
    
    # only want to run t-test for those that arent a factor
    ## surveyglm always uses total sample, not number of non-NA
    N = nrow(data_g)
    X = cbind(1 , treat_g)
    WX = sweep( X , 1 , w.all_g , "*")
    Di = crossprod( WX , X)
    ess1 <- sapply(split(w.all_g,treat_g), function(y){sum(y)^2/sum(y^2)})		
    
    ret.test.numeric.list <- lapply(numeric.vars_g, function(var){
      x = bal.data$bal.data[in.,var]
      index = is.na(x)
      x[index] = 0
      #w = w.all
      #w[!index] = 0
      #WX[!index,]=0
      # should there be an if statement?
      D= solve(Di - crossprod( WX[index,] , X[index,]) ) # solve(crossprod( WX , X)) # solve(t(X)%*%diag(w.all*index)%*%X)
      beta = D %*% crossprod( WX , x ) # D %*% t(WX) %*% x     # D%*%t(X)%*%diag(w.all*index)%*%x
      resid = x - X%*%beta
      resid[index] = 0
      a = sweep(WX,1,resid,"*") # sweep(X,1,resid*w.all*index,"*")
      a = beta[2]/ sqrt(crossprod(a%*%t(D))[2,2]*N/(N-1) ) # sqrt( (D%*%t(a)%*%a%*%t(D)*N/(N-1))[2,2] )
      return(matrix(c(a, 2*pt(-abs(a), df=N-2)), nrow=1, dimnames = list(var, c("stat","pval"))))
    })
    
    ks.p.list <- lapply(numeric.vars_g, function(var){
      x = bal.data$bal.data[in.,var]
      index = is.na(x)
      x[index] = 0
      if (any(index)){
        ess <- sapply(split(w.all_g[!index],treat_g[!index]), function(y){sum(y)^2/sum(y^2)})	 # sapply(split(w,treat), function(y){sum(y)^2/sum(y^2)})		
      }else{
        ess <- ess1
      }
      if ( ifelse( is.null(ks.exact) , prod(ess)<10000 , ks.exact) ){
        #PVAL <- 1 - .Call(C_pSmirnov2x, STATISTIC, n.x, n.y)
        pval <- 1- .Call("pSmirnov2x", as.double(ret.ks[var,"ks"]), as.integer(ess[2]), as.integer(ess[1]))
      }else{
        #n <- n.x * n.y / (n.x + n.y)
        #pval <- 2*sum((-1)^(0:4)*exp(-2*(1:5)^2*(ret.ks[var,"ks"])^2*prod(ess)/sum(ess))) # .C("pkstwo", p=as.double(prod(ess)/sum(ess)*ret.ks[var,"ks"]) , as.integer(1), as.double(1e-06) , PACKAGE="twang")$p
        
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
      return(data.frame(ks.pval=pval, row.names = var))
    })
    
    ret.test <- do.call("rbind", c(ret.test.factor.list, ret.test.numeric.list))
    ks.p <- do.call("rbind", ks.p.list)
    
    f <- function(x, y) {transform(merge(x, y, by="row.names", all.x = TRUE), row.names = Row.names, Row.names = NULL)}
    res <- Reduce(f, list(ret, ret.test, ret.ks, ks.p))
    
    ## factors are given their chi2 p-value -- they are missing at this point
    res$ks.pval[grepl(":",rownames(res))] = res$pval[grepl(":",rownames(res))] 
    
    ##### Calculate stats for factor variables
    # ret[is.fac] <- lapply(data[,vars[is.fac],drop=FALSE], ps.summary.new2,
    #                       t=data[,treat.var], w=w.all, sampw = sampW, 
    #                       get.means=get.means, get.ks=get.ks,
    #                       na.action=na.action,
    #                       collapse.by.var=FALSE, estimand=estimand, multinom = multinom, fillNAs = fillNAs)
    
    
    # # this keeps the variables in the same order as vars
    # n.rows <- sapply(ret,nrow)
    # var.levels <- unlist(sapply(ret, rownames))
    # var.names <- rep(names(ret),n.rows)
    # var.names[var.levels!=""] <- paste(var.names[var.levels!=""],
    #                                    var.levels[var.levels!=""],sep=":")
    # 
    # res <- data.frame(matrix(0,nrow=length(var.names), ncol=ncol(ret[[1]])))
    # names(res) <- colnames(ret[[1]])
    # rownames(res) <- var.names
    # 
    # # populate the results table
    # i.insert <- 1
    # for(i in 1:length(ret))
    # {
    #   res[i.insert:(i.insert+nrow(ret[[i]])-1),] <- ret[[i]]
    #   i.insert <- i.insert+nrow(ret[[i]])
    # }
    
    ## maintain order of variables
    var.order = NULL
    
    for (var in vars_g){
      if (is.factor(data_g[,var])){
        var.order = c(var.order , paste(var, levels(data_g[,var]) , sep=":" ) )
      }else{
        var.order = c(var.order , var)
      }
      if (any(is.na(data_g[,var]))){
        var.order = c(var.order , paste0(var, ":<NA>") )
      }
    }
    out <- data.frame(res[var.order,])
    rownames(out) <- paste0(g_name, rownames(out))
    return(out)
  })
  
  res <- list(results=do.call("rbind", res.list))
  return(res)
}

