calcES <- function(data,treat,w,numeric.vars,estimand,multinom=FALSE,sw=NULL,get.means=FALSE, moderator.var = NULL){

  # with get.means, only one weight can be specified 
  if (get.means){
    if (  ifelse( is.null(ncol(w)) , FALSE , ncol(w)>1) ){
      stop("Only a single weight can be specified with get.means=TRUE")
    }
  }
  
  if (!is.null(moderator.var)) {
    moderator = data[,moderator.var]
    groups <- c("all", levels(moderator))
    data <- data[,!colnames(data) %in% moderator.var, drop = FALSE]
  }
  else {
    groups <- ""
  }
  
  #Get list of std.effects for each moderator group
  std.effect.list <- lapply(groups, function(g) {
    if (g %in% c("", "all")) { #overall sample
      in. <- rep(TRUE, nrow(data))
      data_ <- data
      if (g == "") g_name <- ""
      else g_name <- paste0(g, "_")
    }
    else { #subgroup-specific sample
      in. <- moderator == g
      #Drop moderator variable from group-specific balance
      data_ <- data[,!colnames(data) %in% paste(moderator.var, levels(moderator), sep = ":"), drop = FALSE]
      g_name <- paste0(g, "_")
    }
    
    #in. means "in the group under analysis". With no moderator and for the overall sample, this is TRUE for everyone.
    
    N = sum(in.)
    N1 = sum(treat[in.])
    N0 = N-N1
    
    # separate out treatment group for easy matrix operations
    w1 = w[in. & treat==1,,drop = FALSE]                  # subset to treated weight matrix
    data1 = as.matrix(data_[in. & treat==1,,drop = FALSE]) # subset to treated data
    index.1 = !is.na(data1)            # only numeric variable will have missing data, which are dropped from ES calculation.
    den1 = t(w1) %*% index.1             # sum of weights changes depending on missing data -- only numeric variables should have missing values in data
    data1[!index.1] = 0                # change missing values to 0 so they dont contribute to the sum and dont mess up matrix operations
    m1 = t(w1) %*% data1 / den1          # calculate mean among treated
    
    
    
    # m0 and the SD depends on estimand
    if (estimand=="ATE"){
      # need variance based only on sw
      data_g = as.matrix(data_[in.,,drop = FALSE])
      sw_g <- sw[in.]
      index = !is.na(data_g)
      data_g[!index] = 0         # change missing values to 0 so they dont contribute to the sum and dont mess up matrix operations
      den = sw_g %*% index
      m.sw = sw_g %*% data_g / den
      
      # pooled variance estimate for ATE
      var.p = sw_g %*% data_g ^ 2 / den - (m.sw)^2
      var.p = ifelse(var.p==0 , NA_real_ , var.p)
      var.p[,numeric.vars] =  var.p[,numeric.vars,drop = FALSE] * N / (N-1)
      
      if (get.means){
        # calculate E[X^2] among treated
        m1.2 =  t(w1)%*%(data1^2) / den1 
        var1=(m1.2-m1^2) 
        var1 = ifelse(var1==0 , NA_real_ , var1)
        var1[,numeric.vars] =  var1[,numeric.vars,drop = FALSE] * N1 / (N1-1)
      }
      
      if (multinom){
        # multinom ATE: m0 is full population mean
        m0 = m.sw
        
        if (get.means){
          var0=var.p
        }
        
        # this just handles that m0 and var.p are vectors when using multinom option
        std.effect = sweep(sweep( m1 , 2 , sqrt(var.p) , "/") , 2 , m0/sqrt(var.p) , "-")
      }else{
        # same steps as with treated, but now for controls
        w0 = w[in. & treat==0,,drop = FALSE]
        data0 = as.matrix(data_[in. & treat==0,,drop = FALSE])
        index.0 = !is.na(data0)
        den0 = t(w0)%*%index.0
        data0[!index.0] = 0
        
        # calculate mean among control
        m0 = t(w0)%*%data0 / den0
        
        if (get.means){
          # calculate E[X^2] among control
          m0.2 =  t(w0)%*%(data0^2) / den0 
          var0=(m0.2-m0^2) 
          var0 = ifelse(var0==0 , NA , var0)
          var0[,numeric.vars] =  var0[,numeric.vars,drop = FALSE] * N0 / (N0-1)
        }
        
        # this just handles that var.p is a vector
        std.effect = sweep( m1-m0 , 2 , sqrt(var.p) , "/")
      }
    }else{ #ATT
      # same steps as with treated, but now for controls
      w0 = w[in. & treat==0,,drop = FALSE]
      data0 = as.matrix(data_[in. & treat==0,,drop = FALSE])
      index.0 = !is.na(data0)
      den0 = t(w0)%*%index.0
      data0[!index.0] = 0
      
      # calculate mean among control
      m0 = t(w0)%*%data0 / den0
      
      # calculate E[X^2] among treated
      m1.2 =  t(w1)%*%(data1^2) / den1 
      var1=(m1.2-m1^2) 
      var1 = ifelse(var1==0 , NA , var1)
      var1[,numeric.vars] =  var1[,numeric.vars,drop = FALSE] * N1 / (N1-1)
      # var1[,] =  var1[,,drop = FALSE] * N1 / (N1-1)
      
      if (get.means){
        # calculate E[X^2] among control
        m0.2 =  t(w0)%*%(data0^2) / den0 
        var0=(m0.2-m0^2) 
        var0 = ifelse(var0==0 , NA , var0)
        var0[,numeric.vars] =  var0[,numeric.vars,drop = FALSE] * N0 / (N0-1)
        # var0[,] =  var0[,,drop = FALSE] * N0 / (N0-1)
      }
      
      std.effect = (m1-m0)/sqrt(var1)
    }
    
    
    if (get.means){
      ## surveyglm always uses total sample, not number of non-NA
      # N = nrow(data)
      # X = cbind(1 , treat)
      # index = !is.na(data)
      # t.n <- NULL
      # for (var in numeric.vars){
      #   x = data[,var]
      #   D=solve(t(X)%*%diag(w*index[,var])%*%X)
      #   beta = D%*%t(X)%*%diag(w*index[,var])%*%x
      #   resid = x - X%*%beta
      #   a = sweep(X,1,resid*w*index[,var],"*")
      #   a = beta[2]/sqrt( (D%*%t(a)%*%a%*%t(D)*N/(N-1))[2,2] )
      #   t.n = rbind(t.n ,c(a , 1-2*pt(-abs(a) ,df=N-1 )) )
      # }
      # rownames(t.n) = numeric.vars
      res = t(rbind(tx.mn=m1 , tx.sd=sqrt(var1), ct.mn=m0 , ct.sd=sqrt(var0) , std.eff.sz=std.effect ))
      colnames(res) = c("tx.mn", "tx.sd",  "ct.mn", "ct.sd", "std.eff.sz")
      rownames(res) <- paste0(g_name, colnames(std.effect))
      #res[numeric.vars,c("tx.sd","ct.sd")] =  res[numeric.vars,c("tx.sd","ct.sd")]*sqrt(N/(N-1))
      return(res)
    }else{
      colnames(std.effect) <- paste0(g_name, colnames(std.effect))
      return(std.effect)
    }
    
  })
  
  if (get.means){
    out <- do.call("rbind", std.effect.list)
  }
  else {
    out <- do.call("cbind", std.effect.list)
  }

  return(out)
}