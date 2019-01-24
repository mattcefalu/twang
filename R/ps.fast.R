## require(twang); data(lalonde); ps.lalonde <- ps(treat~age + educ + black + hispan + nodegree + married + re74 + re75, data = lalonde, stop.method = c("es.max", "es.mean"),estimand = "ATT", n.trees = 5000, verbose = FALSE)
ps<-function(formula = formula(data),
             data,                         # data
             n.trees=10000,                 # gbm options
             interaction.depth=3,
             shrinkage=0.01,
             bag.fraction = 1.0,
             perm.test.iters=0,
             print.level=2,                 # direct optimizer options
             iterlim=1000,
             verbose=TRUE,
             estimand="ATE", 
             stop.method = c("ks.mean", "es.mean"), 
             sampw = NULL, multinom = FALSE,
             ks.exact=NULL,
             ...){
             	
	
#	multinom <- FALSE
	
	if(is.null(sampw)) sampW <- rep(1, nrow(data))
	else sampW <- unlist(sampw, use.names=F)
	
	type <- alert <- NULL
	
	dots <- list(...)
	if(!is.null(dots$plots))
	warning("From version 1.2, the plots argument has been removed from ps(). \nPlease use the plot() function instead.")
	
	stop.method[stop.method == "ks.stat.mean"] <- "ks.mean"
	stop.method[stop.method == "es.stat.mean"] <- "es.mean"
	stop.method[stop.method == "ks.stat.max"] <- "ks.max"
	stop.method[stop.method == "es.stat.max"] <- "es.max"
	stop.method[stop.method == "ks.stat.max.direct"] <- "ks.max.direct"
	stop.method[stop.method == "es.stat.max.direct"] <- "es.max.direct"	
	

	             


if(!(estimand %in% c("ATT","ATE"))) stop("estimand must be either \"ATT\" or \"ATE\".")
allowableStopMethods <- c("ks.mean", "es.mean","ks.max", "es.max","ks.max.direct","es.max.direct")

nMethod <- length(stop.method)
methodList <- vector(mode="list", length = nMethod)
for(i in 1:nMethod){
	if(is.character(stop.method[i])){
		if (!(stop.method[i] %in% allowableStopMethods)){
			print(allowableStopMethods)
			stop("Each element of stop.method must either be one of \nthe above character strings, or an object of the stop.method class.")	
		}
		
	methodList[[i]] <- get(stop.method[i])
	methodName <- paste(stop.method[i], ".", estimand, sep="")
	methodList[[i]]$name <- methodName
	}
	else {
		if (class(stop.method[i]) != "stop.method"){
			print(allowableStopMethods)
			stop("Each element of stop.method must either be one of \nthe above character strings, or an object of the stop.method class.")	
		}
		methodList[[i]] <- stop.method[i]	
	}
}

stop.method <- methodList


   terms <- match.call()
   
   # all this is just to extract the variable names
   mf <- match.call(expand.dots = FALSE)
   m <- match(c("formula", "data"), names(mf), 0)
   mf <- mf[c(1, m)]
   mf[[1]] <- as.name("model.frame")
   mf$na.action <- na.pass
   mf$subset <- rep(FALSE, nrow(data)) # drop all the data
   mf <- eval(mf, parent.frame())
   Terms <- attr(mf, "terms")
   var.names <- attributes(Terms)$term.labels
   
   if(length(var.names) < 2) stop("At least two variables are needed in the right-hand side of the formula.\n")
   
   treat.var <- as.character(formula[[2]])

   # create the desc object. This holds information on variable balance
#   if (class(stop.method)=="stop.method"){
#      stop.method <- list(stop.method)}
      
   stop.method.names <- sapply(stop.method,function(x){x$name})
   i <- which(sapply(stop.method,function(x){x$direct}))
   if(length(i)>0) cat(paste("*** WARNING: Stop method",stop.method.names[i],"involves direct\noptimization, which is very time consuming, especially\nfor datasets with more than a few hundred cases. Consider\nusing another stop.method or be prepared for a very long wait. ***\n\n"))
   desc <- vector("list",1+length(stop.method))
   names(desc) <- c("unw",stop.method.names)


   
   # check sampw
#   if(is.null(sampw)) sampw <- rep(1,nrow(data))
 
   # alert.stack collects all the warnings
   alerts.stack <- textConnection("alert","w")

   # fit the propensity score model
   if(verbose) cat("Fitting gbm model\n")
   # need to reformulate formula to use this environment
   form <- paste(deparse(formula, 500), collapse="") 

   gbm1 <-gbm(formula(form),
              data = data,
              weights=sampW,
              distribution = "bernoulli",
              n.trees = n.trees,
              interaction.depth = interaction.depth,
              n.minobsinnode = 10,
              shrinkage = shrinkage,
   ### bag.fraction was 0.5.  revised 101210
              bag.fraction = bag.fraction,
              train.fraction = 1,
              verbose = verbose,
              keep.data = FALSE)

   if(verbose) cat("Diagnosis of unweighted analysis\n")
   
   if(is.factor(data[,treat.var])) stop("Treatment indicator must be numeric, not a factor")
   
   desc$unw <- desc.wts.fast(data=data[,c(treat.var,var.names)],
                        treat.var=treat.var,
                        w=sampW,
                        sampw = rep(1, nrow(data)), 
                        tp="unw",
                        na.action="level",
                        perm.test.iters=perm.test.iters,
                        verbose=verbose,
                        alerts.stack=alerts.stack,
                        estimand=estimand, multinom = multinom,
                        ks.exact=ks.exact)
   desc$unw$n.trees <- NA

    #balance <- matrix(NA, ncol = nMethod, nrow = 25)
    #balance <- matrix(NA, nc = nMethod, nr = 100)

   
   # predict propensity scores for all iterations
   iters = 1:n.trees
   ps = plogis(predict(gbm1 , newdata = data , n.trees=iters))
   
   if (estimand=="ATE"){
     W = 1 / (1-ps) + data[,treat.var]*( 1/ps - 1/(1-ps)  )
   }else{
     W = ps/(1-ps) + data[,treat.var]*( 1 - ps/(1-ps)  )
   }
   
   # adjust for sampling weights 
   W <- sweep(W,1,sampW,"*")
   
   
   bal.data = gen.bal.data(data=data , var.names=var.names )
   factor.vars = bal.data$factor.vars
   numeric.vars = bal.data$numeric.vars
   bal.data = bal.data$bal.data
   
   # 25 point grid of iterations
   iters.25 <- round(seq(1,gbm1$n.trees,length=25))
   
   std.effect = ks.effect = NULL
   if ( any(grepl("es.",stop.method.names)) ){
    if(verbose) cat("Calculating standardized differences\n")
     std.effect = calcES(data=bal.data, w=W , treat=data[,treat.var],numeric.vars = numeric.vars , estimand=estimand , multinom=multinom , sw=sampW)
   }
   if ( any(grepl("ks.",stop.method.names)) ){
     if(verbose) cat("Calculating Kolmogorov–Smirnov statistics\n")
     
     # find the optimal interval for ks.mean and ks.max based on 25 point grid     
     ks.effect = calcKS(data=bal.data,w=W[,iters.25],treat=data[,treat.var] , multinomATE=(estimand=="ATE" & multinom) , sw=sampW)
     
     iters.ks.mean = iters.ks.max = NULL
     if ( any(grepl("ks.mean", stop.method.names)) ){
        i <- which.min(apply(abs(ks.effect),1,mean, na.rm=T)) +c(-1,1)
        i[1] <- max(1,i[1])
        i[2] <- min(length(iters.25),i[2])
        iters.ks.mean = iters.25[i[1]]:iters.25[i[2]]
     }
     if ( any(grepl("ks.max", stop.method.names)) ){
       i <- which.min(apply(abs(ks.effect),1,mean, na.rm=T)) +c(-1,1)
       i[1] <- max(1,i[1])
       i[2] <- min(length(iters.25),i[2])
       iters.ks.max = iters.25[i[1]]:iters.25[i[2]]
     }
     
     # combine the intervals to do a finer search
     iters.ks = unique(c(iters.ks.max , iters.ks.mean))
     
     # save the 25 point grid ks
     balance.ks = ks.effect
     
     ks.effect = calcKS(data=bal.data,w=W[,iters.ks],treat=data[,treat.var])
     colnames(ks.effect) = colnames(bal.data)
   }

   # holds balance stats needed for plots
   balance = NULL
   
   # allocate space for the propensity scores and weights
   p.s        <- data.frame(matrix(NA,nrow=nrow(data),
                                   ncol=length(stop.method)))
   names(p.s) <- stop.method.names
   w          <- data.frame(matrix(NA,nrow=nrow(data),
                                   ncol=length(stop.method)))
   names(w)   <- stop.method.names
    for (i.tp in 1:nMethod){
       tp <- stop.method.names[i.tp]
       if(verbose) cat("Optimizing with",tp,"stopping rule\n")
   
       # find the optimal values based on the precomputed statistics
       # and save the balance stat to meet plotting requirements
       if ( grepl("es.mean",stop.method.names[i.tp]) ) {
        opt = list(minimum= which.min(apply(abs(std.effect),1,mean, na.rm=T)) )
        balance = cbind(balance , apply(abs(std.effect[iters.25,]),1,mean, na.rm=T) )
       }
       if ( grepl("es.max",stop.method.names[i.tp]) ){
         opt = list(minimum= which.min(apply(abs(std.effect),1,max, na.rm=T)) )
         balance = cbind(balance , apply(abs(std.effect[iters.25,]),1,max, na.rm=T) )
       }
       if ( grepl("ks.mean",stop.method.names[i.tp]) ){
         opt = list(minimum= iters.ks[which.min(apply(ks.effect,1,mean, na.rm=T))] )
         balance = cbind(balance ,  apply(abs(balance.ks),1,mean, na.rm=T) )
       }
       if ( grepl("ks.max",stop.method.names[i.tp]) ){
         opt = list(minimum= iters.ks[which.min(apply(ks.effect,1,max, na.rm=T))] )
         balance = cbind(balance ,  apply(abs(balance.ks),1,max, na.rm=T) )
       }
   
   
      if(verbose) cat("   Optimized at",round(opt$minimum),"\n")
      if(gbm1$n.trees-opt$minimum < 100) warning("Optimal number of iterations is close to the specified n.trees. n.trees is likely set too small and better balance might be obtainable by setting n.trees to be larger.")

      # save propensity score weights
      p.s[,i.tp]  <- ps[,opt$minimum]
      # p.s[,i.tp]  <- predict(gbm1,newdata=data,
      #                        n.trees=round(opt$minimum),
      #                        type="response")

      w[,i.tp] <- W[ , opt$minimum]
      # if (estimand=="ATT") {
      # w[,i.tp] <- p.s[,i.tp]/(1-p.s[,i.tp])
      # w[data[,treat.var]==1,i.tp] <- 1
      # }
      # 
      # if (estimand=="ATE"){
      # w[data[,treat.var]==1,i.tp] <- 1/p.s[data[,treat.var]==1,i.tp]
      # w[data[,treat.var]==0,i.tp] <- 1/(1-p.s[data[,treat.var]==0,i.tp])
      # }
      # 
      # # adjust for sampling weights
      # w[,i.tp] <- w[,i.tp]*sampW 

      # Directly optimize the weights if requested
      # if (stop.method[[i.tp]]$direct){
      #    if (verbose) cat("   Proceeding to direct optimization\n")
      # 
      #    obj   <- nlm(match.fun(stop.method[[i.tp]]$metric),
      #                 p           = log(w[data[,treat.var]==0,i.tp]),
      #                 fscale      = 0.1,
      #                 data        = data,
      #                 vars        = var.names,
      #                 treat.var   = treat.var,
      #                 ndigit      = 3,
      #                 rule.summary = match.fun(stop.method[[i.tp]]$rule.summary),
      #                 na.action   = stop.method[[i.tp]]$na.action,
      #                 print.level = print.level,
      #                 iterlim     = iterlim,
      #                 estimand    = estimand)
      #    if(obj$code %in% c(4,5))
      #       warning("Failed to completely optimize metric directly.")
      # 
      #    w[data[,treat.var]==0,i.tp] <- exp(obj$estimate)
      # }

      if(verbose) cat("Diagnosis of",tp,"weights\n")
      desc[[tp]] <- desc.wts.fast(data[,c(treat.var,var.names)],
                             treat.var    = treat.var,
                             w            = w[,i.tp],
                             sampw = sampW, 
                             tp           = type,
                             na.action    = stop.method[[i.tp]]$na.action,
                             perm.test.iters=perm.test.iters,
                             verbose=verbose,
                             alerts.stack = alerts.stack,
                             estimand       = estimand,
                             multinom = multinom,
                             ks.exact=ks.exact)
      desc[[tp]]$n.trees <- 
         ifelse(stop.method[[i.tp]]$direct, NA, round(opt$minimum))
      

   }


   close(alerts.stack)
   if(verbose) cat(alert,sep="\n")

   result <- list(gbm.obj    = gbm1,
                  treat      = data[,treat.var], # needed for plot.ps
                  treat.var  = treat.var,
                  desc       = desc,
                  ps         = p.s,
                  w          = w,
                  sampw      = sampW, 
                  estimand   = estimand,
#                  plot.info  = plot.info,
                  datestamp  = date(),
                  parameters = terms,
                  alerts     = alert,
                  iters = iters.25,
                  balance = balance,
                  es = std.effect,
                  ks = ks.effect,
                  n.trees = n.trees,
                  data = data)
   class(result) <- "ps"
   return(result)
}

