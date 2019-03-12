## require(twang); data(lalonde); ps.lalonde <- ps(treat~age + educ + black + hispan + nodegree + married + re74 + re75, data = lalonde, stop.method = c("es.max", "es.mean"),estimand = "ATT", n.trees = 5000, verbose = FALSE)
ps.fast<-function(formula ,
             data,                         # data
             params=NULL,
             n.trees=10000,                 # gbm options
             interaction.depth=3,
             shrinkage=0.01,
             bag.fraction = 1.0,
             n.minobsinnode = 10,
             perm.test.iters=0,
             print.level=2,                 # direct optimizer options
             #iterlim=1000,
             verbose=TRUE,
             estimand="ATE", 
             stop.method = c("ks.mean", "es.mean"), 
             sampw = NULL, multinom = FALSE,
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
             	
	

   if(is.null(sampw)){
      sampW <- rep(1, nrow(data))
   }else{ 
      sampW <- unlist(sampw, use.names=F) 
   }
   
   if (!is.null(n.grid)){
      n.grid.ks = n.grid
      n.grid.es = n.grid
   }
   
   if ( n.trees/n.keep< max(n.grid.ks,n.grid.es) ){
      stop('n.tress must be at least max(n.grid.ks,n.grid.es) times n.keep')
   }
	
   type <- alert <- NULL
	
	dots <- list(...)
	if(!is.null(dots$plots)){
	   warning("From version 1.2, the plots argument has been removed from ps(). \nPlease use the plot() function instead.")
	}
	
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
   
   # going to let the boosters handle the parsing of the equation. we only need to extract the variable names
   # var.names = all.vars(formula)[-1]
   # treat.var <- all.vars(formula)[1]
   
   if(length(var.names) < 2) stop("At least two variables are needed in the right-hand side of the formula.\n")
   
   treat.var <- as.character(formula[[2]])
      
   stop.method.names <- sapply(stop.method,function(x){x$name})
   i <- which(sapply(stop.method,function(x){x$direct}))
   if(length(i)>0) cat(paste("*** WARNING: Stop method",stop.method.names[i],"involves direct\noptimization, which is very time consuming, especially\nfor datasets with more than a few hundred cases. Consider\nusing another stop.method or be prepared for a very long wait. ***\n\n"))
   desc <- vector("list",1+length(stop.method))
   names(desc) <- c("unw",stop.method.names)
 
   # alert.stack collects all the warnings
   alerts.stack <- textConnection("alert","w")

   # fit the propensity score model
   if(verbose) cat("Fitting boosted model\n")
   
   if (booster=="gbm"){
      # need to reformulate formula to use this environment
      form <- paste(deparse(formula, 500), collapse="") 
   
      gbm1 <-gbm(formula(form),
                 data = data,
                 weights=sampW,
                 distribution = "bernoulli",
                 n.trees = n.trees,
                 interaction.depth = interaction.depth,
                 n.minobsinnode = n.minobsinnode,
                 shrinkage = shrinkage,
      ### bag.fraction was 0.5.  revised 101210
                 bag.fraction = bag.fraction,
                 train.fraction = 1,
                 verbose = verbose,
                 keep.data = FALSE)
      
      # predict propensity scores for all iterations
      iters = (1:n.trees)[(1:n.trees)%%n.keep==0]
      ps = plogis(predict(gbm1 , newdata = data , n.trees=iters))
      
   }else{
      sparse.form = reformulate(c("-1",var.names))
      sparse.data = model.Matrix(sparse.form , model.frame(terms(sparse.form),data=data[,var.names] , na.action=na.pass), drop.unused.levels=T , sparse = T)
      #dtrain <- xgb.DMatrix(train$data, label=data[,treat.var])
      
      if (is.null(params)){
         params = list(eta = shrinkage , max_depth = interaction.depth , subsample = bag.fraction , min_child_weight=n.minobsinnode , objective = "binary:logistic")
      }else{
         if( is.null(params$objective) ){
            params$objective = "binary:logistic"
         }
      }
         # if (save.propensities){
         #    gbm1 <- xgboost(data=sparse.data , label=data[,treat.var], params=params, tree_method = tree_method, 
         #                    nrounds=n.trees , verbose=verbose , 
         #                    callbacks=list(cb.print.evaluation(),save.model(save_period=n.eval.propensity,save_name=file))) 
         #    ps = readRDS(file=file)
         # }else{
            gbm1 <- xgboost(data=sparse.data , label=data[,treat.var], params=params, tree_method = tree_method, 
                            feval=pred.xgboost , nrounds=n.trees , verbose=verbose , weight = sampW, 
                            callbacks=list(cb.print.evaluation() , cb.evaluation.log(n.keep = n.keep)))
            iters = (1:n.trees)[(1:n.trees)%%n.keep==0]
            ps = as.matrix(gbm1$evaluation_log)
            #ps= #do.call(cbind,gbm1$evaluation_log$train_pred)
         #}
      # predict propensity scores for all iterations
      #iters = 1:n.trees
      #ps = plogis(predict(gbm1 , newdata = data , n.trees=iters))
   }

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
   
   if (estimand=="ATE"){
     W = 1 / (1-ps) + data[,treat.var]*( 1/ps - 1/(1-ps)  )
   }else{
     W = ps/(1-ps) + data[,treat.var]*( 1 - ps/(1-ps)  )
   }
   
   # adjust for sampling weights 
   W <- sweep(W,1,sampW,"*")
   
   # setup balance data
   bal.data = gen.bal.data(data=data , var.names=var.names )
   factor.vars = bal.data$factor.vars
   numeric.vars = bal.data$numeric.vars
   bal.data = bal.data$bal.data
   
   # 25 point grid of iterations
   #iters.25 <- round(seq( 1 , length(iters)  ,length=25))
   iters.grid.ks <- round(seq( 1 , length(iters)  ,length=n.grid.ks))
   iters.grid.es <- round(seq( 1 , length(iters)  ,length=ifelse(is.null(n.grid.es) ,length(iters) ,n.grid.es )))
   
   std.effect = ks.effect = NULL
   if ( any(grepl("es.",stop.method.names)) ){
      if(verbose) cat("Calculating standardized differences\n")
      std.effect = calcES(data=bal.data, w=W[,iters.grid.es] , treat=data[,treat.var],numeric.vars = numeric.vars , estimand=estimand , multinom=multinom , sw=sampW)
      
      if (!is.null(n.grid.es)){
         iters.es.mean = iters.es.max = NULL
         if ( any(grepl("es.mean", stop.method.names)) ){
            i <- which.min(apply(abs(std.effect),1,mean, na.rm=T)) +c(-1,1)
            i[1] <- iters.grid.es[max(1,i[1])]
            i[2] <- iters.grid.es[min(length(iters.grid.es),i[2])]
            iters.es.mean =  which(iters <= iters[i[2]] & iters >= iters[i[1]]) #iters.25[i[1]]:iters.25[i[2]]
         }
         if ( any(grepl("es.max", stop.method.names)) ){
            i <- which.min(apply(abs(std.effect),1,max, na.rm=T)) +c(-1,1)
            i[1] <- iters.grid.es[max(1,i[1])]
            i[2] <- iters.grid.es[min(length(iters.grid.es),i[2])]
            iters.es.max =  which(iters <=iters[i[2]] & iters >=iters[i[1]])  #iters.25[i[1]]:iters.25[i[2]]
         }
         
         # combine the intervals to do a finer search
         iters.es = unique(c(iters.es.max , iters.es.mean))
         
         # save the grid
         balance.es = std.effect
         iters.es.plot = iters.grid.es
         
         std.effect = calcES(data=bal.data, w=W[,iters.es] , treat=data[,treat.var],numeric.vars = numeric.vars , estimand=estimand , multinom=multinom , sw=sampW)
         #ks.effect = calcKS(data=bal.data,w=W[,iters.es],treat=data[,treat.var])
         #colnames(std.effect) = colnames(bal.data)
      }else{
         iters.es = iters.grid.es
         iters.es.plot = iters.grid.ks
         balance.es = std.effect[iters.grid.ks,]
      }
   }
   if ( any(grepl("ks.",stop.method.names)) ){
      if(verbose) cat("Calculating Kolmogorov–Smirnov statistics\n")
     
      # find the optimal interval for ks.mean and ks.max based on 25 point grid     
      ks.effect = calcKS(data=bal.data,w=W[,iters.grid.ks],treat=data[,treat.var] , multinomATE=(estimand=="ATE" & multinom) , sw=sampW)
     
      iters.ks.mean = iters.ks.max = NULL
      if ( any(grepl("ks.mean", stop.method.names)) ){
         i <- which.min(apply(abs(ks.effect),1,mean, na.rm=T)) +c(-1,1)
         i[1] <- iters.grid.ks[max(1,i[1])]
         i[2] <- iters.grid.ks[min(length(iters.grid.ks),i[2])]
         iters.ks.mean =  which(iters <= iters[i[2]] & iters >= iters[i[1]]) #iters.25[i[1]]:iters.25[i[2]]
      }
      if ( any(grepl("ks.max", stop.method.names)) ){
         i <- which.min(apply(abs(ks.effect),1,max, na.rm=T)) +c(-1,1)
         i[1] <- iters.grid.ks[max(1,i[1])]
         i[2] <- iters.grid.ks[min(length(iters.grid.ks),i[2])]
         iters.ks.max =  which(iters <=iters[i[2]] & iters >=iters[i[1]])  #iters.25[i[1]]:iters.25[i[2]]
      }
     
      # combine the intervals to do a finer search
      iters.ks = unique(c(iters.ks.max , iters.ks.mean))
      
      # save the 25 point grid ks
      balance.ks = ks.effect
      
      ks.effect = calcKS(data=bal.data,w=W[,iters.ks],treat=data[,treat.var] , multinomATE=(estimand=="ATE" & multinom) , sw=sampW)

      colnames(ks.effect) = colnames(bal.data)
   }

   # holds balance stats needed for plots
   balance = list()
   
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
         opt = list(minimum= iters.es[which.min(apply(abs(std.effect),1,mean, na.rm=T))] )
         balance = c(balance , list(es.mean = data.table(iter=iters[iters.es.plot] , value=apply(abs(balance.es),1,mean, na.rm=T)) )) #cbind(balance , apply(abs(std.effect[iters.25,]),1,mean, na.rm=T) )
      }
      if ( grepl("es.max",stop.method.names[i.tp]) ){
         opt = list(minimum= iters.es[which.min(apply(abs(std.effect),1,max, na.rm=T))] )
         #balance = cbind(balance , apply(abs(std.effect[iters.25,]),1,max, na.rm=T) )
         balance = c(balance , list(es.max = data.table(iter=iters[iters.es.plot] , value=apply(abs(balance.es),1,max, na.rm=T)) )) #cbind(balance , apply(abs(std.effect[iters.25,]),1,mean, na.rm=T) )
      }
      if ( grepl("ks.mean",stop.method.names[i.tp]) ){
         opt = list(minimum= iters.ks[which.min(apply(abs(ks.effect),1,mean, na.rm=T))] )
         balance = c(balance , list(ks.mean = data.table(iter=iters[iters.grid.ks] , value=apply(abs(balance.ks),1,mean, na.rm=T)) )) 
         #balance = cbind(balance ,  apply(abs(balance.ks),1,mean, na.rm=T) )
      }
      if ( grepl("ks.max",stop.method.names[i.tp]) ){
         opt = list(minimum= iters.ks[which.min(apply(abs(ks.effect),1,max, na.rm=T))] )
         balance = c(balance , list(ks.max = data.table(iter=iters[iters.grid.ks] , value=apply(abs(balance.ks),1,max, na.rm=T)) )) 
         #balance = cbind(balance ,  apply(abs(balance.ks),1,max, na.rm=T) )
      }
      
   
   
      if(verbose) cat("   Optimized at",iters[opt$minimum],"\n")
      if(n.trees-iters[opt$minimum] < 100) warning("Optimal number of iterations is close to the specified n.trees. n.trees is likely set too small and better balance might be obtainable by setting n.trees to be larger.")
      
      # save propensity score weights
      p.s[,i.tp]  <- ps[,opt$minimum]
      w[,i.tp] <- W[ , opt$minimum]
      
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
      
      desc[[tp]]$n.trees <- ifelse(stop.method[[i.tp]]$direct, NA, iters[opt$minimum])
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
                  booster = booster , 
                  version = version, 
#                  plot.info  = plot.info,
                  datestamp  = date(),
                  parameters = terms,
                  alerts     = alert,
                  iters.ks = iters.grid.ks,
                  iters.es = iters.grid.es,
                  balance.ks = ks.effect,
                  balance.es = balance.es,
                  balance = balance,
                  es = std.effect,
                  ks = ks.effect,
                  n.trees = n.trees,
                  data = data)
   class(result) <- "ps"
   return(result)
}

