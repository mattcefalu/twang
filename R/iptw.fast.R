
iptw.fast <- function(formula, data, 
                      # iptw options
                      timeInvariant = NULL, 
                      cumulative = TRUE, 
                      timeIndicators = NULL, 
                      ID = NULL, 
                      priorTreatment = TRUE,
                      # boosting options -- gbm version first, then xgboost name
                      n.trees=10000, nrounds=n.trees,
                      interaction.depth=3, max_depth=interaction.depth,
                      shrinkage=0.01, eta=shrinkage , 
                      bag.fraction = 1.0, subsample=bag.fraction,
                      params=NULL,
                      perm.test.iters=0,
                      print.level=2,       
                      verbose=TRUE,
                      stop.method = c("es.max"), 
                      sampw = NULL, 
                      ks.exact=NULL,
                      version="gbm",
                      tree_method="hist",
                      n.keep = 1,
                      n.grid = 25,
                      ...){
   if(!is.list(formula) & is.null(timeIndicators)) stop("\"formula\" must be a list with length equal to the number of time points (wide data format), or timeIndicators must be specified (long data format).")
   
   if(length(stop.method) > 1) stop("Only one stop.method is allowed at one time.") 
   
   isLong <- !is.list(formula)      
   estimand <- NULL
   
   if(isLong){
      if(length(grep(".time.", attr(terms(formula), "term.labels"))) > 0){
         stop("For long-format data specifications, \".time.\" is not allowed in variable names in the \n
              right-hand side of the formula. Please rename: \n", attr(terms(formula), "term.labels")[grep(".time.", attr(terms(formula), "term.labels"))])
      }
      }
   else{
      wideDat <- NULL
   }
   
   if(isLong){
      if(is.null(ID)) warning("Using long data format specification without IDs. The function assumes \n that the ordering of subjects in 'data' is the same at all time points.")
      tvCov <- attr(terms(formula), "term.labels")
      if(!is.null(timeInvariant)) tiCov <- attr(terms(timeInvariant), "term.labels")
      tb <- table(timeIndicators)
      tb <- tb[tb > 0]
      if(!(max(tb) == min(tb))) stop("Outcomes must be observed at all times for all subjects.")
      unqTimes <- sort(unique(timeIndicators))
      formList <- vector(mode = "list", length = length(unqTimes))
      hldDat <- subset(data, timeIndicators == unqTimes[1])
      trtVar <- all.vars(formula)[1]
      hldDat <- hldDat[,names(hldDat) %in% c(tvCov, trtVar)]
      hdNm <- names(hldDat)
      names(hldDat) <- paste(hdNm, ".time.", unqTimes[1], sep = "")
      trtVarLong <- paste(trtVar, ".time.", unqTimes[1], sep = "")
      tvCovLong <- paste(tvCov, ".time.", unqTimes[1], sep = "")
      formList[[1]] <- as.formula(paste(trtVarLong, paste(tvCovLong, collapse = " + "), sep= "~"))
      wideDat <- hldDat
      for(i in 2:length(unqTimes)){
         hldDat <- subset(data, timeIndicators == unqTimes[i])
         hldDat <- hldDat[,names(hldDat) %in% c(tvCov, trtVar)]
         names(hldDat) <- paste(hdNm, ".time.", unqTimes[i], sep = "")
         wideDat <- cbind(wideDat, hldDat)
         trtVarLong <- paste(trtVar, ".time.", unqTimes[i], sep = "")
         tvCovLong <- paste(tvCov, ".time.", unqTimes[i], sep = "")
         formList[[i]] <- as.formula(paste(trtVarLong, paste(tvCovLong, collapse = " + "), sep= "~"))
      }
      dt2 <- wideDat
      
   }
   else{
      unqTimes <- 1:length(formula)
      formList <- formula
      dt2 <- data
   }
   
   if(! is.null(timeInvariant)){
      invTerms <- attr(terms(timeInvariant), "term.labels")
      for(i in 1:length(formula)){
         currTerms <- union(attr(terms(formList[[i]]), "term.labels"), invTerms)
         #formList[[i]] <- as.formula(paste(all.vars(formList[[i]])[1], paste(invTerms, collapse = " + "), sep = "~"))
         formList[[i]] <- as.formula(paste(all.vars(formList[[i]])[1], paste(currTerms, collapse = " + "), sep = "~"))
         
      }
   }
   
   
   nFits <- length(formList)
   
   if(priorTreatment) {
      for(i in 2:nFits){
         oldTerms <- attr(terms(formList[[i]]), "term.labels")
         allTerms <- c(oldTerms, attr(terms(formList[[i-1]]), "variables")[[2]])
         allTerms <- allTerms[!duplicated(allTerms)]
         formList[[i]] <- as.formula(paste(all.vars(formList[[i]])[1], paste(allTerms, collapse = " + "), sep = "~"))
      }
   }		
   
   if(cumulative){
      for(i in 2:nFits){
         oldTerms <- attr(terms(formList[[i-1]]), "term.labels")
         allTerms <- union(oldTerms, attr(terms(formList[[i]]), "term.labels"))
         formList[[i]] <- as.formula(paste(all.vars(formList[[i]])[1], paste(allTerms, collapse = " + "), sep = "~"))
      }
   }
   
   
   
   formula <- formList	
   
   y1 <- model.extract(model.frame(formula[[1]], data = dt2) , component = "response")
   if(all(y1 %in% c(0,1))) fitType <- 1
   else if(is.factor(y1)) fitType <- 2
   else stop("Treatment variable must either be 0/1 or a factor.")
   
   psList <- vector(mode = "list", length = length(formula))
   
   
   for(i in 1:nFits){
      if(fitType == 1) psList[[i]] <- ps.fast(formula[[i]], data = dt2, 
                                              # boosting options
                                              n.trees=nrounds,                 
                                              interaction.depth=max_depth,
                                              shrinkage=eta,
                                              bag.fraction = subsample,
                                              params=params,
                                              perm.test.iters=perm.test.iters,
                                              print.level=print.level,       
                                              verbose=verbose,
                                              stop.method = stop.method, 
                                              sampw = sampw, 
                                              ks.exact=ks.exact,
                                              version=version,
                                              tree_method=tree_method,
                                              n.keep = n.keep,
                                              n.grid = n.grid,
                                              estimand = "ATE", ...)	
      if(fitType == 2) psList[[i]] <- mnps.fast(formula[[i]], data = dt2, 
                                                # boosting options
                                                n.trees=nrounds,                 
                                                interaction.depth=max_depth,
                                                shrinkage=eta,
                                                bag.fraction = subsample,
                                                params=params,
                                                perm.test.iters=perm.test.iters,
                                                print.level=print.level,       
                                                verbose=verbose,
                                                stop.method = stop.method, 
                                                sampw = sampw, 
                                                ks.exact=ks.exact,
                                                version=version,
                                                tree_method=tree_method,
                                                n.keep = n.keep,
                                                n.grid = n.grid,
                                                estimand = "ATE", ...)	
   }
   
   outObj <- list(psList = psList, stopMethods = stop.method, nFits = nFits, 
                  uniqueTimes = unqTimes)
   
   if(fitType ==1) class(outObj) <- "iptw"
   if(fitType == 2) class(outObj) <- "mniptw"
   return(outObj)
   
   }