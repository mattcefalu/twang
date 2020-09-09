#' Compute diagnostics assessing covariates balance.
#'
#' `dx.wts` takes a `ps` object or a set of propensity scores and 
#' computes diagnostics assessing covariates balance.
#'
#' Creates a balance table that compares unweighted and weighted means and 
#' standard deviations, computes effect sizes, and KS statistics to assess the 
#' ability of the propensity scores to balance the treatment and control groups.
#'
#' @param x A data frame, matrix, or vector of propensity score weights or a ps 
#'   object. `x` can also be a data frame, matrix, or vector of 
#'   propensity scores if `x.as.weights=FALSE`.
#' @param data A data frame.
#' @param estimand The estimand of interest: either "ATT" or "ATE".
#' @param vars A vector of character strings naming variables in `data` on 
#'   which to assess balance.
#' @param treat.var A character string indicating which variable in `data`
#'   contains the 0/1 treatment group indicator.
#' @param x.as.weights `TRUE` or `FALSE` indicating whether `x`
#'   specifies propensity score weights or propensity scores. 
#'   Ignored if `x` is a ps object. Default: `TRUE`.
#' @param sampw Optional sampling weights. If `x` is a `ps` object, then the 
#'   sampling weights should have been passed to [ps] and 
#'   not specified here. `dx.wts` will issue a warning if 
#'   `x` is a ps object and `sampw` is also specified.
#' @param perm.test.iters A non-negative integer giving the number of iterations
#'   of the permutation test for the KS statistic. If `perm.test.iters=0`,
#'   then the function returns an analytic approximation to the p-value. This
#'   argument is ignored is `x` is a `ps` object. Setting
#'   `perm.test.iters=200` will yield precision to within 3% if the true
#'   p-value is 0.05. Use `perm.test.iters=500` to be within 2%.
#'
#' @return Returns a list containing
#'   * `treat` The vector of 0/1 treatment assignment indicators.
#'
#' @seealso [ps]
#'
#' @export
dx.wts <- function(x,
                   data,
                   estimand,
                   vars = NULL,
                   treat.var,
                   x.as.weights = TRUE,
                   sampw = NULL,
                   perm.test.iters = 0){


	
	if(is.null(estimand))
	stop("'estimand' must be specified")
	if(!(estimand %in% c("ATE","ATT")))
	stop("'estimand' must be either 'ATT' or 'ATE'")
	
   if (class(x)[1]!="ps"){
      if(!is.vector(x) && !is.matrix(x) && !is.data.frame(x))
         stop("x must be a ps object, a vector, a matrix, or a data frame")
      if(any(x < 0)) stop("x has negative values")
      if(is.null(dim(x))) x <- matrix(x,ncol=1)
      if(nrow(x) != nrow(data)) stop("length(x) != nrow(data)")

      if(x.as.weights){
      	
      	if(is.data.frame(x)) x <- as.matrix(x)
      	
      	if (estimand=="ATT") {
         w <- x
         p.s <- x/(1+x)
         }
         
        if (estimand=="ATE") {
        	w <- x
        	p.s <- x # to create a default length of p.s
        	p.s[data[,treat.var]==0] <- (w[data[,treat.var]==0]-1)/w[data[,treat.var]==0] 
        	p.s[data[,treat.var]==1] <- 1/w[data[,treat.var]==1]
        	}
      }     
  else {
      	if(any(x > 1) && !x.as.weights) stop("x has values greater than 1. With x.as.weights=FALSE x should be a vector of propensity scores.")
        w <- matrix(1,nrow=nrow(x),ncol=ncol(x))
        i <- data[,treat.var]==0
        x <- data.matrix(x)
         	
        if (estimand=="ATT") {
         	w[i,] <- x[i,]/(1-x[i,])
         	p.s <- x
         }
         
         if (estimand=="ATE") {
         	w[i,] <- 1/(1-x[i,])
         	w[!i,] <- 1/x[!i,]
         	p.s <- x
         }   
      }
      
      if(any(is.infinite(w))) stop("Some propensity weights are infinite.")
      # add a column for unweighted analysis
   } 
   else{
      if(!is.null(sampw)) warning("Sampling weights given when x is a ps object. The sampling weights should be utilized when running ps and are probably not needed here as well and the results may not be correct.")
      # extract the propensity scores and weights from the ps object
      if(estimand != x$estimand) warning("A different estimand was specified when fitting this ps object.  Results for
      the original estimand are being returned.")
      p.s  <- x$ps
      w    <- x$w
      desc <- x$desc
      data <- x$data
      treat.var <- x$treat.var
      estimand <- x$estimand
   }
   if(!all(w[,1]==1)) {
      w   <- cbind(unw=rep(1,nrow(w)),w)
      p.s <- cbind(unw=rep(0.5,nrow(p.s)),p.s)
   }
   if(!is.null(sampw)) w <- w*sampw
   if(is.null(vars)) vars <- names(data)[names(data) != treat.var]

#   summary.tab <- alert <- NULL   
  summary.tab <- NULL
#   zz      <- textConnection("alert","a")
#   if(plots) pdf(file=paste(title,".pdf",sep=""))

   n.tp <- ifelse(class(x)[1]=="ps",length(x$desc),ncol(w))
   if(class(x)[1]!="ps"){ 
     desc<-vector("list",ncol(w))
     names(desc) <- colnames(w)
   }

   for(i.tp in 1:n.tp){
      if((class(x)[1]=="ps") & is.null(vars)){
         desc.temp <- x$desc[[i.tp]]
         iter      <- desc.temp$n.trees
         tp        <- names(x$desc)[i.tp]
      } 
      else {
         desc.temp <- desc.wts(data,
                               w=w[,i.tp],
                               sampw = rep(1,nrow(w)),
                               vars=vars,
                               treat.var=treat.var,
                               perm.test.iters=perm.test.iters,
                               verbose=TRUE,
                               estimand=estimand)
         iter <- NA
         tp <- colnames(w)[i.tp]
         desc[[i.tp]] <- desc.temp
      }
      if(is.null(tp)) tp <- "unnamed"

      summary.tab <- rbind(summary.tab,
         with(desc.temp, data.frame(type      = tp,
                                    n.treat   = n.treat,
                                    n.ctrl    = n.ctrl,
                                    ess.treat = ess.treat,
                                    ess.ctrl  = ess.ctrl,
                                    max.es    = max.es,
                                    mean.es   = mean.es,
                                    max.ks    = max.ks,
                                    mean.ks   = mean.ks,
                                    iter      = iter)))

   } 

#   if(plots) dev.off()

#   cat(alert,sep="\n")
   rownames(summary.tab) <- 1:nrow(summary.tab)

   result <- list(treat      = data[, treat.var],
                  desc       = desc,
                  summary.tab = summary.tab,
                  ps         = as.matrix(p.s[,-1]),
                  w          = as.matrix(w[,-1]),
                  datestamp  = date(),
                  parameters = match.call(),
#                  alerts     = alert,
                  varNames   = vars)
   class(result) <- "dxwts"
   return(result)
}
