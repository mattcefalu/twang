#' Extract unstabilized propensity score weights for `iptw` fits
#'
#' Extracts propensity score weights from an `iptw` or `mniptw`  object.
#'
#' Weights are the reciprocal of the product of the probability of receiving
#' the treatment received.
#'
#' @param x An `iptw` or `mniptw` object.
#' @param stop.method The twop method used for the fit of interest.
#' @param withSampW Returns weights with sample weights multiplied in, if they were 
#'   provided in the original `iptw` call. Default: `TRUE`.
#'
#' @return Returns a data.frame of weights.
#'
#' @seealso [iptw]
#'
#' @export
get.weights.unstab <- function(x, stop.method = NULL, withSampW = TRUE){
	if(!((class(x) == "iptw") | class(x) == "mniptw")) stop("\"get.weights.unstab\" is only defined for iptw objects.")
	if(class(x) == "iptw"){
	prodWt <- rep(1, length(x$psList[[1]]$treat))
	for(i in 1:length(x$psList)){
		hdWt <- x$psList[[i]]$ps * x$psList[[i]]$treat + (1-x$psList[[i]]$ps) * (1-x$psList[[i]]$treat)
		prodWt <- hdWt * prodWt
	}	
	return(1/prodWt)
	}
	
	if(class(x) == "mniptw"){
		denom <- 1/get.weights(x$psList[[1]], stop.method = stop.method, estimand = "ATE", withSampW = FALSE)
		for(i in 2:length(x$psList)) denom <- denom * 1/get.weights(x$psList[[i]], stop.method = stop.method, estimand = "ATE", withSampW = FALSE)
		w <- 1/denom
		if(withSampW) w <- w * x$psList[[1]]$sampw
	}
}