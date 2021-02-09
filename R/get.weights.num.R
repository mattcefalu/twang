#' Get numerators to stabilize propensity score weights for `iptw` fits.
#'
#' Forms numerators to stabilize weights for an iptw object.
#'
#' @param iptw An `iptw`` object.
#' @param fitList A list containing objects with an associated "fitted" function.
#'
#' @return Returns numerator of stabilized weights to be used in conjunction with 
#'   `get.weights.unstab`
#'
#' @seealso [iptw]
#'
#' @export
get.weights.num <- function(iptw, fitList){
	hdWt <- rep(1, length(iptw$psList[[1]]$treat))
	for(i in 1:length(iptw$psList)){
		hdPred <- fitted(fitList[[i]]) * iptw$psList[[i]]$treat + (1-fitted(fitList[[i]])) * (1-iptw$psList[[i]]$treat)
		hdWt <- hdWt * hdPred
	}
	return(hdWt)
}