#' Extract propensity score weights.
#'
#' Extracts propensity score weights from a ps or mnps object.
#'
#' Weights for ATT are 1 for the treatment cases and p/(1-p) for the control cases. 
#' Weights for ATE are 1/p for the treatment cases and 1/(1-p) for the control cases.
#'
#' @param ps1 A `ps` or `mnps` object.
#' @param stop.method Indicates which set of weights to retrieve from the `ps` object.
#' @param estimand Indicates whether the weights are for the average treatment effect on
#'   the treated (ATT) or the average treatment effect on the population (ATE). By default,
#'   `get.weights` will use the estimand used to fit the `ps` object.
#' @param withSampw Whether to return weights with sample weights multiplied in, if they were
#'   provided in the original `ps` or `mnps` call. Default: `TRUE`.
#'
#' @return Returns a vector of weights.
#'
#' @seealso [ps]
#'
#' @export
get.weights <- function(ps1, stop.method = NULL, estimand = NULL, withSampW = TRUE)
{
   if(class(ps1)=="ps"){
   if(is.null(estimand)) estimand <- ps1$estimand
   
   if(!(estimand %in% c("ATT","ATE"))) stop("estimand must be either \"ATT\" or \"ATE\".")
   if(estimand != ps1$estimand){
   	warning("Estimand specified for get.weights() differs from the estimand used to fit the ps object.")
   }
   if(length(stop.method)>1) stop("More than one stop.method was selected.")
   if(!is.null(stop.method))
   {
   	stop.method.long <- paste(stop.method, ps1$estimand, sep=".")
      i <- match(stop.method.long, names(ps1$w))
      if(is.na(i)) stop("Weights for stop.method=",stop.method, " and estimand=", estimand, " are not available. Please a stop.method and used when fitting the ps object.") 
#      Available options: ",names(ps1$ps),".")
   } else
   {
   	warning("No stop.method specified.  Using ", names(ps1$ps)[1], "\n")
      i <- 1
   }
  
   if (estimand == "ATT") {
   	w <- with(ps1,  treat + (1- treat) * ps[[i]]/(1-ps[[i]]))
   	if(withSampW) w <- w * ps1$sampw
   	return(w)
   }
   else if (estimand == "ATE")
   { 
      w <- with(ps1, treat/ps[[i]] + (1-treat)/(1-ps[[i]]))
      if(withSampW) w <- w* ps1$sampw
      return(w)
   }
   }
   if(class(ps1) == "mnps"){
   if(is.null(estimand)) estimand <- ps1$estimand
   
   if(!(estimand %in% c("ATT","ATE"))) stop("estimand must be either \"ATT\" or \"ATE\".")
   if(estimand != ps1$estimand){
   	warning("Estimand specified for get.weights() differs from the estimand used to fit the ps object.")
   }
   if(length(stop.method)>1) stop("More than one stop.method was selected.")
   if(!is.null(stop.method))
   {
   	stop.method.long <- paste(stop.method, ps1$estimand, sep=".")
      i <- match(stop.method.long, names(ps1$psList[[1]]$ps))
      if(is.na(i)) stop("Weights for stop.method=",stop.method, " and estimand=", estimand, " are not available. Please a stop.method and used when fitting the mnps object.") 
#      Available options: ",names(ps1$ps),".")
   } else
   {
   	if(length(ps1$stopMethods) > 1) warning("No stop.method specified.  Using ", names(ps1$psList[[1]]$ps)[1], "\n")
      i <- 1
   }
  
   if (estimand == "ATT") {
   	w <- rep(0, nrow(ps1$data))
   	w[ps1$treatVar == ps1$treatATT] <- 1
   	for(j in 1:length(ps1$levExceptTreatATT)){
   		### in mnps, treatATT is the treatment in the individual fits
   		### fill in weights for categories other than treatATT first, using (1-ps)/ps
   		w[ps1$treatVar %in% c(ps1$levExceptTreatATT[j], ps1$treatATT)] <- with(ps1$psList[[j]],  (1-treat) * ps[[i]]/(1-ps[[i]]))
   		}
   		w[ps1$treatVar == ps1$treatATT] <- 1
   		if(withSampW) w <- w * ps1$sampw
   		return(w)
   }
   else if (estimand == "ATE"){
    	w <- with(ps1$psList[[1]],  treat/ps[[i]])
	   	for(j in 2:length(ps1$psList)){
	   		w <- w + with(ps1$psList[[j]],  treat/ps[[i]])
   		}
   		if(withSampW) w <- w * ps1$sampw
   		return(w)
   		}
  	 }	 
#      w <- with(ps1, treat/ps[[i]] + (1-treat)/(1-ps[[i]]))
#      if(withSampW) w <- w* ps1$sampW
#      return(w)
#   }
   	
   
   if(!(class(ps1) %in% c('ps', 'mnps'))) stop("The object 'ps1' must be of class 'ps' or 'mnps'.")
}
