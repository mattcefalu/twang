#' Summarize a `iptw` object
#'
#' Computes summary information about a stored `iptw` object
#'
#' Compresses the information in the `desc` component of the `iptw` object
#' into a short summary table describing the size of the dataset and the quality of
#' the propensity score weights.
#'
#' @param object An `iptw` object.
#' @param ... Additional arguments.
#'
#' @return See [iptw] for details on the returned table.
#'
#' @seealso [iptw]
#' @keywords models
#'
#' @method summary iptw
summary.iptw <- function (object, ...){
	nFits <- object$nFits
	summaryList <- vector(mode = "list", length = nFits)
	for (i in 1:nFits) {
		summaryList[[i]] <- summary(object$psList[[i]], ...)
	}

	retObj <- list(summaryList = summaryList, nFit = object$nFit, uniqueTimes = object$uniqueTimes)

	class(retObj) <- "summary.iptw"
	return(retObj)
}

