#' Summarize a `mniptw` object
#'
#' @param object A `mniptw` object.
#' @param ... Additional arguments.
#'
#' @method summary mniptw
#' @export
#' @md
summary.mniptw <- function (object, ...){
	nFits <- object$nFits
	summaryList <- vector(mode = "list", length = nFits)
	for (i in 1:nFits) {
		summaryList[[i]] <- summary(object$psList[[i]], ...)
	}

	retObj <- list(summaryList = summaryList, nFit = object$nFit, uniqueTimes = object$uniqueTimes)

	class(retObj) <- "summary.mniptw"
	return(retObj)
}

