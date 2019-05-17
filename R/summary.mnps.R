#' Summarize a `mnps` object
#'
#' Computes summary information about a stored `mnps` object
#'
#' Compresses the information in the `desc` component of the `mnps` object
#' into a short summary table describing the size of the dataset and the quality of
#' the propensity score weights.
#'
#' @param object An `mnps` object.
#' @param ... Additional arguments.
#'
#' @return See [mnps] for details on the returned table.
#'
#' @seealso [mnps]
#' @keywords models
#'
#' @method summary mnps
summary.mnps <- function (object, ...){

	nFits <- object$nFits
	summaryList <- vector(mode = "list", length = nFits)
	for (i in 1:nFits) {
		summaryList[[i]] <- summary(object$psList[[i]], ...)
	}
	ns <- NULL
	for (i in 1:length(summaryList)) ns <- c(ns, summaryList[[i]][1,1])
	if(object$estimand == "ATE"){.txs <- object$treatLev}
	else .txs <- object$levExceptTreatATT
	shell <- data.frame(treatment = .txs, n = ns)
	for (i in 1:length(object$stopMethods)) {
		hldESS <- NULL
		for (j in 1:nrow(shell)) {
			hldESS <- c(hldESS, ifelse(object$estimand == "ATE",
				summaryList[[j]][i + 1, "ess.treat"], summaryList[[j]][i + 1, "ess.ctrl"]))
		}
		shell <- data.frame(shell, currESS = hldESS)
		names(shell)[names(shell) == "currESS"] <- paste("ESS:", object$stopMethods[i], sep = "")
	}
	if (object$estimand == "ATT") {
		retObj <- list(summaryList = summaryList, nFit = object$nFit,
			estimand = object$estimand, treatATT = object$treatATT,
			treatLev = object$treatLev, levExceptTreatATT = object$levExceptTreatATT,
			ess = shell)
	}
	else {
		hd1 <- pairwiseComparison(object, collapse.to = "stop.method" , ks.exact = object$ks.exact)
		retObj <- list(comp = hd1, ess = shell, estimand = object$estimand)
	}
	class(retObj) <- "summary.mnps"
	return(retObj)
}

