#' Boxplot for `ps` objects
#'
#' This function produces a collection of diagnostic plots for ps objects.
#'
#' This function produces lattice-style graphics of diagnostic plots. 
#'
#' @param x A `ps` object
#' @param subset If multiple `stop.method` rules were used in the `ps()` call, `subset` 
#'   restricts the plots of a subset of the stopping rules that were employed. This argument
#'   expects a subset of the integers from 1 to k, if k `stop.method`s were used.
#' @param color If `FALSE`, a grayscale figure will be returned.
#' @param time For use with iptw fits.
#' @param ... Additional arguments that are passed to boxplot function, which may bepassed to
#'   the underlying `lattice` package plotting functions.
#'
#' @seealso [ps]
#' @keywords multivariate
#'
#' @references Dan McCaffrey, G. Ridgeway, Andrew Morral (2004). "Propensity
#'   Score Estimation with Boosted Regression for Evaluating Adolescent
#'   Substance Abuse Treatment", *Psychological Methods* 9(4):403-425.
#'
#' @method boxplot ps
#' @export
boxplot.ps <- function(x,
                       subset = NULL,
                       color = TRUE,
                       time = NULL,
                       ...){
	longDat <- matrix(t(as.matrix(x$ps)), ncol = 1)
	nms <- names(x$ps)
	bwDat <- data.frame(ps = longDat, nm = nms, treat = rep(x$treat, each = length(nms)))
	if(is.null(subset)) subset <- 1:length(levels(as.factor(bwDat$nm)))
	ptSymCol <- ifelse(color, "#0080ff", "black")	
	bwCols <- list(col = ptSymCol)
	stripBgCol <- ifelse(color, "#ffe5cc", "transparent")
	
	if(is.null(time)) xlb = "Propensity scores"
	else xlb = paste("Propensity scores (Time ", time, ")", sep = "")

	pt1 <- bwplot(treat~ps|nm, data=bwDat, scales = list(alternating = 1), ylab = "Treatment", xlab=xlb, subset = as.factor(bwDat$nm) %in% levels(as.factor(bwDat$nm))[subset], par.settings = list(strip.background = list(col=stripBgCol), box.rectangle = bwCols, plot.symbol = bwCols, box.umbrella = bwCols), ...)
	return(pt1)
	
}