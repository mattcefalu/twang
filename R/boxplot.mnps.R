#' Boxplot for `mnps` objects
#'
#' This function produces a collection of diagnostic plots for mnps objects.
#'
#' This function produces lattice-style graphics of diagnostic plots. 
#'
#' @param x A `ps` object
#' @param stop.method Only 1 `stop.method` can be presented at a time for `mnps` objects. 
#'   Use a numeric indicator of which `stop.method` (among those specified when fitting
#'   the `mnps` object) should be used.
#' @param color If `FALSE`, a grayscale figure will be returned.
#' @param figureRows The number of rows in the figure. Defaults to the number of panels.
#' @param singlePlot If multiple sets of boxplots are produced, `singlePlot` can be used
#'   to select only one. For example, `singlePlot = 2` would return only the second set
#'   of boxplots.
#' @param multiPage When multiple frames of a figure are produced, `multiPage = TRUE` will
#'   print each frame on a different page. This is intended for situations where the graphical
#'   output is being saved to a file.
#' @param time For use with iptw fits.
#' @param print If `FALSE`, the figure is returned but not printed.
#' @param ... Additional arguments that are passed to boxplot function, which may bepassed to
#'   the underlying `lattice` package plotting functions.
#'
#' @seealso [mnps]
#' @keywords multivariate
#'
#' @references Dan McCaffrey, G. Ridgeway, Andrew Morral (2004). "Propensity
#'   Score Estimation with Boosted Regression for Evaluating Adolescent
#'   Substance Abuse Treatment", *Psychological Methods* 9(4):403-425.
#'
#' @method boxplot mnps
#' @export
boxplot.mnps <- function(x,
                         stop.method = NULL,
                         color = TRUE,
                         figureRows = NULL,
                         singlePlot = NULL,
                         multiPage = FALSE,
                         time = NULL,
                         print = TRUE,
                         ...){
	
	ptSymCol <- ifelse(color, "#0080ff", "black")	
	bwCols <- list(col = ptSymCol)
	stripBgCol <- ifelse(color, "#ffe5cc", "transparent")
	nPlot <- x$nFits
	
	
	ptHld <- vector(mode = "list", length = nPlot)
	
	
	if(is.null(stop.method)) stop.method <- x$stopMethods
	
	if(!is.null(singlePlot)){
		if(!is.numeric(singlePlot) | singlePlot > x$nFits) stop(paste("If specified, the \"singlePlot\" argument must be an integer between 1 and", x$nFits, "for this object."))
		if(round(singlePlot) != singlePlot) stop("If specified, the \"singlePlot\" argument must be a positive integer.")
	}
	
	if(length(stop.method) > 1){ 
		if(is.numeric(stop.method))
		warning("Using only the first stop.method, ",  x$stopMethods[1])
		else warning("Using only the first stop.method, ",  stop.method[1])
		stop.method <- stop.method[1]
		}
			
	
	
	if(is.numeric(stop.method)) stop.method = x$stopMethods[stop.method]
	

	
	if(x$estimand == "ATE"){

		bwDat <- whichResp <- NULL
			
		stopMethLong <- paste(stop.method, ".ATE", sep = "")
			
			
		for(j in 1:nPlot){
			bwDat <- data.frame(ps = x$psList[[j]]$ps[,stopMethLong], treat = x$treatVar, whichResp = x$treatLev[j])
			ptNm <- paste(x$treatLev[j], " propensity scores by Tx group")
			if(!is.null(time)) ptNm <- paste(ptNm, " (time ", time, ")", sep = "")
			pt1 <- bwplot(ps ~ treat, groups = whichResp,  
		xlab = "Treatment", ylab = "Propensity scores", ylim = c(-.1,1.1), data = bwDat, main = ptNm, par.settings = list(strip.background = list(col=stripBgCol), box.rectangle = bwCols, plot.symbol = bwCols, box.umbrella = bwCols), ...)
		
		ptHld[[j]] <- pt1

		}
		}
		
		
	else if(x$estimand == "ATT"){
		
		bwDat <- NULL
		
		stopMethLong <- paste(stop.method, ".ATT", sep = "")
		bwDat <- NULL
		
		for(j in 1:nPlot){
			currCats <- c(x$treatATT, x$levExceptTreatATT[j])
			bwDat <- data.frame(ps = x$psList[[j]]$ps[,stopMethLong], treat = currCats[1 + x$psList[[j]]$treat], respCat = x$levExceptTreatATT[j], attGrp = x$treatATT)
			ptNm <- paste("Propensity score of ", x$levExceptTreatATT[j], " versus ", x$treatATT, ".", sep = "")
			if(!is.null(time)) ptNm <- paste(ptNm, " (time ", time, ")", sep = "")
			pt1 <- bwplot(ps ~ treat, data = bwDat, ylim = c(-.1,1.1), ylab = "Propensity scores", xlab = "Treatment", main = ptNm,par.settings = list(strip.background = list(col=stripBgCol), box.rectangle = bwCols, plot.symbol = bwCols, box.umbrella = bwCols), ...)
			
			ptHld[[j]] <- pt1
	
			
		}		
	}	
	
	if(print) displayPlots(ptHld, figureRows = figureRows, singlePlot = singlePlot, multiPage = multiPage, bxpt = TRUE)
	return(ptHld)


}