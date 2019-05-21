#' Plots for `iptw` objects
#'
#' This function produces a collection of diagnostic plots for `iptw` objects.
#'
#' This function produces lattice-style graphics of diagnostic plots.
#'
#' @param x An `iptw` object.
#' @param plots An indicator of which type of plot is desired. The options are
#' * `"optimize" or 1` A plot of the balance criteria as a function of the GBM 
#'   iteration.
#' * `"boxplot" or 2` Boxplots of the propensity scores for the treatment and 
#'   control cases
#' * `"es" or 3` Plots of the standardized effect size of the pre-treatment 
#'   variables before and after reweighing
#' * `"t" or 4` Plots of the p-values from t-statistics comparing means of 
#'   treated and control subjects for pretreatment variables, before and after 
#'   weighting.
#' * `"ks" or 5` Plots of the p-values from Kolmogorov-Smirnov statistics 
#'   comparing distributions of pretreatment variables of treated and control 
#'   subjects, before and after weighting.
#' @param subset Used to restrict which of the `stop.method`s will be used 
#'   in the figure. For example `subset = c(1,3)` would indicate that the 
#'   first and third `stop.method`s (in alphabetical order of those specified 
#'   in the original call to `iptw`) should be included in the figure.
#' @param color If `color = FALSE`, figures will be gray scale. Default: `TRUE`.
#' @param timePeriods The number of distinct time points. If `NULL`, this is assumed to be the number
#'   of `ps` objects (i.e., the number of distinct time points).
#' @param multiPage When multiple frames of a figure are produced, `multiPage = TRUE` will print each
#'   frame on a different page. This is intended for situations where the graphical output is being
#'   saved to a file. Default: `FALSE`.
#' @param figureRows The figure rows, passed to [displayPlots]. Default: `NULL`.
#' @param hline Arguments passed to `panel.abline`.
#' @param ... Additional arguments.
#'
#' @references Dan McCaffrey, G. Ridgeway, Andrew Morral (2004). "Propensity
#'   Score Estimation with Boosted Regression for Evaluating Adolescent
#'   Substance Abuse Treatment", *Psychological Methods* 9(4):403-425.
#'
#' @seealso [iptw]
#' @keywords models, multivariate
#'
#' @method plot iptw
plot.iptw <- function(x,
                      plots = "optimize",
                      subset = NULL,
                      color = TRUE,
                      timePeriods = NULL,
                      multiPage = FALSE,
                      figureRows = NULL,
                      hline=c(0.1,0.5,0.8),
                      ...)
{
   # Creates diag.plot plots and sends to current device
   # x:     ps object 
   # label: Label added to the plot titles

   # extract the propensity scores and weights from the ps object
   
   singlePlot <- NULL
   
	if(!is.numeric(subset) & !is.null(subset)){
		if(!all(subset %in% x$stopMethods)) stop("The \"subset\" arugment must be NULL, numeric, or one of the stop.methods specified when fitting the mnps object.")
	}
   	if(is.null(subset)) subset <- 1:length(x$stopMethods)
   	if(!is.numeric(subset)){
   		hldLen <- 1:length(x$stopMethods)
   		subset <- hldLen[x$stopMethods %in% subset]
   	}
   	
   	subset <- sort(subset)
   	
   	if(is.null(subset)) stop("The \"subset\" arugment must either be NULL, numeric, or some subset of", print(x$stopMethods))    
   
   if(is.null(timePeriods)) timePeriods <- 1:x$nFits
   
   hdPt <- vector(mode = "list", length = length(timePeriods))
   
   for(i in 1:length(timePeriods)){
   	hdPt[[i]] <- diag.plot.color(x$psList[[i]], plots, subset = subset, color = color, time = timePeriods[i], hline=hline,...)
   }
   
   	if(length(timePeriods) == 1) return(hdPt[[1]])
	else displayPlots(hdPt, figureRows = figureRows, singlePlot = singlePlot, multiPage = multiPage)	
}
