#' Plot `mniptw`
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
#' @param pairwiseMax If `FALSE`, the plots for the underlying `ps` fits 
#'   will be returned.  Otherwise, pairwise maxima will be returned.
#' @param figureRows The figure rows, passed to [displayPlots]. Default: `NULL`.
#' @param color If `color = FALSE`, figures will be gray scale. Default: `TRUE`.
#' @param subset Used to restrict which of the `stop.method`s will be used 
#'   in the figure. For example `subset = c(1,3)` would indicate that the 
#'   first and third `stop.method`s (in alphabetical order of those specified 
#'   in the original call to `iptw`) should be included in the figure.
#' @param treatments Only applicable when `pairwiseMax` is `FALSE` and `plots` 3, 4, and 5.  
#'   If left at `NULL`, panels for all treatment pairs are created.  If one level of the treatment 
#'   variable is specified, plots comparing that treatment to all others are produced.  If two
#'   levels are specified, a comparison for that single pair is produced.
#' @param singlePlot For Plot calls that produce multiple plots, specifying an integer value of 
#'   `singlePlot` will return only the corresponding plot.  E.g., specifying `singlePlot = 2` 
#'   will return the second plot.
#' @param multiPage When multiple frames of a figure are produced, `multiPage = TRUE` will print each
#'   frame on a different page. This is intended for situations where the graphical output is being
#'   saved to a file. Default: `FALSE`.
#' @param timePeriods The number of distinct time points. If `NULL`, this is assumed to be the number
#'   of `ps` objects (i.e., the number of distinct time points).
#' @param hline Arguments passed to `panel.abline`.
#' @param ... Additional arguments.
#'
#' @method plot mniptw
#' @export
#' @md
plot.mniptw <- function(x,
                        plots="optimize",
                        pairwiseMax = TRUE,
                        figureRows = NULL,
                        color = TRUE,
                        subset = NULL,
                        treatments = NULL,
                        singlePlot = NULL,
                        multiPage = FALSE,
                        timePeriods = NULL,
                        hline=c(0.1,0.5,0.8), ...)
{
   # Creates diag.plot plots and sends to current device
   # x:     ps object 
   # label: Label added to the plot titles

   # extract the propensity scores and weights from the ps object
   
#   singlePlot <- NULL
   
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
   		hdPt[[i]] <- plot(x$psList[[timePeriods[i]]], plots = plots, subset = subset, color = color, time = timePeriods[i], print = FALSE, hline=hline, ...)
   	}
   
   	#if(length(timePeriods) == 1) return(hdPt[[1]])
	#else 
	if(class(hdPt[[1]])[1] == "trellis") {
		if(length(hdPt) == 1) print(hdPt[[1]])
		else displayPlots(hdPt, figureRows = figureRows, singlePlot = singlePlot, multiPage = multiPage)
		}	
	else{
		if(length(timePeriods) > 1) warning("Only returning first time point specified.")
		displayPlots(hdPt[[1]], figureRows = figureRows, singlePlot = singlePlot, multiPage = multiPage, bxpt = plots == 2)
		}
		
}
