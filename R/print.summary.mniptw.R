#' Produces a summary table for `mniptw` object 
#'
#' @param x An `mniptw` object
#' @param ... Additional arguments.
#'
#' @method print summary.mniptw
#' @export
print.summary.mniptw <- function(x, ...)
{
      	nSum <- length(x$summaryList)      	
      	for(i in 1:nSum){	
      		cat("Summary for time period ", x$uniqueTimes[i], ": \n")
      		print(x$summaryList[[i]])
      		if (i<nSum) cat("\n")
      		
      	}

      invisible(x)
     }

