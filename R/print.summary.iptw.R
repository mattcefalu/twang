#' Produces a summary table for `iptw` object
#'
#' @param x An `iptw` object
#' @param ... Additional arguments.
#'
#' @method print summary.iptw
#' @export
#' @md
print.summary.iptw <- function(x, ...)
{
      	nSum <- length(x$summaryList)      	
      	for(i in 1:nSum){	
      		cat("Summary for time period ", x$uniqueTimes[i], ": \n")
      		print(x$summaryList[[i]])
      		if (i<nSum) cat("\n")
      	}

      invisible(x)
     }

