#' Produces a summary table for `ps` object
#'
#' @param x An `ps` object
#' @param ... Additional arguments.
#' 
#' @method print summary.ps
#' @export
#' @md
print.summary.ps <- function(x, ...)
{
      dots <- list(...)
      if(!is.null(dots$digits))
      obj <- round(x, digits = digits)
      else
      obj <- x
      class(obj) <- "matrix"
      print(obj)
      }

