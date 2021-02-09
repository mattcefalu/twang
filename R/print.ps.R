#' Default print statement for `ps` class
#'
#' @param x An `ps` object
#' @param ... Additional arguments.
#' 
#' @method print ps
#' @export
print.ps <- function(x, ...)
{
   print(summary(x,...))
}

