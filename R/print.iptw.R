#' Default print statement for `iptw` class
#' 
#' @param x A `iptw` object
#' @param ... Additional arguments.
#'
#' @method print iptw
#' @export
print.iptw <- function(x, ...)
{
   print(summary(x,...))
}

