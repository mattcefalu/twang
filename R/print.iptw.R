#' Default print statement for `iptw` class
#' 
#' @param x A `iptw` object
#' @param ... Additional arguments.
#'
#' @method print iptw
#' @export
#' @md
print.iptw <- function(x, ...)
{
   print(summary(x,...))
}

