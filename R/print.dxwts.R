#' Default print statement for `dxwts` class
#' 
#' @param x A `dxwts` object
#' @param ... Additional arguments.
#'
#' @method print dxwts
#' @export
#' @md
print.dxwts <- function(x,...)
{
   print(x$summary.tab)
}