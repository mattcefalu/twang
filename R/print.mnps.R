#' Default print statement for `mnps` class
#'
#' @param x An `mnps` object
#' @param ... Additional arguments.
#'
#' @method print mnps
#' @export
print.mnps <- function(x, ...)
{
   print(summary(x,...))
}

