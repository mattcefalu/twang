#' Default print statement for `mnps` class
#'
#' @param x A `mnps` object
#' @param ... Additional arguments.
#'
#' @method print mnps
#' @export
#' @md
print.mnps <- function(x, ...)
{
   print(summary(x,...))
}

