#' Default print statement for `mniptw` class
#'
#' @param x A `mniptw` object
#' @param ... Additional arguments.
#'
#' @method print mniptw
#' @export
print.mniptw <- function(x, ...)
{
   print(summary(x,...))
}

