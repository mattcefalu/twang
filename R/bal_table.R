#' Generic for balance table methods
#' 
#' @param x An object
#' @param ... list, optional
#'   Additional arguments.
#'
#' @export
bal_table <- function (x, ...) {
  UseMethod("bal_table", x)
}
