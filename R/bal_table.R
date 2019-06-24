#' Generic for balance table methods
#'
#' @export
bal_table <- function (x, ...) {
  UseMethod("bal_table", x)
}
