#' Describe the effects
#' 
#' Describe the effects, and calculate standard
#' errors and confidence intervals
#' 
#' @param x An object
#' @param ... list, optional
#'   Additional arguments.
#'
#' @export
desc_effects <- function (x, ...) {
   UseMethod("desc_effects", x)
}
