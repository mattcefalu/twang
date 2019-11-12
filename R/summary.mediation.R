#' Displays a useful description of a `mediation` object.
#'
#' @param object mediation object to plot
#' @param ... Additional arguments.
#'
#' @method summary mediation
#' @export
summary.mediation <- function(object, ...) {

  # grab the balance table
  balance_tables <- bal_table.mediation(object)

  # Grab the effects, if they exist
  effects_logical <- grepl('_effects', names(object))
  if (any(effects_logical)) {
    results_table <- object[effects_logical]

  } else {
    results_table <- NULL
  }

  return(list('results_table' = results_table, 'balance_tables' = balance_tables))

}