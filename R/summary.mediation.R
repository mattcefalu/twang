#' Displays a useful description of a `mediation` object.
#'
#' @param object mediation object to plot
#' @param ... Additional arguments.
#'
#' @method summary mediation
#' @export
summary.mediation <- function(object, ...) {

  column_order <- c('model', 'tx.mn', 'tx.sd', 'ct.mn', 'ct.sd',
                    'std.eff.sz', 'stat', 'p', 'ks', 'ks.pval')

  # grab the balance table, reorder the columns, and split into
  # multiple data frames using the `model` column
  balance_tables <- bal_table(object)
  balance_tables <- balance_tables[,column_order]
  balance_tables <- split(balance_tables, balance_tables$model)

  # loop through the effects and put them
  # into a single data frame, then print
  effects_names <- sapply(names(object), function (x) grepl('_effects', x))
  if (any(effects_names)) {

    effects <- object[effects_names]
    results <- matrix(NA, nrow = 7, ncol = length(effects))
    for (effect_idx in 1:length(effects)) {

      if (effect_idx == 1) {
        effects_columns <- names(effects[[effect_idx]])
      }

      for (idx in 1:7) {
        results[idx, effect_idx] <- effects[[effect_idx]][[idx]]
      }

    }
    results_table <- as.data.frame(results)
    colnames(results_table) <- names(effects)
    rownames(results_table) <- effects_columns

    return(list('results_table' = results_table, 'balance_tables' = balance_tables))

  }

}