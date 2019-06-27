#' Displays a useful description of a `mediation` object.
#'
#' @param object mediation object to plot
#' @param ... Additional arguments.
#'
#' @method summary mediation
#' @export
summary.mediation <- function(object, ...) {

  # TODO : We need a separate impelmentation for outcome mediation
  if (object$mediation_type != 'weighted') {
    stop('Summary only implemented for weighted mediation objects.')
  }

  # collect the descriptives for both models,
  # and put these into a named list, `models`
  model_a <- object$model_a$desc
  model_m <- object$model_m$desc
  models <- list("Model A" = model_a,
                 "Model M" = model_m)
  
  balance_tables <- bal_table(object)
  
  # loop through the effects and put them
  # into a single data frame, then print
  effects_names <- sapply(names(object), function (x) grepl('_effects', x))
  if (any(effects_names)) {
    
    effects <- object[effects_names]
    results <- matrix(NA, nrow = 6, ncol = length(effects))
    for (effect_idx in 1:length(effects)) {
      
      if (effect_idx == 1) {
        effects_columns <- names(effects[[effect_idx]])
      }

      for (idx in 1:6) {
        results[idx, effect_idx] <- effects[[effect_idx]][[idx]]
      }
      
    }
    results_table <- as.data.frame(results)
    colnames(results_table) <- names(effects)
    rownames(results_table) <- effects_columns
    
    return(list('results_table' = results_table, 'balance_tables' = balance_tables))
    
  }
  
}