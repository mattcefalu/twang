#' Displays a useful description of a `mediation` object.
#'
#' @param object mediation object to plot
#' @param ... Additional arguments.
#'
#' @method summary mediation
#' @export
summary.mediation <- function(object, ...) {
  
  # TODO : We'll definitely want to fix / clean up this function
  
  # collect the descriptives for both models,
  # and put these into a named list, `models`
  model_a <- object$model_a$desc
  model_m <- object$model_m$desc
  models <- list("Model A" = model_a,
                 "Model M" = model_m)
  
  # loop through through each model and each stopping
  # method, and grab the balance table, then print
  balance_tables = list()
  for (model_name in names(models)) {
    model <- models[[model_name]]
    
    for (idx in 1:length(model)) {
      balance_table <- model[[idx]]$bal.tab$results
      balance_tables[[paste(model_name, names(model)[[idx]], sep = ", ")]] <- balance_table
    }
  }
  balance_tables <- do.call(rbind, balance_tables)
  
  # loop through the effects and put them
  # into a single data frame, then print
  effects_boolean <- sapply(names(object), function (x) grepl('_effects', x))
  if (any(effects_boolean)) {
    
    effects <- object[effects_boolean]
    results <- matrix(NA, nrow = 7, ncol = length(effects))
    value_names = NULL
    for (stopping_method_index in 1:length(effects)) {
      
      if (stopping_method_index == 1) {
        value_names <- names(effects[[stopping_method_index]])
      }
      
      for (value_index in 1:7) {
        results[value_index, stopping_method_index] <- effects[[stopping_method_index]][[value_index]]
      }
      
    }
    results_table <- as.data.frame(results)
    colnames(results_table) <- names(effects)
    rownames(results_table) <- value_names
    
    return(list('results_table' = results_table, 'balance_tables' = balance_tables))
    
  }
  
}