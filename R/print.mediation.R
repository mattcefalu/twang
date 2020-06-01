#' Default print statement for `mediation` class
#'
#' @param object A `mediation` object
#' @param ... Additional arguments.
#' 
#' @method print mediation
#' @export
print.mediation <- function(object, ...)
{
  # Grab the effects, if they exist
  effects_logical <- grepl('_effects', names(object))
  if (any(effects_logical)) {
    estimates_table <- object[effects_logical]
  } else {
    estimates_table <- NULL
  }

  # Get summaries of ps objects
  model_a  <- summary(object$model_a)
  model_m0 <- summary(object$model_m0)
  model_m1 <- summary(object$model_m1)
  ps_tables  <- list(model_a=model_a,model_m0=model_m0,model_m1=model_m1)

  # Get balance tables for NIE_1 and NIE_0
  # to check that weights for the counterfactual 
  # mediator distributions yeild distributions of 
  # mediators that match the target
  mediator_distribution_check <- bal.table.mediation(object)[c("check_counterfactorial_nie_1","check_counterfactorial_nie_0")]

  print(list(estimates_table = estimates_table, ps_summary_tables = ps_tables, mediator_distribution = mediator_distribution_check))
}