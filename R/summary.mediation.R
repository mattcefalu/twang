#' Displays a useful description of a `mediation` object.
#'
#' @param object A `mediation` object 
#' @param ... Additional arguments.
#'
#' @method summary mediation
#' @export
summary.mediation	<- 
function(object,...) 
{
   # Get confidence intervals of effects 
    if(!is.null(object$y_outcome)) {
      desc_effects  <- desc.effects.mediation(object)
    }
    else {
      desc_effects  <- NULL
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

    return(list(results_table = desc_effects, ps_summary_tables = ps_tables, mediator_distribution  = mediator_distribution_check ))
}
