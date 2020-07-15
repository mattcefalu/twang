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
    if(object$method=="ps") {
      model_a  <- summary(object$model_a)
      model_m0 <- summary(object$model_m0)
      model_m1 <- summary(object$model_m1)
    } 
    if(object$method!="ps") {
      data <- object$data 

      wts_a <- attr(object,"w_11")
      wts_a[is.na(wts_a)] <- attr(object,"w_00")[is.na(wts_a)]
      dx_a <- dx.wts(wts_a, data = data, 
          vars = object$covariate_names, treat.var = object$a_treatment, x.as.weights = TRUE, 
          estimand = "ATT")      
      dx_a$desc[[1]]["iter"] <- NA
      dx_a$desc[[2]]["iter"] <- NA
      names(dx_a$desc)[2] <- object$method
      model_a <- twang:::summary.ps(dx_a)

      data$trt0 <- 1-data[,object$a_treatment]
      if(object$method=="logistic") {
        model_m_preds <- predict(object$model_m0,type="link")
      }      
        else {
        best.iter <- gbm.perf(object$model_m0, method="cv",plot.it=FALSE)
        model_m_preds <- predict(object$model_m0, n.trees=best.iter, newdata=data, type="link")
      }
      wts_m0 <- ifelse(data[,object$a_treatment]==0,1,1/exp(model_m_preds))
      dx_m0 <- dx.wts(wts_m0, data = data, 
          vars = c(object$mediator_names,object$covariate_names), treat.var = "trt0", x.as.weights = TRUE, 
          estimand = "ATT")
      dx_m0$desc[[1]]["iter"] <- NA
      dx_m0$desc[[2]]["iter"] <- NA
      names(dx_m0$desc)[2] <- object$method
      model_m0 <- twang:::summary.ps(dx_m0)

      wts_m1 <- ifelse(data[,object$a_treatment]==0,exp(model_m_preds),1)
      dx_m1 <- dx.wts(wts_m1, data = data, 
        vars =  c(object$mediator_names,object$covariate_names), treat.var = object$a_treatment, x.as.weights = TRUE, 
        estimand = "ATT")
      dx_m1$desc[[1]]["iter"] <- NA
      dx_m1$desc[[2]]["iter"] <- NA
      names(dx_m1$desc)[2] <- object$method
      model_m1 <- twang:::summary.ps(dx_m1)
    }

    ps_tables  <- list(model_a=model_a,model_m0=model_m0,model_m1=model_m1)

  # Get balance tables for NIE_1 and NIE_0
  # to check that weights for the counterfactual 
  # mediator distributions yeild distributions of 
  # mediators that match the target
  mediator_distribution_check <- bal.table.mediation(object)[c("check_counterfactorial_nie_1","check_counterfactorial_nie_0")]

  print(list(estimates_table = estimates_table, ps_summary_tables = ps_tables, mediator_distribution = mediator_distribution_check))
}