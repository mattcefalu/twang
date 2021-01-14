#' Compute the balance table for mediation object.
#'
#' @param x A `mediation` object 
#' @param ... Additional arguments.
#' @method bal.table mediation
#' @export
bal.table.mediation <- 
function(x,digits=3, ...) 
{

  # we extract what we need from the mediation object
  data <- x$data
  model_a <- x$model_a
  model_m0 <- x$model_m0
  model_m1 <- x$model_m1
  stopping_methods <- x$stopping_methods
  w_00 <- attr(x, "w_00")
  w_10 <- attr(x, "w_10")
  w_01 <- attr(x, "w_01")
  w_11 <- attr(x, "w_11")
  column_names <- colnames(data)
  m_and_x_names <- c(x$mediator_names,x$covariate_names)

  ##results if method == ps
  if(x$method=="ps") {

    # get the balance table for Model A
    balance_a <- do.call(rbind, bal.table(model_a,digits=digits))
    balance_a['model'] <- 'Model A'

    # get the balance table for Model M0
    balance_m0 <- do.call(rbind, bal.table(model_m0,digits=digits))
    balance_m0["model"] <- "Model M0"
  
    # get the balance table for Model M1
    balance_m1 <- do.call(rbind, bal.table(model_m1,digits=digits))
    balance_m1["model"] <- "Model M1"
  }
  ##results if method = "logistic" or "crossval"
  if(x$method!="ps") {
    # get the balance table for Model A

    if(x$method=="logistic") {
      model_a_preds <- predict(x$model_a,type="response")
    }
   if(x$method=="crossval") {     
      best.iter <- gbm.perf(x$model_a, method="cv",plot.it=FALSE)
      model_a_preds <- predict(x$model_a, n.trees=best.iter, newdata=data, type="response")
    }
    wts_a <- ifelse(data[,x$a_treatment]==1,1/model_a_preds,1/(1-model_a_preds))

    tmp_a <- bal.table(dx.wts(wts_a, data = data, 
        vars = x$covariate_names, treat.var = x$a_treatment, x.as.weights = TRUE, 
        estimand = "ATE"),digits=digits)
    names(tmp_a)[2] <- x$method
    balance_a <- do.call(rbind, tmp_a)
    balance_a['model'] <- 'Model A'

    # get the balance table for Model M0
    data$trt0 <- 1-data[,x$a_treatment]
    if(x$method=="logistic") {
      model_m_preds <- predict(model_m0,type="link")
    }
    else {
      best.iter <- gbm.perf(model_m0, method="cv",plot.it=FALSE)
      model_m_preds <- predict(model_m0, n.trees=best.iter, newdata=data, type="link")
    }
    wts_m0 <- ifelse(data[,x$a_treatment]==0,1,1/exp(model_m_preds))
    tmp_m0 <- bal.table(dx.wts(wts_m0, data = data, 
        vars = m_and_x_names, treat.var = "trt0", x.as.weights = TRUE, 
        estimand = "ATT"),digits=digits)
    names(tmp_m0)[2] <- x$method
    balance_m0 <- do.call(rbind, tmp_m0)
    balance_m0["model"] <- "Model M0"
 
    # get the balance table for Model M1 
    wts_m1 <- ifelse(data[,x$a_treatment]==0,exp(model_m_preds),1)
    tmp_m1 <- bal.table(dx.wts(wts_m1, data = data, 
        vars = m_and_x_names, treat.var = x$a_treatment, x.as.weights = TRUE, 
        estimand = "ATT"),digits=digits)
    names(tmp_m1)[2] <- x$method
    balance_m1<- do.call(rbind, tmp_m1)
    balance_m1["model"] <- "Model M1"
  }
    
  
  # Compare the weighted distribution of M(1)=m weighted to have the distribution of M(0)=m on the entire pop'n 
  # Note: Using ATE here (instead of ATT)
  nie_1_wts <- data.frame(matrix(nrow = nrow(w_00), ncol = ncol(w_00)))
  for (i in 1:ncol(w_00)) {
       nie_1_wts[, i] <- ifelse(!is.na(w_00[, i]), w_00[, i], w_10[, i])
   }
  names(nie_1_wts) <- paste(stopping_methods, "ATE", sep = ".")
  balance_nie_1 <- bal.table(dx.wts(nie_1_wts, data = data, 
        vars = x$mediator_names, treat.var = x$a_treatment, x.as.weights = TRUE, 
        estimand = "ATE"),digits=digits)
  balance_nie_1 <- do.call(rbind, balance_nie_1)
  balance_nie_1["model"] = "NIE_1"
  names(balance_nie_1) <- gsub("tx","cntfact",gsub("ct","target",names(balance_nie_1)))
   
  # Compare the weighted distribution of M(0)=m weighted to have the distribution of M(1)=m on the entire pop'n 
  nie_0_wts <- data.frame(matrix(nrow = nrow(w_11), ncol = ncol(w_11)))
  for (i in 1:ncol(w_11)) {
        nie_0_wts[, i] <- ifelse(!is.na(w_11[, i]), w_11[, i], w_01[, i])
  }
  names(nie_0_wts) <- paste(stopping_methods, "ATE", sep = ".")

  balance_nie_0 <- bal.table(dx.wts(nie_0_wts, data = data, 
      vars = x$mediator_names, treat.var = x$a_treatment, x.as.weights = TRUE, 
      estimand = "ATE"),digits=digits)
  balance_nie_0 <- do.call(rbind, balance_nie_0)
  balance_nie_0["model"] = "NIE_0"
  names(balance_nie_0) <- gsub("tx","target",gsub("ct","cntfact",names(balance_nie_0)))
  balance_nie_0 <- balance_nie_0[,c(3,4,1,2,5:10)]
  balance_nie_0$std.eff.sz <- -1*balance_nie_0$std.eff.sz
  balance_nie_0$stat <- -1*balance_nie_0$stat

  # Remove ks.pvalue from all tables 
  balance_a <- balance_a[,-which(colnames(balance_a)=="ks.pval")]
  balance_m0 <- balance_m0[,-which(colnames(balance_m0)=="ks.pval")]
  balance_m1 <- balance_m1[,-which(colnames(balance_m1)=="ks.pval")]
  balance_nie_1 <- balance_nie_1[,-which(colnames(balance_nie_1)=="ks.pval")]
  balance_nie_0 <- balance_nie_0[,-which(colnames(balance_nie_0)=="ks.pval")]

  # Remove model from all tables 
  balance_a <- balance_a[,-which(colnames(balance_a)=="model")]
  balance_m0 <- balance_m0[,-which(colnames(balance_m0)=="model")]
  balance_m1 <- balance_m1[,-which(colnames(balance_m1)=="model")]
  balance_nie_1 <- balance_nie_1[,-which(colnames(balance_nie_1)=="model")]
  balance_nie_0 <- balance_nie_0[,-which(colnames(balance_nie_0)=="model")]

  cat("**********************************************************\nNotes: Treatment and control are switched for model m0.\nModel m0 is used for NDE_0 and NIE_1 effects.\nModel m1 is used for NDE_1 and NIE_0 effects.\n**********************************************************\n")

  return(list(balance_a = balance_a, balance_m0 = balance_m0,balance_m1 = balance_m1, 
        check_counterfactual_nie_1 = balance_nie_1, check_counterfactual_nie_0 = balance_nie_0))
  
}
