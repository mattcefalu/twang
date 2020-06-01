bal.table.mediation <- function(x, ...) {

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
  x_and_m_names <- c(x$covariate_names, x$mediator_names)

  # get the balance table for Model A
  balance_a <- do.call(rbind, bal.table(model_a))
  balance_a['model'] <- 'Model A'

  # get the balance table for Model M0
  balance_m0 <- do.call(rbind, bal.table(model_m0))
  balance_m0["model"] <- "Model M0"
  
  # get the balance table for Model M1
  balance_m1 <- do.call(rbind, bal.table(model_m1))
  balance_m1["model"] <- "Model M1"
  
  # Compare the weighted distribution of M(1)=m weighted to have the distribution of M(0)=m on the entire pop'n 
  nie_1_wts <- data.frame(matrix(nrow = nrow(w_00), ncol = ncol(w_00)))
  for (i in 1:ncol(w_00)) {
       nie_1_wts[, i] <- ifelse(!is.na(w_00[, i]), w_00[, i], w_10[, i])
   }
  names(nie_1_wts) <- paste(stopping_methods, "ATT", sep = ".")
  balance_nie_1 <- bal.table(dx.wts(nie_1_wts, data = data, 
        vars = x$mediator_names, treat.var = x$a_treatment, x.as.weights = TRUE, 
        estimand = "ATT"))
  balance_nie_1 <- do.call(rbind, balance_nie_1)
  balance_nie_1["model"] = "NIE_1"
  names(balance_nie_1) <- gsub("tx","cntfact",gsub("ct","target",names(balance_nie_1)))
   
  # Compare the weighted distribution of M(0)=m weighted to have the distribution of M(1)=m on the entire pop'n 
  nie_0_wts <- data.frame(matrix(nrow = nrow(w_11), ncol = ncol(w_11)))
  for (i in 1:ncol(w_11)) {
        nie_0_wts[, i] <- ifelse(!is.na(w_11[, i]), w_11[, i], w_01[, i])
  }
  names(nie_0_wts) <- paste(stopping_methods, "ATT", sep = ".")

  balance_nie_0 <- bal.table(dx.wts(nie_0_wts, data = data, 
      vars = x$mediator_names, treat.var = x$a_treatment, x.as.weights = TRUE, 
      estimand = "ATT"))
  balance_nie_0 <- do.call(rbind, balance_nie_0)
  balance_nie_0["model"] = "NIE_0"
  names(balance_nie_0) <- gsub("tx","target",gsub("ct","cntfact",names(balance_nie_0)))
  balance_nie_0 <- balance_nie_0[,c(3,4,1,2,5:10)]
  balance_nie_0$std.eff.sz <- -1*balance_nie_0$std.eff.sz
  balance_nie_0$stat <- -1*balance_nie_0$stat
  return(list(balance_a = balance_a, balance_m0 = balance_m0, balance_m1 = balance_m1, 
        check_counterfactorial_nie_1 = balance_nie_1, check_counterfactorial_nie_0 = balance_nie_0))
  
}
