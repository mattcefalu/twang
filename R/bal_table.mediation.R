#' Compute the balance table for mediation object.
#'
#' @param x A `mediation` object
#' @param ... list, optional
#'   Additional arguments.
#'
#' @method bal_table mediation
#' @export
bal_table.mediation <- function(x, ...) {

  # we extract what we need from the mediation object
  data <- x$data
  model_a <- x$model_a
  model_m <- x$model_m0
  stopping_methods <- x$stopping_methods
  
  # to get the balance table for NIE, we also need to create a
  # composite weight, which is w_10 for treatment, and w_00
  # for the control group, so we get the weight attributes
  w_00 <- attr(x, 'w_00')
  w_10 <- attr(x, 'w_10')
  
  # then we get rid of the mediation object, since it's so large
  rm(x)
  
  # we find the column names with the `X.` prefix,
  # and use these as the X variables for the model
  column_names <- colnames(data)
  x_names <- column_names[grepl("^X.", column_names)]
  x_and_m_names <- c(x_names, 'M')

  # get the balance table for Model A
  balance_a <- do.call(rbind, bal.table(model_a))
  balance_a['model'] <- 'Model A'
  row.names(balance_a) <- sub(".X.", ".", row.names(balance_a))

  # get the balance table for Model M
  balance_m <- do.call(rbind, bal.table(model_m))
  balance_m['model'] <- 'Model M'
  row.names(balance_m) <- sub(".X.", ".", row.names(balance_m))
  
  nie_wts <- data.frame(matrix(nrow=nrow(w_00), ncol=ncol(w_00)))
  for (i in 1:ncol(w_00)) {nie_wts[,i] <- ifelse(!is.na(w_00[,i]), w_00[,i], w_10[,i])}
  names(nie_wts) <- paste(stopping_methods, 'ATT', sep='.')

  # get the balance table for NIE
  balance_nie <- bal.table(dx.wts(nie_wts,
                                  data = data,
                                  vars = x_and_m_names,
                                  treat.var = 'A',
                                  x.as.weights = TRUE,
                                  estimand = 'ATT'))
  
  balance_nie <- do.call(rbind, balance_nie)
  balance_nie['model'] = 'NIE'

  # keep only the mediator rows that end with '.M'
  possible_prefixes <- paste(stopping_methods, 'ATT', 'M', sep = '.')
  possible_prefixes <- paste(possible_prefixes, collapse = '|')
  rows_to_keep <- grep(possible_prefixes, rownames(balance_nie))
  balance_nie <- balance_nie[rows_to_keep,]
  
  return(list(balance_a=balance_a, balance_m=balance_m, balance_nie=balance_nie))
  
}
