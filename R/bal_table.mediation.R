#' Compute the balance table for mediation object.
#'
#' @param x A `mediation` object
#' @param ... list, optional
#'   Additional arguments.
#'
#' @method bal_table mediation
#' @export
bal_table.mediation <- function(x, ...) {

  # we find the column names with the `X.` prefix,
  # and use these as the X variables for the model
  column_names <- colnames(x$data)
  x_names <- column_names[grepl("^X.", column_names)]
  x_and_m_names <- c(x_names, 'M')

  # get the balance table for Model A
  balance_a <- do.call(rbind, bal.table(x$model_a))
  balance_a['model'] <- 'Model A'
  row.names(balance_a) <- gsub(".X.", ".", row.names(balance_a))

  # get the balance table for Model M
  balance_m <- do.call(rbind, bal.table(x$model_m0))
  balance_m['model'] <- 'Model M'
  row.names(balance_m) <- gsub(".X.", ".", row.names(balance_m))
  
  # to get the balance table for NIE, we need to create a
  # composite weight, which is w_10 for treatment, and w_00
  # for the control group
  nie_wts <- data.frame(matrix(nrow=nrow(x$w_00), ncol=ncol(x$w_00)))
  for (i in 1:ncol(x$w_00)) {nie_wts[,i] <- ifelse(!is.na(x$w_00[,i]), x$w_00[,i], x$w_10[,i])}
  names(nie_wts) <- paste(x$stopping_methods, 'ATT', sep='.')

  # get the balance table for NIE
  balance_nie <- bal.table(dx.wts(nie_wts,
                                  data = x$data,
                                  vars = x_and_m_names,
                                  treat.var = 'A',
                                  x.as.weights = TRUE,
                                  estimand = 'ATT'))
  
  balance_nie <- do.call(rbind, balance_nie)
  balance_nie['model'] = 'NIE'

  # keep only the mediator rows that end with '.M'
  rows_to_keep <- endsWith(row.names(balance_nie), '.M')
  balance_nie <- balance_nie[rows_to_keep,]
  
  return(list(balance_a=balance_a, balance_m=balance_m, balance_nie=balance_nie))
  
}
