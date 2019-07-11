#' Compute the balance table for mediation object.
#'
#' @param x A `mediation` object
#' @param ... list, optional
#'   Additional arguments.
#'
#' @method bal_table mediation
#' @export
bal_table.mediation <- function(x, ...) {

  fix_med_balance_names <- function(obj) {

    obj <- obj[, order(names(obj))]
    obj_names <- names(obj)

    if ('pval' %in% obj_names) {
      obj_names[which(obj_names == 'pval')] <- 'p'
      names(obj) <- obj_names
    }
    return(obj)
  }

  # we find the column names with the `X.` prefix,
  # and use these as the X variables for the model
  column_names <- colnames(x$data)
  x_names <- column_names[grepl("^X.", column_names)]
  x_and_m_names <- c(x_names, 'M')

  # get the balance table for Model A
  balance_a <- bal.table(dx.wts(x$model_a_wts,
                                data = x$data,
                                vars = x_names,
                                treat.var = 'A',
                                x.as.weights = TRUE,
                                estimand = 'ATE'))
  balance_a <- do.call(rbind, balance_a)
  balance_a['model'] <- 'Model A'

  # get the balance table for Model M
  balance_m <- bal.table(dx.wts(x$model_m_wts,
                                data = x$data,
                                vars = x_and_m_names,
                                treat.var = 'A',
                                x.as.weights = TRUE,
                                estimand = 'ATT'))
  balance_m <- do.call(rbind, balance_m)
  balance_m['model'] <- 'Model M'

  # get the balance table for NDE
  balance_nde <- bal.table(dx.wts(x$natural_direct_wts,
                                  data = x$data,
                                  vars = x_and_m_names,
                                  treat.var = 'A',
                                  x.as.weights = TRUE,
                                  estimand = 'ATT'))
  balance_nde <- do.call(rbind, balance_nde)
  balance_nde['model'] = 'NDE'

  # make sure the balance tables have the same column names
  balance_a <- fix_med_balance_names(balance_a)
  balance_m <- fix_med_balance_names(balance_m)
  balance_nde <- fix_med_balance_names(balance_nde)

  rbind(balance_a, balance_m, balance_nde)
}
