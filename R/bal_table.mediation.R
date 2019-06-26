#' Compute the balance table for mediation object.
#'
#' @name bal_table.mediation
#' @rdname bal_table.mediation
#' 
#' @param x A `mediation` object
#' @param ... list, optional
#'   Additional arguments.
#'
#' @export
bal_table.mediation <- function(x, ...) {

  # TODO : We need a separate impelmentation for outcome mediation
  if (class(x) != 'weighted') {
    stop('Balance table only implemented for weighted mediation objects.')
  }

  fix_med_balance_names <- function(obj) {

    obj <- obj[, order(names(obj))]
    obj_names <- names(obj)

    if ('pval' %in% obj_names) {
      obj_names[which(obj_names == 'pval')] <- 'p'
      names(obj) <- obj_names
    }
    return(obj)
  }

  # if we don't actually have Model A, we need to calculate the
  # balance table using `dx.wts()` instead
  if (is.null(x$model_a)) {

    # we find the column names with the `X.` prefix,
    # and use these as the X variables for the model
    column_names <- colnames(x$data)
    x_names <- column_names[grepl("^X.", column_names)]

    # use `dx.wts()` to extract the balance table
    balance_a <- bal.table(dx.wts(x$model_a_weights,
                                  data = x$data,
                                  vars = x_names,
                                  treat.var = 'A',
                                  x.as.weights = TRUE,
                                  estimand = 'ATE'))
   } else {
    balance_a <- bal.table(x$model_a)
   }
  balance_a <- do.call(rbind, balance_a)
  balance_a['model'] <- 'Model A'

  # get the balance table for Model M
  balance_m <- bal.table(x$model_m)
  balance_m <- do.call(rbind, balance_m)
  balance_m['model'] <- 'Model M'

  # use `dx.wts()` to extract the balance table
  balance_nde <- bal.table(dx.wts(x$natural_direct_weights,
                                  data = x$data,
                                  vars = 'M',
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
