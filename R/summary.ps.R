#' Summarize a `ps` object
#'
#' Computes summary information about a stored `ps` object
#'
#' Compresses the information in the `desc` component of the `ps` object
#' into a short summary table describing the size of the dataset and the quality of
#' the propensity score weights.
#'
#' @param object An `ps` object.
#' @param ... Additional arguments.
#'
#' @return See [ps] for details on the returned table.
#'
#' @seealso [ps]
#' @keywords models
#'
#' @method summary ps
#' @export
summary.ps <- function(object,...){
      summary.tab <- NULL
      typ <- NULL   
   
      n.tp <- length(object$desc)
      for(i.tp in 1:n.tp){
         desc.temp <- object$desc[[i.tp]]
         iter      <- desc.temp$n.trees
         tp        <- names(object$desc)[i.tp]

		summary.tab <- rbind(summary.tab,
            with(desc.temp, c(n.treat,n.ctrl,ess.treat,
                                       ess.ctrl,
                                       max.es,
                                       mean.es,
                                       max.ks,
                                       max.ks.p,
                                       mean.ks,
                                       iter)))

                                       
                                       typ <- c(typ, tp)
      }
      

summary.tab <- matrix(summary.tab, nrow = n.tp)
      rownames(summary.tab) <- typ
      colnames(summary.tab) <- c("n.treat", "n.ctrl", "ess.treat", "ess.ctrl", "max.es", "mean.es", "max.ks", "max.ks.p","mean.ks","iter")
      
       
      class(summary.tab) <- "summary.ps"
      return(summary.tab)
}