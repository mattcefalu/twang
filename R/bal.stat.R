#' Calculate weighted balance statistics
#'
#' `bal.stat` compares the treatment and control subjects by means, standard 
#' deviations, effect size, and KS statistics
#'
#' `bal.stat` calls  auxiliary functions for each variable and assembles
#' the results in a table.
#'
#' @param data A data frame containing the data
#' @param vars A vector of character strings with the names of the variables 
#'   on which the function will assess the balance
#' @param treat.var The name of the treatment variable
#' @param w.all Oobservation weights (e.g. propensity score weights, sampling 
#'   weights, or both)
#' @param sampw Sampling weights. These are passed in addition to `w.all` because
#'   the "unweighted" results shoud be adjusted for sample weights (though not propensity
#'   score weights).            
#' @param get.means logical. If `TRUE` then `bal.stat` will compute means
#'   and variances
#' @param get.ks logical. If `TRUE` then `bal.stat` will compute KS statistics
#' @param na.action A character string indicating how `bal.stat` should 
#'   handle missing values. Current options are "level", "exclude", or "lowest"
#' @param  estimand Either "ATT" or "ATE"
#' @param multinom logical. `TRUE` if used for multinomial propensity scores.
#' @param fillNAs logical. If `TRUE`, fills in zeros for missing values.
#'
#' @return `get.means` and `get.ks` manipulate the inclusion of certain 
#'   columns in the returned result.
#'
#' @seealso The example for [ps] contains an example of the use of [bal.table]
#' @keywords multivariate
#'
#' @references Dan McCaffrey, G. Ridgeway, Andrew Morral (2004). "Propensity
#'   Score Estimation with Boosted Regression for Evaluating Adolescent
#'   Substance Abuse Treatment", *Psychological Methods* 9(4):403-425.
#'
#' @export
bal.stat <- function(data,
                     vars = NULL,
                     treat.var,
                     w.all,
                     sampw, 
                     get.means = TRUE,
                     get.ks = TRUE,
                     na.action = "level",
                     estimand,
                     multinom,
                     fillNAs = FALSE)
{
   if(is.null(vars)) vars<-names(data)[names(data)!=treat.var]

   is.fac   <- sapply(data[,vars,drop=FALSE],is.factor)
   fac      <- vars[is.fac]
   not.fac  <- vars[!is.fac]

   ret <- vector("list",length(vars))
   names(ret) <- vars
   
#   multinom <- FALSE

sampW <- sampw

   ##### Calculate stats for numeric variables
   ret[!is.fac] <- lapply(data[,vars[!is.fac],drop=FALSE], ps.summary.new2,
                          t=data[,treat.var], w=w.all, sampw = sampW,
                          get.means=get.means, get.ks=get.ks,
                          na.action=na.action,
                          collapse.by.var=FALSE, estimand=estimand, multinom = multinom, fillNAs = fillNAs)

   ##### Calculate stats for factor variables
   ret[is.fac] <- lapply(data[,vars[is.fac],drop=FALSE], ps.summary.new2,
                         t=data[,treat.var], w=w.all, sampw = sampW, 
                         get.means=get.means, get.ks=get.ks,
                         na.action=na.action,
                         collapse.by.var=FALSE, estimand=estimand, multinom = multinom, fillNAs = fillNAs)


   # this keeps the variables in the same order as vars
   n.rows <- sapply(ret,nrow)
   var.levels <- unlist(sapply(ret, rownames))
   var.names <- rep(names(ret),n.rows)
   var.names[var.levels!=""] <- paste(var.names[var.levels!=""],
                                      var.levels[var.levels!=""],sep=":")

   res <- data.frame(matrix(0,nrow=length(var.names), ncol=ncol(ret[[1]])))
   names(res) <- colnames(ret[[1]])
   rownames(res) <- var.names

   # populate the results table
   i.insert <- 1
   for(i in 1:length(ret))
   {
      res[i.insert:(i.insert+nrow(ret[[i]])-1),] <- ret[[i]]
      i.insert <- i.insert+nrow(ret[[i]])
   }

   res <- list(results=data.frame(res))
   return(res)
}

