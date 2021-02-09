#' Stop methods (e.g. "es.mean", "ks.mean", etc.) object, used only for backward compatibility
#'
#' In older versions of twang, the `ps` function specified the `stop.method` in a different
#' manner. This `stop.methods` object is used to ensure backward compatibility; new twang
#' users should not make use of it.
#'
#' This is merely a vector with the names of the stopping rules.
#'
#' @export
stop.methods <- matrix(c("es.mean","ks.mean","es.max","ks.max","ks.max.direct","es.max.direct"), nr = 1)
names(stop.methods) <- c("es.stat.mean","ks.stat.mean", "es.stat.max", "ks.stat.max",
"ks.max.direct", "es.max.direct")


ks.mean <- list(metric=ksStat,
                rule.summary=mean,
                direct=FALSE,
                estimand = NULL,
                na.action="level",
                name=NULL)
                
class(ks.mean) <- "stop.method"
                
es.mean <- list(metric=esStat,
                rule.summary=mean,
                direct=FALSE,
                estimand = NULL,
                na.action="level",
                name=NULL)
                
class(es.mean) <- "stop.method"
                
ks.max <- list(metric=ksStat,
				rule.summary=max,
				direct=FALSE,
				estimand = NULL,
				na.action="level",
				name=NULL)
				
class(ks.max) <- "stop.method"
				
es.max <- list(metric=esStat,
				rule.summary=max,
				direct=FALSE,
				estimand = NULL,
				na.action="level",
				name=NULL)
				
class(es.max) <- "stop.method"
				
ks.max.direct <- list(metric=ksStat,
				rule.summary=max,
				direct=TRUE,
				estimand= NULL,
				na.action="level",
				name=NULL)
				
class(ks.max.direct) <- "stop.method"
				
es.max.direct <- list(metric=esStat,
				rule.summary=max,
				direct=TRUE,
				estimand= NULL,
				na.action="level",
				name=NULL)
				
class(es.max.direct) <- "stop.method"