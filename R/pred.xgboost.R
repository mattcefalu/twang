pred.xgboost <- function(preds, dtrain) {
   return(list(metric = "pred", value = (preds) ))
}