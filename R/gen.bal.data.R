gen.bal.data <- function(data,var.names){

  # keeps track of factor variables and numeric variables
  factor.vars = var.names[sapply(data[,var.names],class)=="factor"]
  numeric.vars = var.names[sapply(data[,var.names],class)!="factor"] # numeric or integer
  
  # construct a dataset that is used to optimize balance 
  # numeric vars must be first columns to ensure names are not altered due to factor expansion
  bal.data = data[,numeric.vars,drop=F]
  
  # expand factor variables into series of indicators
  for (vars in factor.vars){
    # na.pass doesnt drop observations with missing values
    form = as.formula(paste0("~-1+",vars))
    temp = model.matrix(form, model.frame(form,data=data, na.action=na.pass))
    # control naming convention of columns
    colnames(temp) = paste(vars , gsub(vars,"",colnames(temp)) , sep=":")
    # set missingness to 0. for factors, they are treated as own level.
    temp[is.na(temp)] = 0
    bal.data = cbind( bal.data , temp )
    rm(temp)
  }
  
  # create missing data indicators 
  NAdata = is.na(data[,c(numeric.vars,factor.vars)])
  colnames(NAdata) = paste0(colnames(NAdata),":<NA>")
  bal.data = cbind(bal.data , NAdata[,apply(NAdata,2,any),drop=F])
  
  # handles factors or NAs that have names that expand to the same name as a numeric
  colnames(bal.data) = make.unique(colnames(bal.data))
  
  return(list(bal.data=bal.data , numeric.vars=numeric.vars , factor.vars=factor.vars))
}
