surveytable <- function(formula , data ){
  rows<-formula[[2]][[2]]
  cols<-formula[[2]][[3]]
  rowvar<-unique(data[,as.character(rows)])
  colvar<-unique(data[,as.character(cols)])

  nr<-length(rowvar)
  nc<-length(colvar)
  
  fsat<-eval(bquote(~interaction(factor(.(rows)),factor(.(cols)))-1))
  mm<-model.matrix(fsat,model.frame(fsat, data,na.action=na.pass))
  N<-nrow(mm)
  nu <- N-1
  
  
  p = xtabs(eval(bquote(w~.(rows)+.(cols))),data=data)*N/sum(data$w)
  pearson <- chisq.test( p , correct=F)
  
  mf1<-expand.grid(rows=1:nr,cols=1:nc)
  X1<-model.matrix(~factor(rows)+factor(cols),mf1)
  X12<-model.matrix(~factor(rows)*factor(cols),mf1)
  
  

  mean2 = as.vector(data$w%*%mm / sum(data$w))
  score = (sweep(mm , 2  , mean2 , "-")) / sum(data$w)
  a = sweep(score,1,data$w,"*")
  a = sweep(a,2, colSums(a)/N , "-")
  V = t(a)%*%a*N/(N-1)
  
  Cmat<-qr.resid(qr(X1),X12[,-(1:(nr+nc-1)),drop=FALSE])
  Dmat <- diag(mean2)
  iDmat<- diag(ifelse(mean2==0,0,1/mean2))
  Vsrs <- (Dmat - outer(mean2,mean2))/N
  denom<- t(Cmat) %*% (iDmat/N) %*% Cmat
  numr<-t(Cmat)%*% iDmat %*% V %*% iDmat %*% Cmat
  Delta<-solve(denom,numr)
  d0<- sum(diag(Delta))^2/(sum(diag(Delta%*%Delta)))
  
  pearson$statistic<-pearson$statistic/sum(diag(Delta))
  pearson$p.value<-pf(pearson$statistic, d0, d0*nu, lower.tail=FALSE)
  attr(pearson$statistic,"names")<-"F"
  pearson$parameter<-c(ndf=d0,ddf=d0*nu)
  pearson$method<-"Pearson's X^2: Rao & Scott adjustment"
  pearson$p <- p
  
  return(pearson)
}