upsilon <- function(x, mediator, dv, conf.level = 0.95, bootstrap = TRUE, B = 1000, boot.type = "ordinary", ...){
Data <- data.frame(x=x,mediator=mediator,dv=dv)
med.model <- paste0(paste0(colnames(Data)[2],'~',colnames(Data)[1]),' \n ',paste0(colnames(Data)[3],'~',colnames(Data)[1],'+',colnames(Data)[2]))
med.res <- lavaan::sem(med.model,Data)
upsilon <- lavaan::coef(med.res)[1]**2*lavaan::coef(med.res)[3]**2*var(Data[,colnames(Data)[1]],na.rm=TRUE)/var(Data[,colnames(Data)[3]],na.rm=TRUE)
adj.upsilon <- (lavaan::coef(med.res)[1]**2-lavaan::vcov(med.res)[1,1])*(lavaan::coef(med.res)[3]**2-lavaan::vcov(med.res)[3,3])*var(Data[,colnames(Data)[1]],na.rm=TRUE)/var(Data[,colnames(Data)[3]],na.rm=TRUE)
if(bootstrap==TRUE){
  med.boot.fun <- function(out){
data <- lavaan::lavInspect(out,what = "data")
upsilon.boot <- lavaan::coef(out)[1]**2*lavaan::coef(out)[3]**2*var(data[,1],na.rm=TRUE)/var(data[,3],na.rm=TRUE)
adj.upsilon.boot <- (lavaan::coef(out)[1]**2-lavaan::vcov(out)[1,1])*(lavaan::coef(out)[3]**2-lavaan::vcov(out)[3,3])*var(data[,1],na.rm=TRUE)/var(data[,3],na.rm=TRUE)
as.vector(c(upsilon.boot,adj.upsilon.boot))
  }
  boot.med.res <- lavaan::bootstrapLavaan(med.res,R=B,type=boot.type,FUN=med.boot.fun)
  upsES.out <- data.frame(est=c(upsilon,adj.upsilon),
  lcl=c(ups.LCL=quantile(boot.med.res[,1],probs=(1-conf.level)/2),
adj.ups.UCL=quantile(boot.med.res[,2],probs=(1-conf.level)/2)),
  ucl=c(ups.LCL=quantile(boot.med.res[,1],probs=1-(1-conf.level)/2),
adj.ups.UCL=quantile(boot.med.res[,2],probs=1-(1-conf.level)/2)),
  row.names=c('Upsilon','Adj Upsilon'))
  colnames(upsES.out) <- c('Estimate',paste0(conf.level*100,'% LCL'),paste0(conf.level*100,'% UCL'))
  return(upsES.out)
} else {
  upsES.out <- data.frame(Est=c(upsilon,adj.upsilon),row.names=c('Upsilon','Adj Upsilon'))
  colnames(upsES.out) <- 'Estimate'
  return(upsES.out)
}
}
