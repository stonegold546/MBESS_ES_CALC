ss.aipe.rc.sensitivity <-function(True.Var.Y=NULL, True.Cov.YX=NULL, True.Cov.XX=NULL, 
Estimated.Var.Y=NULL, Estimated.Cov.YX=NULL, Estimated.Cov.XX=NULL, Specified.N=NULL, 
which.predictor=1, w=NULL, Noncentral=FALSE, Standardize=FALSE, conf.level=.95, 
degree.of.certainty=NULL, G=1000, print.iter=TRUE)
{
return(ss.aipe.reg.coef.sensitivity(True.Var.Y=True.Var.Y,True.Cov.YX=True.Cov.YX, True.Cov.XX=True.Cov.XX, 
Estimated.Var.Y=Estimated.Var.Y, Estimated.Cov.YX=Estimated.Cov.YX, Estimated.Cov.XX=Estimated.Cov.XX, Specified.N=Specified.N, 
which.predictor=which.predictor, w=w, Noncentral=Noncentral, Standardize=Standardize, conf.level=conf.level, 
degree.of.certainty=degree.of.certainty, G=G, print.iter=print.iter))
}
