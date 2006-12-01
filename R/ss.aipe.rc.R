ss.aipe.rc <- function(Rho2.Y_X=NULL, Rho2.k_X.without.k=NULL, K=NULL, b.k=NULL, width, which.width="Full", sigma.Y=1, sigma.X.k=1, RHO.XX=NULL, Rho.YX=NULL, which.predictor=NULL, alpha.lower=NULL, alpha.upper=NULL, conf.level=1-.05, degree.of.certainty=NULL, Suppress.Statement=FALSE)
{
return(ss.aipe.reg.coef(Rho2.Y_X=Rho2.Y_X, Rho2.j_X.without.j=Rho2.k_X.without.k, p=K, b.j=b.k, width=width, which.width=which.width, sigma.Y=sigma.Y, sigma.X=sigma.X.k, RHO.XX=RHO.XX, Rho.YX=Rho.YX, which.predictor=which.predictor, Noncentral=FALSE, alpha.lower=alpha.lower, alpha.upper=alpha.upper, conf.level=conf.level, degree.of.certainty=degree.of.certainty, Suppress.Statement=Suppress.Statement))
}
