

ci.rc <- function(b.k, SE.b.k=NULL, s.Y=NULL, s.X=NULL, N,K , R2.Y_X=NULL, R2.k_X.without.k=NULL, conf.level=.95, R2.Y_X.without.k=NULL, t.value=NULL, alpha.lower=NULL, alpha.upper=NULL, Noncentral=FALSE, Suppress.Statement=FALSE, ...)
{
return(ci.reg.coef(b.j=b.k, SE.b.j=SE.b.k, s.Y=s.Y, s.X=s.X, N=N, p=K, R2.Y_X=R2.Y_X, R2.j_X.without.j=R2.k_X.without.k, conf.level=conf.level, R2.Y_X.without.j=R2.Y_X.without.k, t.value=t.value, alpha.lower=alpha.lower, alpha.upper=alpha.upper, Noncentral=Noncentral, Suppress.Statement=Suppress.Statement, ...))
}
