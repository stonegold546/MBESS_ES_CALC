ci.R2 <-
function(R2=NULL, df.1=NULL, df.2=NULL, conf.level=NULL, F.value=NULL, N=NULL, p=NULL, alpha.lower=NULL, alpha.upper=NULL, tol=1e-9)
{
if((!is.null(N) | !is.null(p)) & (!is.null(df.1) | !is.null(df.2))) stop("Either specify \'df.1\' and \'df.2\' or \'N\' and \'p,\' but not both combinations.")

if(!is.null(N) & !is.null(p) & is.null(df.1) & is.null(df.2))
{
df.1 <- p
df.2 <- N-p-1
}

if(!is.null(df.1) & !is.null(df.2) & is.null(N) & is.null(p))
{
N <- df.1 + df.2 + 1
p <- df.1
}

if(is.null(F.value))
{
F.value <- Rsquare2F(R2=R2, df.1=df.1, df.2=df.2, p=p, N=N)
}

Limits <- conf.limits.ncf(F.value=F.value, df.1=df.1, df.2=df.2, conf.level=conf.level, tol=tol, alpha.lower=alpha.lower, alpha.upper=alpha.upper)

if(length(Limits)==4)
{
LL <- Lambda2Rsquare(Limits$Lower.Limit, N=N)
Prob.LL <- Limits$Prob.Less.Lower
UL <- Lambda2Rsquare(Limits$Upper.Limit, N=N)
Prob.UL <- Limits$Prob.Greater.Upper

if(is.na(Limits$Lower.Limit))
{
LL <- 0
Prob.LL <-0
}

if(is.na(Limits$Upper.Limit))
{
UL <- 1
Prob.UL <- 0
}
return(list(Lower.Conf.Limit.R2=LL, Prob.Less.Lower=Prob.LL, Upper.Conf.Limit.R2=UL, Prob.Greater.Upper=Prob.UL))
}

if(length(Limits)==2 & is.null(Limits$Upper.Limit))
{
return(list(Lower.Conf.Limit.R2=Lambda2Rsquare(Limits$Lower.Limit, N=N), Prob.Less.Lower=1-Limits$Prob.Less.Lower))
}

if(length(Limits)==2 & is.null(Limits$Lower.Limit))
{
return(list(Upper.Conf.Limit.R2=Lambda2Rsquare(Limits$Upper.Limit, N=N), Prob.Greater.Upper=1-Limits$Prob.Greater.Upper))
}
}
