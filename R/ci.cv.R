ci.cv <- function(mean=NULL, sd=NULL, n=NULL, data=NULL, cv=NULL, conf.level=.95, alpha.lower=NULL, alpha.upper=NULL, ...)
{

if(is.null(conf.level))
{
if(alpha.lower>=1 | alpha.lower<0) stop("\'alpha.lower\' is not correctly specified.")
if(alpha.upper>=1 | alpha.upper<0) stop("\'alpha.upper\' is not correctly specified.")
}

if(!is.null(conf.level))
{
if(!is.null(alpha.lower) | !is.null(alpha.upper)) stop("Since \'conf.level\' is specified, \'alpha.lower\' and \'alpha.upper\' should be \'NULL\'.")

alpha.lower <- (1-conf.level)/2
alpha.upper <- (1-conf.level)/2

}

if(is.null(data) & is.null(cv))
{
if(is.null(mean)) stop("Either input the whole data set using \'data\' or specify the sample mean.")
if(mean<=0) stop("The sample mean must be some non-zero positive value (Does taking the absolute value of the mean make sense?).")

if(is.null(sd)) stop("Either input the whole data set using \'data\' or specify the sample standard deviation (using \'n\'-1 in the denominator).")
if(is.null(n)) stop("Either input the whole data set using \'data\' or specify the sample size")

k <- sd/mean
ncp.estimate <- sqrt(n)/k # or (sqrt(n)*x.bar)/sd

CI.NCP <- conf.limits.nct(ncp=ncp.estimate, df=n-1, alpha.lower=alpha.lower, alpha.upper=alpha.upper)

Low.Lim <- CI.NCP$Lower.Limit*((k^2)/sqrt(n))
Up.Lim <- CI.NCP$Upper.Limit*((k^2)/sqrt(n))

Result <- c(Lower.Limit.CofV=Low.Lim, Prob.Less.Lower=alpha.lower, Upper.Limit.CofV=Up.Lim, Prob.Greater.Upper=alpha.upper, C.of.V=k)

return(Result)
}

if(!is.null(data) & is.null(cv))
{

if(!is.null(mean)) stop("Since \'data\' is specified, do not specify the \'mean\'.")
if(!is.null(sd)) stop("Since \'data\' is specified, do not specify the \'sd\'.")
if(!is.null(n)) stop("Since \'data\' is specified, do not specify the \'n\'.")

n <- length(data)
sd.data <- (var(data))^.5
mean.data <- mean(data)
k <- sd.data/mean.data

print(n)

ncp.estimate <- sqrt(n)/k # or (sqrt(n)*x.bar)/sd

CI.NCP <- conf.limits.nct(ncp=ncp.estimate, df=n-1, alpha.lower=alpha.lower, alpha.upper=alpha.upper)

Low.Lim <- CI.NCP$Lower.Limit*((k^2)/sqrt(n))
Up.Lim <- CI.NCP$Upper.Limit*((k^2)/sqrt(n))

Result <- c(Lower.Limit.CofV=Low.Lim, Prob.Less.Lower=alpha.lower, Upper.Limit.CofV=Up.Lim, Prob.Greater.Upper=alpha.upper, C.of.V=k)

return(Result)
}

if(!is.null(cv))
{

k <- cv

if(is.null(n)) stop("Since you specified the coefficient of variation directly (\'cv\'), you must specify the sample size.")

if(!is.null(data)) stop("Since you specified the coefficient of variation (\'cv\') directly, do not include the raw data.")
if(!is.null(mean)) stop("Since you specified the coefficient of variation (\'cv\') directly, do not specify the mean (\'mean\').")
if(!is.null(sd)) stop("Since you specified the coefficient of variation (\'cv\') directly, do not specify the standard deviation (\'sd\').")

ncp.estimate <- sqrt(n)/k

CI.NCP <- conf.limits.nct(ncp=ncp.estimate, df=n-1, alpha.lower=alpha.lower, alpha.upper=alpha.upper)

Low.Lim <- CI.NCP$Lower.Limit*((k^2)/sqrt(n))
Up.Lim <- CI.NCP$Upper.Limit*((k^2)/sqrt(n))

Result <- c(Lower.Limit.CofV=Low.Lim, Prob.Less.Lower=alpha.lower, Upper.Limit.CofV=Up.Lim, Prob.Greater.Upper=alpha.upper, C.of.V=k)

return(Result)
}

}
