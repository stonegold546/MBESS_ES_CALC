ci.sc <- function(means=NULL, error.variance=NULL, c.weights=NULL, n=NULL, N=NULL,
Psi=NULL, ncp=NULL, conf.level=.95, alpha.lower=NULL, alpha.upper=NULL, df.error=NULL, ...)
{

if(length(n)==1) 
{
n <- rep(n, length(means))
}

if(length(n)!=length(c.weights)) stop("The lengths of 'n' and 'c.weights' differ, which should not be the case.")
if(!(sum(c.weights)==0)) stop("The sum of the contrast weights ('c.weights') should equal zero.")

part.of.se <- sqrt(sum((c.weights^2)/n))

if(!is.null(Psi))
{
if(!is.null(means)) stop("Since the contrast effect ('Psi') was specified, you should not specify the vector of means ('means').")
if(!is.null(ncp)) stop("Since the contrast effect ('Psi') was specified, you should not specify the noncentral parameter ('ncp').")
if(is.null(error.variance)) stop("You must specify the error variance ('error.variance').")
if(is.null(n)) stop("You must specify the vector per group/level sample size ('n').")
if(is.null(c.weights)) stop("You must specify the vector of contrast weights ('c.weights').")

psi <- Psi/sqrt(error.variance)
lambda <- psi*part.of.se
}

if(!is.null(ncp))
{
if(!is.null(means)) stop("Since the noncentral parameter was specified directly, you should not specify the vector of means ('means').")
if(!is.null(Psi)) stop("Since the noncentral parameter was specified directly, you should not specify the the contrast effect ('Psi').")
if(is.null(error.variance)) stop("You must specify the error variance ('error.variance').")
if(is.null(n)) stop("You must specify the vector per group/level sample size ('n').")
if(is.null(c.weights)) stop("You must specify the vector of contrast weights ('c.weights'.")

lambda <- ncp
}

if(!is.null(means))
{
Psi <- sum(c.weights*means)
psi <- Psi/sqrt(error.variance)

lambda <- psi/part.of.se

}

if(is.null(alpha.lower) & is.null(alpha.upper))
{
alpha.lower <- (1-conf.level)/2
alpha.upper <- (1-conf.level)/2
}

if(is.null(N)) stop("You must specify the total sample size ('N').")
if(is.null(df.error)) df.2 <- N - length(means)

Lims <- conf.limits.nct(ncp=lambda, df=df.2, conf.level = NULL, alpha.lower = alpha.lower, 
        alpha.upper = alpha.upper, sup.int.warns=TRUE, method = "all", ...)

Result <- list(Lower.Conf.Limit.Standardized.Contrast = Lims$Lower.Limit*part.of.se, Standardized.contrast = psi, 
        Upper.Conf.Limit.Standardized.Contrast = Lims$Upper.Limit*part.of.se)
        
print(paste("The", 1 - (alpha.lower + alpha.upper), "confidence limits for the standardized contrast are given as:"))

return(Result)
}
