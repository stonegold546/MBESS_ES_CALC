ci.c <- function(means=NULL, error.variance=NULL, c.weights=NULL, n=NULL, N=NULL,
Psi=NULL, conf.level=.95, alpha.lower=NULL, alpha.upper=NULL, df.error=NULL, ...)
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
}


if(!is.null(means))
{
Psi <- sum(c.weights*means)
}

if(is.null(alpha.lower) & is.null(alpha.upper))
{
alpha.lower <- (1-conf.level)/2
alpha.upper <- (1-conf.level)/2
}

if(is.null(N)) stop("You must specify the total sample size ('N').")
if(is.null(df.error)) df.2 <- N - length(means)

CV.up <- qt(1-alpha.upper, df=df.2, ncp = 0, lower.tail = TRUE, log.p = FALSE)
CV.low <- qt(alpha.lower, df=df.2, ncp = 0, lower.tail = TRUE, log.p = FALSE) 

Result <- list(Lower.Conf.Limit.Contrast = Psi + CV.low*part.of.se*sqrt(error.variance),
Contrast = Psi, Upper.Conf.Limit.Contrast = Psi + CV.up*part.of.se*sqrt(error.variance))
        
print(paste("The", 1 - (alpha.lower + alpha.upper), "confidence limits for the contrast are given as:"))
return(Result)
}
