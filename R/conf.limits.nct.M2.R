conf.limits.nct.M2 <- function(ncp=ncp, df=df, conf.level=.95, alpha.lower=NULL, alpha.upper=NULL, tol=1e-9, sup.int.warns=TRUE, ...)
{
if(sup.int.warns==TRUE) Orig.warn <- options()$warn; options(warn=-1)

if(ncp < 0) stop("\'conf.limits.nct.M2\' should not be given negative values for the noncentral parameter. Use \'conf.limits.nct\' instead.")

# No longer necessary, given the above restriction.
ncp.orig <- ncp
ncp <- abs(ncp)

if(df <= 0) stop("The degrees of freedom must be some positive value.", call.=FALSE)

if((!is.null(alpha.lower) | !is.null(alpha.upper)) & !is.null(conf.level)) stop("You must choose either to use \'conf.level\' or define the \'lower.alpha\' and \'upper.alpha\' values, but not both", call.=FALSE)

if(!is.null(conf.level))
{
alpha.lower <- (1-conf.level)/2
alpha.upper <- (1-conf.level)/2
}




if(ncp.orig < 0)
{
.tmp.MBESS <- alpha.lower
alpha.lower <- alpha.upper
alpha.upper <- .tmp.MBESS
}


# Internal function for upper limit.
###########################
.ci.nct.lower <- function(val.of.interest, ...)
{
(qt(p=alpha.lower, df=df, ncp=val.of.interest, lower.tail = FALSE, log.p = FALSE) - ncp)^2
}
###########################

# Internal function for lower limit.
###########################
.ci.nct.upper <- function(val.of.interest, ...)
{
(qt(p=alpha.upper, df=df, ncp=val.of.interest, lower.tail = TRUE, log.p = FALSE) - ncp)^2
}

if(sup.int.warns==TRUE)
{
Low.Lim <- suppressWarnings(nlm(f=.ci.nct.lower, p=ncp, ...))
Up.Lim <- suppressWarnings(nlm(f=.ci.nct.upper, p=ncp, ...))
}

if(sup.int.warns==FALSE)
{
Low.Lim <- nlm(f=.ci.nct.lower, p=ncp, ...)
Up.Lim <- nlm(f=.ci.nct.upper, p=ncp, ...)
}

if(ncp.orig >= 0)
{
if(alpha.lower==0) Result <- list(Lower.Limit=-Inf, Prob.Less.Lower=0, Upper.Limit=Up.Lim$estimate, Prob.Greater.Upper=pt(q=ncp, ncp=Up.Lim$estimate, df=df))
if(alpha.upper==0) Result <- list(Lower.Limit=Low.Lim$estimate, Prob.Less.Lower=pt(q=ncp, ncp=Low.Lim$estimate, df=df, lower.tail=FALSE), Upper.Limit=Inf, Prob.Greater.Upper=0)
if(alpha.lower!=0 & alpha.upper!=0) Result <- list(Lower.Limit=Low.Lim$estimate, Prob.Less.Lower=pt(q=ncp, ncp=Low.Lim$estimate, df=df, lower.tail=FALSE), Upper.Limit=Up.Lim$estimate, Prob.Greater.Upper=pt(q=ncp, ncp=Up.Lim$estimate, df=df))
}

if(ncp.orig < 0)
{
# Recall that the alpha values were reversed for negative noncentral parameters.
if(alpha.upper==0) Result <- list(Lower.Limit=-Inf, Prob.Less.Lower=0, Upper.Limit=-Low.Lim$estimate, Prob.Greater.Upper=pt(q=ncp, ncp=Low.Lim$estimate, df=df, lower.tail=FALSE))
if(alpha.lower==0) Result <- list(Lower.Limit=-Up.Lim$estimate, Prob.Less.Lower=pt(q=ncp, ncp=Up.Lim$estimate, df=df), Upper.Limit=Inf, Prob.Greater.Upper=0)
if(alpha.lower!=0 & alpha.upper!=0) Result <- list(Lower.Limit=-Up.Lim$estimate, Prob.Less.Lower=pt(q=ncp, ncp=Up.Lim$estimate, df=df), Upper.Limit=-Low.Lim$estimate, Prob.Greater.Upper=pt(q=ncp, ncp=Low.Lim$estimate, df=df, lower.tail=FALSE))
}
if(sup.int.warns==TRUE) options(warn=Orig.warn)

return(Result)
}
