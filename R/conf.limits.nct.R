conf.limits.nct <- function(ncp, df, conf.level=.95, alpha.lower=NULL, alpha.upper=NULL, tol=1e-9, sup.int.warns=TRUE, method="all", ...)
{

if(!(method=="all" | method=="ALL" | method=="All" | method==1 | method==2 | method==3)) stop("You need to specify the method to use; the default is \'all\'")

if(df <= 0) stop("The degrees of freedom must be some positive value.", call.=FALSE)
if((!is.null(alpha.lower) | !is.null(alpha.upper)) & !is.null(conf.level)) stop("You must choose either to use \'conf.level\' or define the \'lower.alpha\' and \'upper.alpha\' values, but not both", call.=FALSE)

if(!is.null(conf.level))
{
alpha.lower <- (1-conf.level)/2
alpha.upper <- (1-conf.level)/2
}

if(method=="all" | method=="ALL" | method=="All")
{
Res.M1 <- NULL
Res.M2 <- NULL
Res.M3 <- NULL

try(Res.M1 <- conf.limits.nct.M1(ncp=ncp, df=df, conf.level=NULL, alpha.lower=alpha.lower, alpha.upper=alpha.upper, tol=tol), silent=TRUE)
if(length(Res.M1)!=4) Res.M1 <- NULL

try(Res.M2 <- conf.limits.nct.M2(ncp=ncp, df=df, conf.level=NULL, alpha.lower=alpha.lower, alpha.upper=alpha.upper, tol=tol), silent=TRUE)
if(length(Res.M2)!=4) Res.M2 <- NULL

try(Res.M3 <- conf.limits.nct.M3(ncp=ncp, df=df, conf.level=NULL, alpha.lower=alpha.lower, alpha.upper=alpha.upper, tol=tol), silent=TRUE)
if(length(Res.M3)!=4) Res.M3 <- NULL

Low.M1 <- Res.M1$Lower.Limit
Prob.Low.M1 <- Res.M1$Prob.Less.Lower
Upper.M1 <- Res.M1$Upper.Limit
Prob.Upper.M1 <- Res.M1$Prob.Greater.Upper

Low.M2 <- Res.M2$Lower.Limit
Prob.Low.M2 <- Res.M2$Prob.Less.Lower
Upper.M2 <- Res.M2$Upper.Limit
Prob.Upper.M2 <- Res.M2$Prob.Greater.Upper

Low.M3 <- Res.M3$Lower.Limit
Prob.Low.M3 <- Res.M3$Prob.Less.Lower
Upper.M3 <- Res.M3$Upper.Limit
Prob.Upper.M3 <- Res.M3$Prob.Greater.Upper

Min.for.Best.Low <- min((c(Prob.Low.M1, Prob.Low.M2, Prob.Low.M3)-alpha.lower)^2)

if(!is.null(Res.M1)){if(Min.for.Best.Low==(Prob.Low.M1-alpha.lower)^2) Best.Low <- 1}
if(!is.null(Res.M2)){if(Min.for.Best.Low==(Prob.Low.M2-alpha.lower)^2) Best.Low <- 2}
if(!is.null(Res.M3)){if(Min.for.Best.Low==(Prob.Low.M3-alpha.lower)^2) Best.Low <- 3}

Min.for.Best.Up <- min((c(Prob.Upper.M1, Prob.Upper.M2, Prob.Upper.M3)-alpha.upper)^2)

if(!is.null(Res.M1)){if(Min.for.Best.Up==(Prob.Upper.M1-alpha.upper)^2) Best.Up <- 1}
if(!is.null(Res.M2)){if(Min.for.Best.Up==(Prob.Upper.M2-alpha.upper)^2) Best.Up <- 2}
if(!is.null(Res.M3)){if(Min.for.Best.Up==(Prob.Upper.M3-alpha.upper)^2) Best.Up <- 3}

if(is.null(Res.M1)) {Low.M1 <- NA; Prob.Low.M1 <- NA; Upper.M1 <- NA; Prob.Upper.M1 <- NA}
if(is.null(Res.M2)) {Low.M2 <- NA; Prob.Low.M2 <- NA; Upper.M2 <- NA; Prob.Upper.M2 <- NA}
if(is.null(Res.M3)) {Low.M3 <- NA; Prob.Low.M3 <- NA; Upper.M3 <- NA; Prob.Upper.M3 <- NA}

Result <- list(Lower.Limit=c(Low.M1, Low.M2, Low.M3)[Best.Low], Prob.Less.Lower=c(Prob.Low.M1, Prob.Low.M2, Prob.Low.M3)[Best.Low], Upper.Limit=c(Upper.M1, Upper.M2, Upper.M3)[Best.Up], Prob.Greater.Upper=c(Prob.Upper.M1, Prob.Upper.M2, Prob.Upper.M3)[Best.Up])
}

if(method==1)
{
Res.M1 <- NULL

try(Res.M1 <- conf.limits.nct.M1(ncp=ncp, df=df, conf.level=NULL, alpha.lower=alpha.lower, alpha.upper=alpha.upper, tol=tol), silent=TRUE)
if(length(Res.M1)!=4) Res.M1 <- NULL

if(is.null(Res.M1)) stop("There was a problem with Method 1 in this situation; try another method.")

Result <- list(Lower.Limit=Res.M1$Lower.Limit, Prob.Less.Lower=Res.M1$Prob.Less.Lower, Upper.Limit=Res.M1$Upper.Limit, Prob.Greater.Upper=Res.M1$Prob.Greater.Upper)
}

if(method==2)
{
Res.M2 <- NULL

try(Res.M2 <- conf.limits.nct.M2(ncp=ncp, df=df, conf.level=NULL, alpha.lower=alpha.lower, alpha.upper=alpha.upper, tol=tol), silent=TRUE)
if(length(Res.M2)!=4) Res.M2 <- NULL

if(is.null(Res.M2)) stop("There was a problem with Method 2 in this situation; try another method.")

Result <- list(Lower.Limit=Res.M2$Lower.Limit, Prob.Less.Lower=Res.M2$Prob.Less.Lower, Upper.Limit=Res.M2$Upper.Limit, Prob.Greater.Upper=Res.M2$Prob.Greater.Upper)
}

if(method==3)
{
Res.M3 <- NULL

try(Res.M3 <- conf.limits.nct.M3(ncp=ncp, df=df, conf.level=NULL, alpha.lower=alpha.lower, alpha.upper=alpha.upper, tol=tol), silent=TRUE)
if(length(Res.M3)!=4) Res.M3 <- NULL

if(is.null(Res.M3)) stop("There was a problem with Method 3 in this situation; try another method.")

Result <- list(Lower.Limit=Res.M3$Lower.Limit, Prob.Less.Lower=Res.M3$Prob.Less.Lower, Upper.Limit=Res.M3$Upper.Limit, Prob.Greater.Upper=Res.M3$Prob.Greater.Upper)
}

return(Result)
}
