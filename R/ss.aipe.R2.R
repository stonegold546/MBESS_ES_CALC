ss.aipe.R2 <-
function(Population.R2=NULL, conf.level=.95, width=NULL, Random.Predictors=TRUE, which.width="Full", p=NULL, degree.of.certainty=NULL, Tol=1e-9, start.N=NULL, ...)
{
char.expand(which.width, c("Full", "Lower", "Upper"), nomatch = stop("Problems with 'which.width' specification. You must choose either 'Full', 'Lower', or 'Upper'.", call.=FALSE))
if(is.null(p)) stop("You need to specify \'p\', the number of predictors.")
if(is.null(Population.R2)) stop("You need to specify the population squared multiple correlation coefficient, \'Population.R2\'.")

if(!is.null(degree.of.certainty)){if(degree.of.certainty < .50 | degree.of.certainty >= 1) stop("The \'degree.of.certainty\' must be greater than .50 and less than 1.")}

# Expected value of R^2 given Population.R2, N, and p.
Expected.R2 <- function(Population.R2, N, p)
{
# From Kendall's Advanced Theory (1999, p. 531).
Value <- 1 - ((N-p-1)/(N-1))*(1-Population.R2)*hyperg_2F1(1, 1, .5*(N+1), Population.R2)
Value <- max(0, Value)
return(Value)
}
#######

alpha.lower <- alpha.upper <- (1-conf.level)/2

#Common for all methods in order to obtain a starting sample size.
######################################################################
N.0 <- p + 1 + p # p denominator DF to begin.
if(!is.null(start.N)) N.0 <- max(N.0, ceiling(start.N))

Continue <- TRUE
while(Continue==TRUE)
{
#Using Random.Predictors=FALSE here to get an estimate for sample size.
CI.0 <- ci.R2(R2=Expected.R2(Population.R2=Population.R2, N=N.0, p=p), alpha.lower=alpha.lower, alpha.upper=alpha.upper, N=N.0, p=p, Random.Predictors=FALSE, ...)
Continue <- sum(is.na(CI.0))>0
N.0 <- N.0 + 1
}
######################################################################
if(which.width=="Full")
{


N.1 <- N.0 + 1
CI.1 <- ci.R2(R2=Expected.R2(Population.R2=Population.R2, N=N.1, p=p), alpha.lower=alpha.lower, alpha.upper=alpha.upper, N=N.1, p=p, Random.Predictors=FALSE)
w.F <- CI.1$Upper.Conf.Limit.R2 - CI.1$Lower.Conf.Limit.R2
Diff <- w.F - width
while(Diff > Tol | is.na(Diff))
{
N.1 <- N.1 + 1
CI.1 <- ci.R2(R2=Expected.R2(Population.R2=Population.R2, N=N.1, p=p), alpha.lower=alpha.lower, alpha.upper=alpha.upper, N=N.1, p=p, Random.Predictors=FALSE)
w.F <- CI.1$Upper.Conf.Limit.R2 - CI.1$Lower.Conf.Limit.R2
Diff <- w.F - width
}

if(Random.Predictors==TRUE)
{
# Starting value is from the method where the predictors were treated as fixed (from above procedure).
CI.1 <- ci.R2(R2=Expected.R2(Population.R2=Population.R2, N=N.1, p=p), alpha.lower=alpha.lower, alpha.upper=alpha.upper, N=N.1, p=p, Random.Predictors=TRUE)
w.F <- CI.1$Upper.Conf.Limit.R2 - CI.1$Lower.Conf.Limit.R2
Diff <- w.F - width
while(Diff > Tol | is.na(Diff))
{
N.1 <- N.1 + 1
CI.1 <- ci.R2(R2=Expected.R2(Population.R2=Population.R2, N=N.1, p=p), alpha.lower=alpha.lower, alpha.upper=alpha.upper, N=N.1, p=p, Random.Predictors=TRUE)
w.F <- CI.1$Upper.Conf.Limit.R2 - CI.1$Lower.Conf.Limit.R2
Diff <- w.F - width
}
}

Result.Full <- list(Required.Sample.Size=N.1, Expected.Width=w.F)

if(!is.null(degree.of.certainty)) 
{
Conf.Limit.Desired.Certainty.Upper <- ci.R2(R2=Expected.R2(Population.R2=Population.R2, N=N.1, p=p), alpha.lower=0, alpha.upper=1-degree.of.certainty, N=N.1, p=p, Random.Predictors=Random.Predictors)$Upper.Conf.Limit.R2
Conf.Limit.Desired.Certainty.Lower <- ci.R2(R2=Expected.R2(Population.R2=Population.R2, N=N.1, p=p), alpha.lower=1-degree.of.certainty, alpha.upper=0, N=N.1, p=p, Random.Predictors=Random.Predictors)$Lower.Conf.Limit.R2

N.2 <- N.1 + 1
CI.2 <- ci.R2(R2=Expected.R2(Population.R2=Conf.Limit.Desired.Certainty.Upper, N=N.2, p=p), alpha.lower=alpha.lower, alpha.upper=alpha.upper, N=N.2, p=p, Random.Predictors=Random.Predictors)
w.F <- CI.2$Upper.Conf.Limit.R2 - CI.2$Lower.Conf.Limit.R2
Diff <- w.F - width
while(Diff > Tol | is.na(Diff))
{
N.2 <- N.2 + 1
CI.2 <- ci.R2(R2=Expected.R2(Population.R2=Conf.Limit.Desired.Certainty.Upper, N=N.2, p=p), alpha.lower=alpha.lower, alpha.upper=alpha.upper, N=N.2, p=p, Random.Predictors=Random.Predictors)
w.F <- CI.2$Upper.Conf.Limit.R2 - CI.2$Lower.Conf.Limit.R2
Diff <- w.F - width
}
N.Upper.Conf.Lim <- N.2
Ex.Width.Upper <- w.F

N.3 <- N.1 + 1
CI.3 <- ci.R2(R2=Expected.R2(Population.R2=Conf.Limit.Desired.Certainty.Lower, N=N.3, p=p), alpha.lower=alpha.lower, alpha.upper=alpha.upper, N=N.3, p=p, Random.Predictors=Random.Predictors)
w.F <- CI.3$Upper.Conf.Limit.R2 - CI.3$Lower.Conf.Limit.R2
Diff <- w.F - width
while(Diff > Tol | is.na(Diff))
{
N.3 <- N.3 + 1
CI.3 <- ci.R2(R2=Expected.R2(Population.R2=Conf.Limit.Desired.Certainty.Lower, N=N.3, p=p), alpha.lower=alpha.lower, alpha.upper=alpha.upper, N=N.3, p=p, Random.Predictors=Random.Predictors)
w.F <- CI.3$Upper.Conf.Limit.R2 - CI.3$Lower.Conf.Limit.R2
Diff <- w.F - width
}
N.Lower.Conf.Lim <- N.3
Ex.Width.Lower <- w.F

if(N.Upper.Conf.Lim >= N.Lower.Conf.Lim)
{
Result.Full <- list(Required.Sample.Size=N.Upper.Conf.Lim, Expected.Width=Ex.Width.Upper)
}

if(N.Lower.Conf.Lim > N.Upper.Conf.Lim)
{
Result.Full <- list(Required.Sample.Size=N.Lower.Conf.Lim, Expected.Width=Ex.Width.Lower)
}


# To get maximum expected CI width to ensure adequate sample size is used.
CI.WIDTH.R2 <- function(R2, conf.level, N, p, Random.Predictors)
{
Lims <- ci.R2(R2=Expected.R2(Population.R2=R2, N, p), conf.level=conf.level, N=N, p=p, Random.Predictors=Random.Predictors)
Lims$Upper - Lims$Lower
}
Optimize.Result <- optimize(f=CI.WIDTH.R2, interval = c(Conf.Limit.Desired.Certainty.Lower, Conf.Limit.Desired.Certainty.Upper), maximum = TRUE, tol = .Machine$double.eps^0.5, N=N.1, p=p, conf.level=conf.level, Random.Predictors=Random.Predictors)$maximum

# The following is used if and only if the maximum confidence interval width is between the confidence limits.
if((round(Conf.Limit.Desired.Certainty.Lower,4) < round(Optimize.Result,4)) & (round(Optimize.Result,4) < round(Conf.Limit.Desired.Certainty.Upper,4)))
{

Prob.F.Max.Width <- pf(q=Rsquare2F(R2=Optimize.Result, p=p, N=N.1), df1=p, df2=N.1-p-1, ncp=Rsquare2Lambda(R2=Population.R2, N=N.1))

low.lim.F <- qf(p=max(0, (Prob.F.Max.Width-(1-degree.of.certainty)/2)), df1=p, df2=N.1-p-1, ncp=Rsquare2Lambda(R2=Population.R2, N=N.1))
up.lim.F <- qf(p=min((Prob.F.Max.Width+(1-degree.of.certainty)/2), 1), df1=p, df2=N.1-p-1, ncp=Rsquare2Lambda(R2=Population.R2, N=N.1))

low.R2 <- F2Rsquare(F.value=low.lim.F, df.1 = p, df.2 = N.1-p-1)
up.R2 <- F2Rsquare(F.value=up.lim.F, df.1 = p, df.2 = N.1-p-1)

# Here SS is determined for the upper limit from the optimize method.
if(up.R2>0 & up.R2<1)
{
N.4 <- N.1 + 1
CI.4 <- ci.R2(R2=Expected.R2(Population.R2=up.R2, N=N.4, p=p), alpha.lower=alpha.lower, alpha.upper=alpha.upper, N=N.4, p=p, Random.Predictors=Random.Predictors)
w.F <- CI.4$Upper.Conf.Limit.R2 - CI.4$Lower.Conf.Limit.R2
Diff <- w.F - width
while(Diff > Tol | is.na(Diff))
{
N.4 <- N.4 + 1
CI.4 <- ci.R2(R2=Expected.R2(Population.R2=up.R2, N=N.4, p=p), alpha.lower=alpha.lower, alpha.upper=alpha.upper, N=N.4, p=p, Random.Predictors=Random.Predictors)
w.F <- CI.4$Upper.Conf.Limit.R2 - CI.4$Lower.Conf.Limit.R2
Diff <- w.F - width
}
n.up.R2 <- N.4
Ex.Width.Upper.after.optim <- w.F
}

# Here SS is determined for the lower limit from the optimize method.
if(low.R2>0 & low.R2<1)
{
N.5 <- N.1 + 1
CI.5 <- ci.R2(R2=Expected.R2(Population.R2=low.R2, N=N.5, p=p), alpha.lower=alpha.lower, alpha.upper=alpha.upper, N=N.5, p=p, Random.Predictors=Random.Predictors)
w.F <- CI.5$Upper.Conf.Limit.R2 - CI.5$Lower.Conf.Limit.R2
Diff <- w.F - width
while(Diff > Tol | is.na(Diff))
{
N.5 <- N.5 + 1
CI.5 <- ci.R2(R2=Expected.R2(Population.R2=low.R2, N=N.5, p=p), alpha.lower=alpha.lower, alpha.upper=alpha.upper, N=N.5, p=p, Random.Predictors=Random.Predictors)
w.F <- CI.5$Upper.Conf.Limit.R2 - CI.5$Lower.Conf.Limit.R2
Diff <- w.F - width
}
n.low.R2 <- N.5
Ex.Width.Lower.after.optim <- w.F
}

# To keep function from failing when an impossible low.R2 arises.
if(low.R2==0 | low.R2==1) n.low.R2 <- N.1
if(up.R2==0 | up.R2==1) n.up.R2 <- N.1

if((n.low.R2 >= n.up.R2) & (n.low.R2 > Result.Full$Required.Sample.Size))
{
Result.Full <- list(Required.Sample.Size=n.low.R2, Expected.Width=Ex.Width.Lower.after.optim)
}

if((n.up.R2 > n.low.R2) & (n.up.R2 > Result.Full$Required.Sample.Size))
{
Result.Full <- list(Required.Sample.Size=n.up.R2, Expected.Width=Ex.Width.Upper.after.optim)
}
}

}

# Get the expected width.
tmp.MBESS <- ci.R2(R2=Expected.R2(Population.R2=Population.R2, N=Result.Full$Required.Sample.Size, p=p),
alpha.lower=alpha.lower, alpha.upper=alpha.upper, N=Result.Full$Required.Sample.Size, p=p, 
Random.Predictors=Random.Predictors)
Result.Full <- list(Required.Sample.Size=Result.Full$Required.Sample.Size, Expected.Width=tmp.MBESS$Upper - tmp.MBESS$Lower)
return(Result.Full)
}


if(which.width=="Lower")
{
N.1 <- N.0 + 1
CI.1 <- ci.R2(R2=Expected.R2(Population.R2=Population.R2, N=N.1, p=p), alpha.lower=alpha.lower, alpha.upper=alpha.upper, N=N.1, p=p, Random.Predictors=FALSE)
w.L <- Expected.R2(Population.R2=Population.R2, N=N.1, p=p) - CI.1$Lower.Conf.Limit.R2
Diff <- w.L - width
while(Diff > Tol | is.na(Diff))
{
N.1 <- N.1 + 1
CI.1 <- ci.R2(R2=Expected.R2(Population.R2=Population.R2, N=N.1, p=p), alpha.lower=alpha.lower, alpha.upper=alpha.upper, N=N.1, p=p, Random.Predictors=FALSE)
w.L <- Expected.R2(Population.R2=Population.R2, N=N.1, p=p) - CI.1$Lower.Conf.Limit.R2
Diff <- w.L - width
}

if(Random.Predictors==TRUE)
{
# Use N.1 as the starting point from above (which was based on fixed predictors).
CI.1 <- ci.R2(R2=Expected.R2(Population.R2=Population.R2, N=N.1, p=p), alpha.lower=alpha.lower, alpha.upper=alpha.upper, N=N.1, p=p, Random.Predictors=TRUE)
w.L <- Expected.R2(Population.R2=Population.R2, N=N.1, p=p) - CI.1$Lower.Conf.Limit.R2
Diff <- w.L - width
while(Diff > Tol | is.na(Diff))
{
N.1 <- N.1 + 1
CI.1 <- ci.R2(R2=Expected.R2(Population.R2=Population.R2, N=N.1, p=p), alpha.lower=alpha.lower, alpha.upper=alpha.upper, N=N.1, p=p, Random.Predictors=TRUE)
w.L <- Expected.R2(Population.R2=Population.R2, N=N.1, p=p) - CI.1$Lower.Conf.Limit.R2
Diff <- w.L - width
}
}

if(is.null(degree.of.certainty)) 
{
Result.Low <- list(Required.Sample.Size=N.1, Expected.Width=w.L)
print("This sample size should be regarded as approximate; the methods have only been throughly evaluated for Full widths.")
return(Result.Low)
}

if(!is.null(degree.of.certainty)) 
{
Conf.Limit.Desired.Certainty <- ci.R2(R2=Expected.R2(Population.R2=Population.R2, N=N.1, p=p), alpha.lower=0, alpha.upper=1-degree.of.certainty, N=N.1, p=p, Random.Predictors=Random.Predictors)$Upper.Conf.Limit.R2
N.2 <- N.1
CI.2 <- ci.R2(R2=Expected.R2(Population.R2=Conf.Limit.Desired.Certainty, N=N.2, p=p), alpha.lower=alpha.lower, alpha.upper=alpha.upper, N=N.2, p=p, Random.Predictors=Random.Predictors)
w.L <- Conf.Limit.Desired.Certainty - CI.2$Lower.Conf.Limit.R2
Diff <- w.L - width

while(Diff > Tol)
{
N.2 <- N.2 + 1
CI.2 <- ci.R2(R2=Expected.R2(Population.R2=Conf.Limit.Desired.Certainty, N=N.2, p=p), alpha.lower=alpha.lower, alpha.upper=alpha.upper, N=N.2, p=p, Random.Predictors=Random.Predictors)
w.L <- Conf.Limit.Desired.Certainty - CI.2$Lower.Conf.Limit.R2
Diff <- w.L - width
}
w.L.M.2 <- Expected.R2(Population.R2=Population.R2, N=N.2, p=p) - ci.R2(R2=Expected.R2(Population.R2=Population.R2, N=N.2, p=p), alpha.lower=alpha.lower, alpha.upper=alpha.upper, N=N.2, p=p, Random.Predictors=Random.Predictors)$Lower.Conf.Limit.R2

Result.Low.2 <- list(Required.Sample.Size=N.2, Expected.Width=w.L.M.2)
################
Conf.Limit.Desired.Certainty <- ci.R2(R2=Expected.R2(Population.R2=Population.R2, N=N.1, p=p), alpha.lower=1-degree.of.certainty, alpha.upper=0, N=N.1, p=p, Random.Predictors=Random.Predictors)$Lower.Conf.Limit.R2
N.3 <- N.1
CI.2 <- ci.R2(R2=Expected.R2(Population.R2=Conf.Limit.Desired.Certainty, N=N.3, p=p), alpha.lower=alpha.lower, alpha.upper=alpha.upper, N=N.3, p=p, Random.Predictors=Random.Predictors)
w.L <- Conf.Limit.Desired.Certainty - CI.2$Lower.Conf.Limit.R2
Diff <- w.L - width

while(Diff > Tol)
{
N.3 <- N.3 + 1
CI.2 <- ci.R2(R2=Expected.R2(Population.R2=Conf.Limit.Desired.Certainty, N=N.3, p=p), alpha.lower=alpha.lower, alpha.upper=alpha.upper, N=N.3, p=p, Random.Predictors=Random.Predictors)
w.L <- Conf.Limit.Desired.Certainty - CI.2$Lower.Conf.Limit.R2
Diff <- w.L - width
}

w.L.M.3 <- Expected.R2(Population.R2=Population.R2, N=N.3, p=p) - ci.R2(R2=Expected.R2(Population.R2=Population.R2, N=N.3, p=p), alpha.lower=alpha.lower, alpha.upper=alpha.upper, N=N.3, p=p, Random.Predictors=Random.Predictors)$Lower.Conf.Limit.R2
Result.Low.3 <- list(Required.Sample.Size=N.3, Expected.Width=w.L.M.3)

if(N.2 >= N.3) Result.Low <- Result.Low.2 # This is to see which one sided bound yields the largest sample size.
if(N.3 > N.2) Result.Low <- Result.Low.3
###
print("This sample size should be regarded as approximate; the methods have only been throughly evaluated for Full widths. You should consider evaluating the sample size with \'ss.aipe.R2.sensitivity\'.")
return(Result.Low)
}
}

if(which.width=="Upper")
{
N.1 <- N.0 + 1
CI.1 <- ci.R2(R2=Expected.R2(Population.R2=Population.R2, N=N.1, p=p), alpha.lower=alpha.lower, alpha.upper=alpha.upper, N=N.1, p=p, Random.Predictors=FALSE)

w.U <- CI.1$Upper.Conf.Limit.R2 - Expected.R2(Population.R2=Population.R2, N=N.1, p=p)
Diff <- w.U - width
while(Diff > Tol | is.na(Diff))
{
N.1 <- N.1 + 1
CI.1 <- ci.R2(R2=Expected.R2(Population.R2=Population.R2, N=N.1, p=p), alpha.lower=alpha.lower, alpha.upper=alpha.upper, N=N.1, p=p, Random.Predictors=FALSE)
w.U <- CI.1$Upper.Conf.Limit.R2 - Expected.R2(Population.R2=Population.R2, N=N.1, p=p)
Diff <- w.U - width
}

if(Random.Predictors==TRUE)
{
# Use N.1 from above (which was based on fixed predictors).
CI.1 <- ci.R2(R2=Expected.R2(Population.R2=Population.R2, N=N.1, p=p), alpha.lower=alpha.lower, alpha.upper=alpha.upper, N=N.1, p=p, Random.Predictors=TRUE)

w.U <- CI.1$Upper.Conf.Limit.R2 - Expected.R2(Population.R2=Population.R2, N=N.1, p=p)
Diff <- w.U - width
while(Diff > Tol | is.na(Diff))
{
N.1 <- N.1 + 1
CI.1 <- ci.R2(R2=Expected.R2(Population.R2=Population.R2, N=N.1, p=p), alpha.lower=alpha.lower, alpha.upper=alpha.upper, N=N.1, p=p, Random.Predictors=TRUE)
w.U <- CI.1$Upper.Conf.Limit.R2 - Expected.R2(Population.R2=Population.R2, N=N.1, p=p)
Diff <- w.U - width
}
}

if(is.null(degree.of.certainty)) 
{
Result.Up <- list(Required.Sample.Size=N.1, Expected.Width=w.U)
print("This sample size should be regarded as approximate; the methods have only been throughly evaluated for Full widths. You should consider evaluating the sample size with \'ss.aipe.R2.sensitivity\'.")
return(Result.Up)
}
#################
if(!is.null(degree.of.certainty)) 
{

############################################################################################################################
Conf.Limit.Desired.Certainty <- ci.R2(R2=Expected.R2(Population.R2=Population.R2, N=N.1, p=p), alpha.lower=0, alpha.upper=1-degree.of.certainty, N=N.1, p=p, Random.Predictors=Random.Predictors)$Upper.Conf.Limit.R2
N.2 <- N.1
CI.2 <- ci.R2(R2=Expected.R2(Population.R2=Conf.Limit.Desired.Certainty, N=N.2, p=p), alpha.lower=alpha.lower, alpha.upper=alpha.upper, N=N.2, p=p, Random.Predictors=Random.Predictors)
w.U <- CI.2$Upper.Conf.Limit.R2 - Conf.Limit.Desired.Certainty
Diff <- w.U - width

while(Diff > Tol)
{
N.2 <- N.2 + 1
CI.2 <- ci.R2(R2=Expected.R2(Population.R2=Conf.Limit.Desired.Certainty, N=N.2, p=p), alpha.lower=alpha.lower, alpha.upper=alpha.upper, N=N.2, p=p, Random.Predictors=Random.Predictors)
w.U <- CI.2$Upper.Conf.Limit.R2 - Conf.Limit.Desired.Certainty
Diff <- w.U - width
}
w.U.M.2 <- ci.R2(R2=Expected.R2(Population.R2=Population.R2, N=N.2, p=p, Random.Predictors=Random.Predictors), alpha.lower=alpha.lower, alpha.upper=alpha.upper, N=N.2, p=p)$Upper.Conf.Limit.R2 - Expected.R2(Population.R2=Population.R2, N=N.2, p=p)
Result.Up.2 <- list(Required.Sample.Size=N.2, Expected.Width=w.U.M.2)
#######################################################
Conf.Limit.Desired.Certainty <- ci.R2(R2=Expected.R2(Population.R2=Population.R2, N=N.1, p=p), alpha.lower=1-degree.of.certainty, alpha.upper=0, N=N.1, p=p, Random.Predictors=Random.Predictors)$Lower.Conf.Limit.R2
N.3 <- N.1
CI.2 <- ci.R2(R2=Expected.R2(Population.R2=Conf.Limit.Desired.Certainty, N=N.3, p=p), alpha.lower=alpha.lower, alpha.upper=alpha.upper, N=N.3, p=p, Random.Predictors=Random.Predictors)
w.U <- CI.2$Upper.Conf.Limit.R2 - Conf.Limit.Desired.Certainty
Diff <- w.U - width

while(Diff > Tol)
{
N.3 <- N.3 + 1
CI.2 <- ci.R2(R2=Expected.R2(Population.R2=Conf.Limit.Desired.Certainty, N=N.3, p=p), alpha.lower=alpha.lower, alpha.upper=alpha.upper, N=N.3, p=p, Random.Predictors=Random.Predictors)
w.U <- CI.2$Upper.Conf.Limit.R2 - Conf.Limit.Desired.Certainty
Diff <- w.U - width
}
w.U.M.3 <- ci.R2(R2=Expected.R2(Population.R2=Population.R2, N=N.3, p=p), alpha.lower=alpha.lower, alpha.upper=alpha.upper, N=N.3, p=p)$Upper.Conf.Limit.R2 - Expected.R2(Population.R2=Population.R2, N=N.3, p=p, Random.Predictors=Random.Predictors)
Result.Up.3 <- list(Required.Sample.Size=N.3, Expected.Width=w.U.M.3)

#################
if(N.2 >= N.3) Result.Up <- Result.Up.2
if(N.3 > N.2) Result.Up <- Result.Up.3
print("This sample size should be regarded as approximate; the methods have only been throughly evaluated for Full widths. You should consider evaluating the sample size with \'ss.aipe.R2.sensitivity\'.")
return(Result.Up)
}
}
}
