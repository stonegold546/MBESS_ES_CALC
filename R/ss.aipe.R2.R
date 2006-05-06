ss.aipe.R2 <-
function(Population.R2=NULL, conf.level=.95, width=NULL, which.width="Full", p=NULL, degree.of.certainty=NULL, Tol=1e-9, ...)
{
char.expand(which.width, c("Full", "Lower", "Upper"), nomatch = stop("Problems with 'which.width' specification. You must choose either 'Full', 'Lower', or 'Upper'.", call.=FALSE))

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

Continue <- TRUE
while(Continue==TRUE)
{
CI.0 <- ci.R2(R2=Expected.R2(Population.R2=Population.R2, N=N.0, p=p), alpha.lower=alpha.lower, alpha.upper=alpha.upper, N=N.0, p=p, ...)
Continue <- sum(is.na(CI.0))>0
N.0 <- N.0 + 1
}

######################################################################
if(which.width=="Full")
{
N.1 <- N.0 + 1
CI.1 <- ci.R2(R2=Expected.R2(Population.R2=Population.R2, N=N.1, p=p), alpha.lower=alpha.lower, alpha.upper=alpha.upper, N=N.1, p=p)
w.F <- CI.1$Upper.Conf.Limit.R2 - CI.1$Lower.Conf.Limit.R2
Diff <- w.F - width
while(Diff > Tol | is.na(Diff))
{
N.1 <- N.1 + 1
CI.1 <- ci.R2(R2=Expected.R2(Population.R2=Population.R2, N=N.1, p=p), alpha.lower=alpha.lower, alpha.upper=alpha.upper, N=N.1, p=p)
w.F <- CI.1$Upper.Conf.Limit.R2 - CI.1$Lower.Conf.Limit.R2
Diff <- w.F - width
}

if(is.null(degree.of.certainty)) 
{
Result.Full <- list(Required.Sample.Size=N.1, Expected.Width=w.F)
return(Result.Full)
}

if(!is.null(degree.of.certainty)) 
{
Conf.Limit.Desired.Certainty <- ci.R2(R2=Expected.R2(Population.R2=Population.R2, N=N.1, p=p), alpha.lower=0, alpha.upper=1-degree.of.certainty, N=N.1, p=p)$Upper.Conf.Limit.R2
Conf.Limit.Desired.Certainty.Upper <- Conf.Limit.Desired.Certainty
N.2 <- N.1
CI.2 <- ci.R2(R2=Expected.R2(Population.R2=Conf.Limit.Desired.Certainty, N=N.2, p=p), alpha.lower=alpha.lower, alpha.upper=alpha.upper, N=N.2, p=p)
w.F <- CI.2$Upper.Conf.Limit.R2 - CI.2$Lower.Conf.Limit.R2
Diff <- w.F - width
while(Diff > Tol)
{
N.2 <- N.2 + 1
CI.2 <- ci.R2(R2=Expected.R2(Population.R2=Conf.Limit.Desired.Certainty, N=N.2, p=p), alpha.lower=alpha.lower, alpha.upper=alpha.upper, N=N.2, p=p)
w.F <- CI.2$Upper.Conf.Limit.R2 - CI.2$Lower.Conf.Limit.R2
Diff <- w.F - width
}
CI.2.M <- ci.R2(R2=Expected.R2(Population.R2=Population.R2, N=N.2, p=p), alpha.lower=alpha.lower, alpha.upper=alpha.upper, N=N.2, p=p)
w.F.M <- CI.2.M$Upper.Conf.Limit.R2 - CI.2.M$Lower.Conf.Limit.R2
Result.Full.2 <- list(Required.Sample.Size=N.2, Expected.Width=w.F.M)
######## 
Conf.Limit.Desired.Certainty <- ci.R2(R2=Expected.R2(Population.R2=Population.R2, N=N.1, p=p), alpha.lower=1-degree.of.certainty, alpha.upper=0, N=N.1, p=p)$Lower.Conf.Limit.R2
Conf.Limit.Desired.Certainty.Lower <- Conf.Limit.Desired.Certainty
N.3 <- N.1
CI.2 <- ci.R2(R2=Expected.R2(Population.R2=Conf.Limit.Desired.Certainty, N=N.3, p=p), alpha.lower=alpha.lower, alpha.upper=alpha.upper, N=N.3, p=p)
w.F <- CI.2$Upper.Conf.Limit.R2 - CI.2$Lower.Conf.Limit.R2
Diff <- w.F - width
while(Diff > Tol)
{
N.3 <- N.3 + 1
CI.2 <- ci.R2(R2=Expected.R2(Population.R2=Conf.Limit.Desired.Certainty, N=N.3, p=p), alpha.lower=alpha.lower, alpha.upper=alpha.upper, N=N.3, p=p)
w.F <- CI.2$Upper.Conf.Limit.R2 - CI.2$Lower.Conf.Limit.R2
Diff <- w.F - width
}
CI.3.M <- ci.R2(R2=Expected.R2(Population.R2=Population.R2, N=N.3, p=p), alpha.lower=alpha.lower, alpha.upper=alpha.upper, N=N.3, p=p)
w.F.M.3 <- CI.2.M$Upper.Conf.Limit.R2 - CI.2.M$Lower.Conf.Limit.R2
Result.Full.3 <- list(Required.Sample.Size=N.3, Expected.Width=w.F.M.3)



if(N.2 >= N.3) Result.Full <- Result.Full.2 # This is to see which one sided bound yields the largest sample size.
if(N.3 > N.2) Result.Full <- Result.Full.3

# To get maximum expected CI width to ensure adequate sample size is used.
CI.WIDTH.R2 <- function(R2, conf.level, N, p)
{
Lims <- ci.R2(R2=Expected.R2(Population.R2=R2, N, p), conf.level=conf.level, N=N, p=p)
Lims$Upper - Lims$Lower
}

Optimize.Result <- optimize(f=CI.WIDTH.R2, interval = c(0, 1), maximum = TRUE, tol = .Machine$double.eps^0.25, N=Result.Full$Required.Sample.Size, p=p, conf.level=conf.level)$maximum

if((Conf.Limit.Desired.Certainty.Lower <= Optimize.Result) & (Optimize.Result <= Conf.Limit.Desired.Certainty.Upper))
{
N.4 <- Result.Full$Required.Sample.Size
CI.4 <- ci.R2(R2=Expected.R2(Population.R2=Optimize.Result, N=N.4, p=p), alpha.lower=alpha.lower, alpha.upper=alpha.upper, N=N.4, p=p)
w.F <- CI.4$Upper.Conf.Limit.R2 - CI.4$Lower.Conf.Limit.R2
Diff <- w.F - width
while(Diff > Tol)
{
N.4 <- N.4 + 1
CI.4 <- ci.R2(R2=Expected.R2(Population.R2=Optimize.Result, N=N.4, p=p), alpha.lower=alpha.lower, alpha.upper=alpha.upper, N=N.4, p=p)
w.F <- CI.4$Upper.Conf.Limit.R2 - CI.4$Lower.Conf.Limit.R2
Diff <- w.F - width
}
CI.4.M <- ci.R2(R2=Expected.R2(Population.R2=Optimize.Result, N=N.4, p=p), alpha.lower=alpha.lower, alpha.upper=alpha.upper, N=N.4, p=p)
w.F.M.4 <- CI.4.M$Upper.Conf.Limit.R2 - CI.4.M$Lower.Conf.Limit.R2
Result.Full.4 <- list(Required.Sample.Size=N.4, Expected.Width=w.F.M.4)

if(Result.Full$Required.Sample.Size < Result.Full.4$Required.Sample.Size)
{
warning("The sample size may be overestimated to some degree. Consider using the \'ss.aipe.R2.sensitivity\' function.")
Result.Full <- Result.Full.4
} 
}



return(Result.Full)
}
}

if(which.width=="Lower")
{
N.1 <- N.0 + 1
CI.1 <- ci.R2(R2=Expected.R2(Population.R2=Population.R2, N=N.1, p=p), alpha.lower=alpha.lower, alpha.upper=alpha.upper, N=N.1, p=p)
w.L <- Expected.R2(Population.R2=Population.R2, N=N.1, p=p) - CI.1$Lower.Conf.Limit.R2
Diff <- w.L - width
while(Diff > Tol | is.na(Diff))
{
N.1 <- N.1 + 1
CI.1 <- ci.R2(R2=Expected.R2(Population.R2=Population.R2, N=N.1, p=p), alpha.lower=alpha.lower, alpha.upper=alpha.upper, N=N.1, p=p)
w.L <- Expected.R2(Population.R2=Population.R2, N=N.1, p=p) - CI.1$Lower.Conf.Limit.R2
Diff <- w.L - width
}

if(is.null(degree.of.certainty)) 
{
Result.Low <- list(Required.Sample.Size=N.1, Expected.Width=w.L)
return(Result.Low)
}

if(!is.null(degree.of.certainty)) 
{
Conf.Limit.Desired.Certainty <- ci.R2(R2=Expected.R2(Population.R2=Population.R2, N=N.1, p=p), alpha.lower=0, alpha.upper=1-degree.of.certainty, N=N.1, p=p)$Upper.Conf.Limit.R2
N.2 <- N.1
CI.2 <- ci.R2(R2=Expected.R2(Population.R2=Conf.Limit.Desired.Certainty, N=N.2, p=p), alpha.lower=alpha.lower, alpha.upper=alpha.upper, N=N.2, p=p)
w.L <- Conf.Limit.Desired.Certainty - CI.2$Lower.Conf.Limit.R2
Diff <- w.L - width

while(Diff > Tol)
{
N.2 <- N.2 + 1
CI.2 <- ci.R2(R2=Expected.R2(Population.R2=Conf.Limit.Desired.Certainty, N=N.2, p=p), alpha.lower=alpha.lower, alpha.upper=alpha.upper, N=N.2, p=p)
w.L <- Conf.Limit.Desired.Certainty - CI.2$Lower.Conf.Limit.R2
Diff <- w.L - width
}
w.L.M.2 <- Expected.R2(Population.R2=Population.R2, N=N.2, p=p) - ci.R2(R2=Expected.R2(Population.R2=Population.R2, N=N.2, p=p), alpha.lower=alpha.lower, alpha.upper=alpha.upper, N=N.2, p=p)$Lower.Conf.Limit.R2

Result.Low.2 <- list(Required.Sample.Size=N.2, Expected.Width=w.L.M.2)
################
Conf.Limit.Desired.Certainty <- ci.R2(R2=Expected.R2(Population.R2=Population.R2, N=N.1, p=p), alpha.lower=1-degree.of.certainty, alpha.upper=0, N=N.1, p=p)$Lower.Conf.Limit.R2
N.3 <- N.1
CI.2 <- ci.R2(R2=Expected.R2(Population.R2=Conf.Limit.Desired.Certainty, N=N.3, p=p), alpha.lower=alpha.lower, alpha.upper=alpha.upper, N=N.3, p=p)
w.L <- Conf.Limit.Desired.Certainty - CI.2$Lower.Conf.Limit.R2
Diff <- w.L - width

while(Diff > Tol)
{
N.3 <- N.3 + 1
CI.2 <- ci.R2(R2=Expected.R2(Population.R2=Conf.Limit.Desired.Certainty, N=N.3, p=p), alpha.lower=alpha.lower, alpha.upper=alpha.upper, N=N.3, p=p)
w.L <- Conf.Limit.Desired.Certainty - CI.2$Lower.Conf.Limit.R2
Diff <- w.L - width
}

w.L.M.3 <- Expected.R2(Population.R2=Population.R2, N=N.3, p=p) - ci.R2(R2=Expected.R2(Population.R2=Population.R2, N=N.3, p=p), alpha.lower=alpha.lower, alpha.upper=alpha.upper, N=N.3, p=p)$Lower.Conf.Limit.R2
Result.Low.3 <- list(Required.Sample.Size=N.3, Expected.Width=w.L.M.3)

if(N.2 >= N.3) Result.Low <- Result.Low.2 # This is to see which one sided bound yields the largest sample size.
if(N.3 > N.2) Result.Low <- Result.Low.3
###
return(Result.Low)
}
}

if(which.width=="Upper")
{
N.1 <- N.0 + 1
CI.1 <- ci.R2(R2=Expected.R2(Population.R2=Population.R2, N=N.1, p=p), alpha.lower=alpha.lower, alpha.upper=alpha.upper, N=N.1, p=p)

w.U <- CI.1$Upper.Conf.Limit.R2 - Expected.R2(Population.R2=Population.R2, N=N.1, p=p)
Diff <- w.U - width
while(Diff > Tol | is.na(Diff))
{
N.1 <- N.1 + 1
CI.1 <- ci.R2(R2=Expected.R2(Population.R2=Population.R2, N=N.1, p=p), alpha.lower=alpha.lower, alpha.upper=alpha.upper, N=N.1, p=p)
w.U <- CI.1$Upper.Conf.Limit.R2 - Expected.R2(Population.R2=Population.R2, N=N.1, p=p)
Diff <- w.U - width
}


if(is.null(degree.of.certainty)) 
{
Result.Up <- list(Required.Sample.Size=N.1, Expected.Width=w.U)
return(Result.Up)
}
#################
if(!is.null(degree.of.certainty)) 
{

############################################################################################################################
Conf.Limit.Desired.Certainty <- ci.R2(R2=Expected.R2(Population.R2=Population.R2, N=N.1, p=p), alpha.lower=0, alpha.upper=1-degree.of.certainty, N=N.1, p=p)$Upper.Conf.Limit.R2
N.2 <- N.1
CI.2 <- ci.R2(R2=Expected.R2(Population.R2=Conf.Limit.Desired.Certainty, N=N.2, p=p), alpha.lower=alpha.lower, alpha.upper=alpha.upper, N=N.2, p=p)
w.U <- CI.2$Upper.Conf.Limit.R2 - Conf.Limit.Desired.Certainty
Diff <- w.U - width

while(Diff > Tol)
{
N.2 <- N.2 + 1
CI.2 <- ci.R2(R2=Expected.R2(Population.R2=Conf.Limit.Desired.Certainty, N=N.2, p=p), alpha.lower=alpha.lower, alpha.upper=alpha.upper, N=N.2, p=p)
w.U <- CI.2$Upper.Conf.Limit.R2 - Conf.Limit.Desired.Certainty
Diff <- w.U - width
}
w.U.M.2 <- ci.R2(R2=Expected.R2(Population.R2=Population.R2, N=N.2, p=p), alpha.lower=alpha.lower, alpha.upper=alpha.upper, N=N.2, p=p)$Upper.Conf.Limit.R2 - Expected.R2(Population.R2=Population.R2, N=N.2, p=p)
Result.Up.2 <- list(Required.Sample.Size=N.2, Expected.Width=w.U.M.2)
#######################################################
Conf.Limit.Desired.Certainty <- ci.R2(R2=Expected.R2(Population.R2=Population.R2, N=N.1, p=p), alpha.lower=1-degree.of.certainty, alpha.upper=0, N=N.1, p=p)$Lower.Conf.Limit.R2
N.3 <- N.1
CI.2 <- ci.R2(R2=Expected.R2(Population.R2=Conf.Limit.Desired.Certainty, N=N.3, p=p), alpha.lower=alpha.lower, alpha.upper=alpha.upper, N=N.3, p=p)
w.U <- CI.2$Upper.Conf.Limit.R2 - Conf.Limit.Desired.Certainty
Diff <- w.U - width

while(Diff > Tol)
{
N.3 <- N.3 + 1
CI.2 <- ci.R2(R2=Expected.R2(Population.R2=Conf.Limit.Desired.Certainty, N=N.3, p=p), alpha.lower=alpha.lower, alpha.upper=alpha.upper, N=N.3, p=p)
w.U <- CI.2$Upper.Conf.Limit.R2 - Conf.Limit.Desired.Certainty
Diff <- w.U - width
}
w.U.M.3 <- ci.R2(R2=Expected.R2(Population.R2=Population.R2, N=N.3, p=p), alpha.lower=alpha.lower, alpha.upper=alpha.upper, N=N.3, p=p)$Upper.Conf.Limit.R2 - Expected.R2(Population.R2=Population.R2, N=N.3, p=p)
Result.Up.3 <- list(Required.Sample.Size=N.3, Expected.Width=w.U.M.3)

#################
if(N.2 >= N.3) Result.Up <- Result.Up.2
if(N.3 > N.2) Result.Up <- Result.Up.3
return(Result.Up)
}
}
}
