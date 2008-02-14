`ss.aipe.rmsea` <-
function(RMSEA, df, width, conf.level=.95)
{

if(conf.level > 50 & conf.level<100) conf.level=conf.level/100
if(conf.level< .50 | conf.level>.9999) stop("The value of 'conf.level' must be between .5 and .999")

######################################################################
n.0 <- 2

F.ML <- df*RMSEA^2

F.hat.ML.0 <- F.ML + df/(n.0-1)
Expected.RMSEA.0 <- sqrt(F.hat.ML.0/df - 1/n.0)

ci.0 <- ci.rmsea(rmsea=Expected.RMSEA.0, df=df, N=n.0, conf.level=conf.level)
w.0 <- ci.0$Upper.Conf.Limit - ci.0$Lower.Conf.Limit

n.i <- n.0
w.i <- w.0
while(w.i > width)
{
n.i <- n.i + 1
F.hat.ML.i <- F.ML + df/(n.i-1)
Expected.RMSEA.i <- sqrt(F.hat.ML.i/df - 1/n.i)

ci.i <- ci.rmsea(rmsea=Expected.RMSEA.i, df=df, N=n.i, conf.level=conf.level)
w.i <- ci.i$Upper.Conf.Limit - ci.i$Lower.Conf.Limit
}


print(paste("Necessary sample size so that the expected width of the ", conf.level*100, "% confidence interval is no greater than ", width, ", given a population RMSEA of ", RMSEA, ", is:", sep=""))
return(n.i)
}

