ci.reliability.bs<- function(data, model="Congeneric", type="Factor Analytic", conf.level=.95, B=100, ...)
{
	
if(!require(sem)) stop("This function depends on the 'sem' package. Please install the 'sem' package first") 
if(!require(boot)) stop("This function depends on the 'boot' package. Please install the 'boot' package first")
 
    # This is an internal MBESS function not meant to be used directly. 
.bs.reliability.ci <- function(data, g, model=model, type=type, conf.level=conf.level)
        {
        ci.results<-ci.reliability(data=data[g,], model=model, type=type, conf.level=conf.level)
        return(ci.results$Estimated.reliability)
        }

    # Output of bootstrap
boot.out <- boot(data=data, statistic=.bs.reliability.ci, R=B, stype="i", model=model, type=type, conf.level=conf.level)

    # Uses the output of bootstrap to get CIs.
BCA.Output <- boot.ci(boot.out=boot.out, conf=conf.level, type ="bca")[4:5]
Percentile.Output <- boot.ci(boot.out=boot.out, conf=conf.level, type ="perc")[4:5]
        
BCA.Output <- BCA.Output$bca
Percentile.Output <- Percentile.Output$perc
        
Output <- rbind(Percentile.Output, BCA.Output)  

rownames(Output) <- c("Percentile.Method", "BCa.Method")
colnames(Output) <- c("Desired.Conf.Level", "Lower.Limit.Index", "Upper.Limit.Index", "Lower.Conf.Limit", "Upper.Conf.Limit")
        
return(Output)
} # end of ci.reliability.bs<-function()
