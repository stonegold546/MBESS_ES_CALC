conf.limits.nct <- function(ncp, df, conf.level=NULL, alpha.lower=NULL, alpha.upper=NULL, tol=1e-9, ...)
{       
tval <- abs(ncp)
###################################################################################################################
if(df <= 0) stop("The degrees of freedom must be some positive value.", call.=FALSE)
if(tol > 1e-7) warning("Your specified tolerance may not be restrictive enough (i.e., \'tol\' may be too large); your solutions may not be accurate.", call.=FALSE)
if((!is.null(alpha.lower) | !is.null(alpha.upper)) & !is.null(conf.level)) stop("You must choose either to use \'conf.level\' or define the \'lower.alpha\' and \'upper.alpha\' values, but not both", call.=FALSE)

############################This section determines the lower bound for the confidence interval####################

# To solve a tricky problem when the 'ncp' is zero and 'alpha.lower' and 'alpha.upper' are null.
if(!is.null(conf.level) & ncp==0)
{
alpha.lower <- (1-conf.level)/2
alpha.upper <- (1-conf.level)/2
conf.level <- NULL
}

if(is.null(conf.level)) ulim <- 1 - alpha.lower
if(is.null(alpha.lower) & is.null(alpha.upper)) ulim <- 1 - (1-conf.level)/2
# This first part finds a lower value from which to start.
if(tval!=0)
{
lc <- c(-tval,tval/2,tval)
    while(pt(tval,df,lc[1])<ulim)    
         {
         lc <- c(lc[1]-tval,lc[1],lc[3])
         }

# This next part finds the lower limit for the ncp.
diff <- 1
    while(diff > tol)
         {
         if(pt(tval,df,lc[2])<ulim)
         lc <- c(lc[1],(lc[1]+lc[2])/2,lc[2])
         else lc <- c(lc[2],(lc[2]+lc[3])/2,lc[3])
         if(abs(lc[2]) > 37.62) print("The noncentrality parameter of the noncentral t-distribution has exceeded 37.62 in magnitude (R's limitation for accurate probabilities from the noncentral t-distribution) in the function's iterative search for the appropriate value(s). The results are likely fine, but they also may be somewhat inaccurate.")
         diff <- abs(pt(tval,df,lc[2]) - ulim)
         ucdf <- pt(tval,df,lc[2])
         }
res.1 <- ifelse(ncp >= 0,lc[2],-lc[2])
}

############################This section determines the upper bound for the confidence interval#####################

if(is.null(conf.level)) llim <- alpha.upper
if(is.null(alpha.lower) & is.null(alpha.upper)) llim <- (1-conf.level)/2
# This first part finds an upper value from which to start.
if(tval!=0)
{
uc <- c(tval,1.5*tval,2*tval)
   while(pt(tval,df,uc[3])>llim)   
        {
        uc <- c(uc[1],uc[3],uc[3]+tval)
        }

# This next part finds the upper limit for the ncp.
diff <- 1
    while(diff > tol)         
         {
         if(pt(tval,df,uc[2])<llim) uc <- c(uc[1],(uc[1]+uc[2])/2,uc[2])
         else uc <- c(uc[2],(uc[2]+uc[3])/2,uc[3])
         if(abs(uc[3]) > 37.62) print("The noncentrality parameter of the noncentral t-distribution has exceeded 37.62 in magnitude (R's limitation for accurate probabilities from the noncentral t-distribution) in the function's iterative search for the appropriate value(s). The results are likely fine, but they also may be somewhat inaccurate.")
         diff <- abs(pt(tval,df,uc[2]) - llim)
         lcdf <- pt(tval,df,uc[2])
         }
res <- ifelse(ncp >= 0,uc[2],-uc[2])
}

#################################This part Compiles the results into a vector#######################################
# Note that 1-ucdf is used rather than ucdf because ucdf gives the probability of greater than the ncp for a 
# distribution whose noncentrality parameter is Lower.Limit. Likewise, 1-lcdf is used rather than lcdf because
# ucdf gives the probability of greater than 'ncp' for a t-distribution whose noncentrality parameter is Upper.Limit.
# We are interested in the confidence limits for ncp, thus taking 1-ucdf and 1-lcdf
if(!is.null(conf.level) | (!is.null(alpha.lower) & !is.null(alpha.upper) & tval!=0))
{
Result <- list(Lower.Limit=min(res,res.1), Prob.Less.Lower=1-ucdf, Upper.Limit=max(res,res.1), Prob.Greater.Upper=lcdf)
}

if(is.null(conf.level))
{
if(!is.null(alpha.lower) & alpha.upper!=0 & alpha.lower==0)
{
if(tval==0) Result <- list(Lower.Limit=-Inf, Prob.Less.Lower=0, Upper.Limit=qt(1-llim, df), Prob.Greater.Upper=llim)

if(tval!=0) Result <- list(Lower.Limit=-Inf, Prob.Less.Lower=0, Upper.Limit=max(res,res.1), Prob.Greater.Upper=lcdf)
}
}

if(is.null(conf.level))
{
if(!is.null(alpha.upper) & alpha.lower!=0 & alpha.upper==0)
{
if(tval==0) Result <- list(Lower.Limit=qt(1-ulim, df), Prob.Less.Lower=1-ulim, Upper.Limit=Inf, Prob.Greater.Upper=0)

if(tval!=0) Result <- list(Lower.Limit=min(res,res.1), Prob.Less.Lower=1-ucdf, Upper.Limit=Inf, Prob.Greater.Upper=0)
}
}


if(is.null(conf.level))
{
if(alpha.lower!=0 & alpha.upper!=0)
{
if(tval==0) Result <- list(Lower.Limit=qt(1-ulim, df), Prob.Less.Lower=1-ulim, Upper.Limit=qt(1-llim, df), Prob.Greater.Upper=llim)
}
}

return(Result)
}
