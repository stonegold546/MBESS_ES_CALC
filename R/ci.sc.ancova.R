`ci.sc.ancova` <-
function(Psi=NULL, means=NULL, s.anova=NULL, s.ancova, standardizer="s.ancova", c.weights, n, x.bar, SSwithin.x, conf.level=.95)
{options(warn=-1)
if(standardizer=="s.anova")
    {if(is.null(s.anova)) stop("'s.anova' is needed to standardized the contrast")}
if(standardizer!="s.ancova" & standardizer!="s.anova") stop("The standardizer must be either 's.anova' or 's.ancova'")

if(is.null(Psi) & is.null(means) ) stop("Input either 'Psi' or 'means'")
if(!is.null(Psi) & !is.null(means) ) stop("Do not input both 'Psi' and 'means'")

if(is.null(Psi)) Psi<- sum(means*c.weights)

if(sum(c.weights)!=0) stop("The sum of the coefficients must be zero")
if(sum(c.weights[c.weights>0])>1) stop("Please use fractions to specify the contrast weights")

J<-length(c.weights)
if(length(n)>1)
    {if (length(n)!=J) stop("'c.weights' and 'n' imply different numbers of groups.")}
if(length(n)==1) n<-rep(n, J)

if(length(x.bar)!=J) stop("'c.weights' and 'x.bar' imply different numbers of groups.")

f.x.numerater<- (sum(c.weights*x.bar))^2
f.x.denominator<- SSwithin.x
sample.size.weighted<- sum(c.weights^2 / n )
ratio<- s.ancova/s.anova
alpha<-1-conf.level
nu<-sum(n)-J-1

if(standardizer=="s.ancova")
    {psi<- Psi/s.ancova
    lambda.obs<- psi/sqrt(sample.size.weighted+(f.x.numerater/f.x.denominator))
    lambda.limits<-conf.limits.nct(ncp=lambda.obs, df=nu, conf.level=1-alpha)
    
    psi.limit.upper<- lambda.limits$Upper.Limit*sqrt(sample.size.weighted+(f.x.numerater/f.x.denominator))
    psi.limit.lower<- lambda.limits$Lower.Limit*sqrt(sample.size.weighted+(f.x.numerater/f.x.denominator))
    }

if(standardizer=="s.anova")
    {psi<-Psi/s.anova
    lambda.obs<-psi/(ratio*sqrt(sample.size.weighted+(f.x.numerater/f.x.denominator)))
    lambda.limits<-conf.limits.nct(ncp=lambda.obs, df=nu, conf.level=1-alpha)
    
    psi.limit.upper<- lambda.limits$Upper.Limit*ratio*sqrt(sample.size.weighted+(f.x.numerater/f.x.denominator))
    psi.limit.lower<- lambda.limits$Lower.Limit*ratio*sqrt(sample.size.weighted+(f.x.numerater/f.x.denominator))
    }

list(stadardizer=standardizer, psi.limit.lower=psi.limit.lower, psi.limit.upper=psi.limit.upper)
}

