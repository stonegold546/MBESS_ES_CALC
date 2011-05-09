ci.c.ancova<-function(Psi=NULL, means=NULL, error.var.ancova=NULL, c.weights, n, x.bar, SSwithin.x, conf.level=.95, ...)
{if (length(x.bar)!=length(c.weights) ) stop("The input 'x.bar' and 'c.weights' imply different number of groups")

if(is.null(Psi) & is.null(means) ) stop("Input either 'Psi' or 'means'")
if(!is.null(Psi) & !is.null(means) ) stop("Do not input both 'Psi' and 'means' at the same time")

if(is.null(Psi)) Psi<- sum(means*c.weights)
J<- length(c.weights)
if(length(n)==1) n<-rep(n, J)
if(length(n)>1 & length(n)!=length(c.weights)) stop("The input 'n' and 'c.weights' imply different number of groups ")
########################################################################
f.x.numerater<- ( sum(c.weights*x.bar) )^2
f.x.denominator<- SSwithin.x
sample.size.weighted<- sum(c.weights^2 / n)

se.Psi2<- error.var.ancova*(sample.size.weighted + f.x.numerater/f.x.denominator)
se.Psi<- (sqrt(se.Psi2))

alpha<- 1-conf.level
nu<- sum(n)-J-1    
t.value<- qt(1-alpha/2, df=nu) 

list(lower.limit=Psi- t.value*se.Psi, upper.limit=Psi+ t.value*se.Psi)
}
