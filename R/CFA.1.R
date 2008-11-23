CFA.1 <- function(S, N, equal.loading=FALSE, equal.error=FALSE) 
{current.package<- search()
lib <- library()
if( sum(current.package=="package:sem")!=1 ) {
    if( sum(lib$results[,1]=="sem")==1 ) library(sem)
    else stop("This function depends on the 'sem' package. Please install the 'sem' package 
    as you installed the 'MBESS' package")
    }

if(!isSymmetric(S, tol=1e-5)) stop ("Input a symmetric covariance or correlation matrix 'S'")

q<- nrow(S)
x<- matrix(NA, nrow=q, ncol=1)
x<- paste("x", row(x), sep="")

if(equal.loading)
    {lamda <- matrix(rep("lamda", q), nrow=q, ncol=1)    
    }
else 
    {lamda<-matrix(NA, nrow=q, ncol=1)
    lamda<-paste("lamda", row(lamda), sep="")
    }

if(equal.error)
    {psi.sq<- matrix(rep("psi.sq", q), nrow=q, ncol=1 )
    }
else
    {psi.sq<- matrix(NA, nrow=q, ncol=1)
    psi.sq<- paste("psi.sq", row(psi.sq), sep="")
    }

model.1<-cbind(paste("ksi", "->", x), lamda)
model.2<-cbind(paste(x, "<->", x), psi.sq)
model<-rbind(model.1, model.2)
start<- matrix(NA, nrow=nrow(model), ncol=1)
model<- cbind(model, start)
model<-rbind(model, c(paste("ksi", "<->", "ksi"), NA, 1))
class(model)<-"mod"

rownames(S)<- x
colnames(S)<- x
model.fitted<- sem(model, N=N, S=S) 

if (equal.loading) k<- 1 else k<-q
Factor.Loadings<- model.fitted$coeff[1:k]
if (equal.error) m<-1 else m<-q
Indicator.var<- model.fitted$coeff[k+1:m]
Parameter.cov<- model.fitted$cov

result<-list(Model=model, Factor.Loadings=Factor.Loadings, Indicator.var=Indicator.var, Parameter.cov=Parameter.cov)
return(result)
}
