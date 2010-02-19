`ss.aipe.sem.path.sensitiv` <-
function(est.model, est.Sigma, true.Sigma=est.Sigma, which.path,
desired.width, conf.level=.95, assurance=NULL, G=1000, ...){
cur.package<- search()
lib <- library()
if(sum(cur.package=="package:MASS") !=1 ) {
    if( sum(lib$results[,1]=="MASS")==1 ) library(MASS)
    else stop("This function depends on the 'MASS' package. Please install the 'MASS' package.")
    }

if(sum(cur.package=="package:sem") !=1 ) {
    if( sum(lib$results[,1]=="sem")==1 ) library(sem)
    else stop("This function depends on the 'sem' package. Please install the 'sem' package.") 
    }

result.plan <- ss.aipe.sem.path(model=est.model, Sigma=est.Sigma, desired.width=desired.width, 
which.path=which.path, conf.level=conf.level, assurance=assurance, ...)

obs.vars<- result.plan$obs.vars
J <- length(result.plan$parameter)
alpha <- 1-conf.level
p <- dim(est.Sigma)[1]
N<- result.plan$sample.size
j <- result.plan$path.index

Data<- matrix(NA, N,p)
S<- matrix(NA, p,p)
theta.hat.j <- rep(NA, G)
#theta.hat <- array(NA, c(J,1,G))
SE.theta.hat.j <- rep(NA, G)
cov.theta.hat <- matrix(NA, J,J)

for(g in 1:G){
gc()
    p.rep <- seq(from=1, to=G, by=100)
    if (any(g==p.rep)) cat("Replication", g, "\n")
    #singular<- TRUE
    #while(singular){
        Data <- mvrnorm(n = N, mu=rep(0,p), Sigma=true.Sigma)  
        S <- var(Data)
    #    singular <- ifelse (determinant(S, logarithm = FALSE)$modulus <.04, TRUE, FALSE)
    #    }
    
    colnames(S) <- rownames(S)<- obs.vars
    S.fit <- sem(ram=est.model, S=S, N=N, gradtol=0.0001, ...)    
    
    cov.theta.hat <- S.fit$cov
    SE.theta.hat.j[g] <- ifelse (any(diag(cov.theta.hat)< 0), NA, sqrt(cov.theta.hat[j,j]))

    #theta.hat[,,g] <- S.fit$coeff
    theta.hat.j[g] <- S.fit$coeff[j]
    }# end of for(g in 1:G)
CI.upper <- theta.hat.j+ qnorm(1-alpha/2)*SE.theta.hat.j
CI.lower <- theta.hat.j- qnorm(1-alpha/2)*SE.theta.hat.j
w <- CI.upper - CI.lower

true.fit <- sem(ram=est.model, S=true.Sigma, N=50000)
theta.j <- true.fit$coeff[j]

CI.upper <- na.omit(CI.upper)
CI.lower <- na.omit(CI.lower)
w <- na.omit(w)
G<- length(w)

percent.narrower <- sum(w<=desired.width)/G
alpha.emp.upper <- sum(theta.j>CI.upper)/G
alpha.emp.lower <- sum(theta.j<CI.lower)/G

result<- list()
#result$SE.theta.j.hat <- SE.theta.hat.j
result$w <- w
#result$theta.j.hat <- theta.hat.j
result$sample.size <- N
#result$var.theta.j <- result.plan$var.theta.j
result$path.of.interest <- which.path
result$desired.width <- desired.width
result$mean.width <- mean(w)
result$median.width <- median(w)
result$quantile.width <- quantile(w, c(.99,.95,.90,.85,.80, .75, .70, .60))
result$width.less.than.desired <- percent.narrower
result$Type.I.err.upper <- alpha.emp.upper
result$Type.I.err.lower <- alpha.emp.lower
result$Type.I.err <- alpha.emp.upper+alpha.emp.lower
result$conf.level <- conf.level
result$rep <- G
#result$CI.upper<- CI.upper
#result$CI.lower<- CI.upper
#result$w <- w
return(result)  
}# end of function()

