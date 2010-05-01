ss.aipe.rmsea.sensitivity <-
function(width, model.est, 
Sigma.est, Sigma.true=Sigma.est, N=NULL, conf.level=0.95, G=1000, ...){

cur.package<- search()
lib <- library()
if(sum(cur.package=="package:MASS") !=1 ) {
    if( sum(lib$results[,1]=="MASS")==1 ) library(MASS)
    else stop("This function depends on the 'MASS' package. Please install the 'MASS' package.")
    }

if(sum(cur.package=="package:sem") !=1 ) {
    if( sum(lib$results[,1]=="sem")==1 ) library(sem)
    else stop("This function depends on the 'sem' package. Please install the 'sem' package first.") 
    }

M.fit <- sem(model.est, Sigma.true, 810000)
rmsea <- summary(M.fit)$RMSEA[1]
df <- summary(M.fit)$df

if (is.null(N)) N <- ss.aipe.rmsea(RMSEA=rmsea, df=df, width=width, conf.level =conf.level) 
p<- dim(Sigma.true)[1]
Data <- matrix(NA, N, p)
rmsea.hat <- rep(NA, G)
CI.upper <- rep(NA, G)
CI.lower <- rep(NA, G)

for (g in 1:G){
    gc()
    p.rep <- seq(from=1, to=G, by=100)
    if (any(g==p.rep)) cat("Replication", g, "\n")
    
    Data <- mvrnorm(n = N, mu=rep(0,p), Sigma=Sigma.true)  
    S <- var(Data)
    colnames(S) <- rownames(S)<- rownames(Sigma.est)
    m.fit <- sem(model.est, S, N, gradtol=0.0001)
    if (any(is.na(m.fit$cov))){
        rmsea.hat[g]<- NA
        CI.upper[g] <- NA
        CI.lower[g] <- NA
        }
    else {
        rmsea.hat[g] <- summary(m.fit)$RMSEA[1]
        CI <- ci.rmsea(rmsea.hat[g], df=df, N=N)
        CI.lower[g] <- CI$Lower.Conf.Limit
        CI.upper[g] <- CI$Upper.Conf.Limit
        }
    }# end of for (g in 1:G)

rmsea.hat<- rmsea.hat
CI.upper <- CI.upper
CI.lower <- CI.lower
w <- CI.upper - CI.lower
suc.rep<-length(na.omit(rmsea.hat))

result <- list()
result$w <- w
result$RMSEA.hat <- rmsea.hat
result$sample.size <- N
result$df <- df
result$RMSEA.pop <- rmsea
result$desired.width <- width
result$mean.width <- mean(w, na.rm=TRUE)
result$median.width <- median(w, na.rm=TRUE)
result$assurance <- sum(w< width, na.rm=TRUE)
result$quantile.width <- quantile(w, c(.99, .97, .95, .90, .80, .70, .60), na.rm=TRUE)
result$alpha.upper <- sum(CI.upper< rmsea, na.rm=TRUE)/suc.rep
result$alpha.lower <- sum(CI.lower> rmsea, na.rm=TRUE)/suc.rep
result$alpha <- result$alpha.upper+result$alpha.lower 
result$conf.level <- conf.level
result$suc.rep <- suc.rep
return(result)
}# end of ss.aipe.rmsea.sensitiviti <- function()

