"Expected.R2" <-
function(Population.R2, N, p)
{current.package<- search()
lib <- library()
if( sum(current.package=="package:gsl")!=1 ) {
    if( sum(lib$results[,1]=="gsl")==1 ) library(gsl)
    else stop("This function depends on the 'gsl' package. Please install the 'gsl' package 
    as you installed the 'MBESS' package")
    
    if( sum(lib$results[,1]=="MASS")==1 ) library(MASS)
    else stop("This function depends on the 'MASS' package. Please install the 'MASS' package 
    as you installed the 'MBESS' package")
    }
Value <- 1 - ((N-p-1)/(N-1))*(1-Population.R2)*hyperg_2F1(1, 1, .5*(N+1), Population.R2)
Value <- max(0, Value)
Value <- min(Value, 1)
return(Value)
}
