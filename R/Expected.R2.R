"Expected.R2" <-
function(Population.R2, N, p)
{
if(!require(gsl)) stop("This function depends on the 'gsl' package. Please install the 'gsl' package first")
if(!require(MASS)) stop("This function depends on the 'MASS' package. Please install the 'MASS' package first")

Value <- 1 - ((N-p-1)/(N-1))*(1-Population.R2)*hyperg_2F1(1, 1, .5*(N+1), Population.R2)
Value <- max(0, Value)
Value <- min(Value, 1)
return(Value)
}
