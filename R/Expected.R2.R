"Expected.R2" <-
function(Population.R2, N, p)
{
Value <- 1 - ((N-p-1)/(N-1))*(1-Population.R2)*hyperg_2F1(1, 1, .5*(N+1), Population.R2)
Value <- max(0, Value)
Value <- min(Value, 1)
return(Value)
}
