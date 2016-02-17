\name{upsilon}
\alias{upsilon}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function for estimating the mediation effect size as discussed in Lachowicz, Preacher, & Kelley (submitted). 
}
\description{
This function implements the upsilon effect size statistic for general mediation from Lachowicz, Preacher, & Kelley (submitted).  
}
\usage{
upsilon(data, x, m, y, covs = NULL, parallel.med = NULL, seq.med = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{\code{data} is the data that contains the variable that are to be used in the mediation model.}
  \item{x}{The \code{x} is the independent variable or variables (listed in a vector) that are included in \code{data}. That is, this is a string identifying the \code{x} variables. 
}
  \item{m}{
The \code{m} is the mediator variable or variables (listed in a vector) that are included in \code{data}. That is, this is a string identifying the \code{m} variables. 
}
  \item{y}{
The \code{y} is the outcome or dependent variable or variables (listed in a vector) that are included in \code{data}. That is, this is a string identifying the \code{y} variables. 
}
  \item{covs}{
Covariates to include in the mediation model are identified in \code{covs}. That is, this is a string identifying the \code{covs} (covariates) to inlucde. 
}
  \item{parallel.med}{
Identifies which of the mediators are parallel mediators. Not that the parallel mediators listed in \code{parallel.med} are also included in \code{m}.
}
  \item{seq.med}{
Identifies which of the mediators are sequential and in the order in a \code{list}.  
That is, \code{seq.med} identifying the sequential mediators and their order. The order is such that \code{list('m2~m1','m3~m2')} \code{m2} is caused by \code{m1} and \code{m3} is caused by \code{m2}.
}
}

\details{
See the examples below for example applications of the function. 
}
\value{
\item{IE}{is the value of the indirect effect for the path noted. Note that the total, specific, and unconditional paths will be produced (depending on the model)}
\item{Upsilon}{is the value of the effect size for the path noted. Note that the total, specific, unconditional and unique paths will be produced (depending on the model)}
%% ...
}
\references{
Lachowicz, M. J., Preacher, K. J., & Kelley, K. (submitted). A novel measure of effect size for mediation analysis. Submited for publication. 

Preacher, K. J., & Kelley, K. (2011). Effect size measures for mediation models: quantitative strategies for communicating indirect effects. \emph{Psychological Methods}, \emph{16}, 93--115. 

Wen, Z., & Fan, X. (2015). Monotonicity of effect sizes: Questioning kappa-squared as mediation effect size measure. 
\emph{Psychological Methods}, \emph{20}, 193--203. 

}
\author{
Lachowicz
Mark J. Lachowicz (Vanderbilt University; \email{Mark.J.Lachowicz@Vanderbilt.edu})
}
\note{
Note that this function overcomes some limitations of other effects for mediation models, such as those discussed in Preacher and Kelley (2012) and Wen and Fan (2015).
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\link{mediation}
}
\examples{

# To generate the multivariate data for the examples. 
require(MASS)

# Generate data for example 1 and 2. 
X<-matrix(c(1,.4,.2,.3,
            .4,1,.4,.1,
            .2,.4,1,.3,
            .3,.1,.3,1),4,4,byrow=TRUE)
data<-mvrnorm(500,c(0,0,0,0),X,empirical=TRUE)
colnames(data)<-c('x','m1','m2','y')

data <- as.data.frame(data)

#3 Example 1; three variable mediation; simple mediation model. 
x<-'x'
m<-'m1'
y<-'y'

upsilon(data,x,m,y)

## Example 2; four variable mediation; 2 mediators 
data<-mvrnorm(100,c(0,0,0,0),X,empirical=TRUE)
colnames(data)<-c('x1','x2','m','y')
data <- as.data.frame(data)

x <- c('x1','x2')
m <- 'm'
y <- 'y'

B <- upsilon(data,x,m,y)
B


# Generate data for example 3.
X1<-matrix(c(1,.4,.3,.3,.1,
             .4,1,-.2,.1,.2,
             .3,-.2,1,-.1,.3,
             .3,.1,-.1,1,.4,
             .1,.2,.3,.4,1),5,5,byrow=TRUE)

data<-mvrnorm(100,c(0,0,0,0,0),X1,empirical=TRUE)
colnames(data)<-c('x1','x2','x3','m','y')
data <- as.data.frame(data)

## Example 3; 3 predictors, 1 mediator, 1 outcome 
x <- c('x1','x2','x3')
m <- 'm'
y <- 'y'

C<-upsilon(data,x,m,y)
C

# Generate data for example 4.
X1<-matrix(c(1,.2,.3,.3,.1,
             .2,1,-.2,.3,.2,
             .3,-.2,1,.3,.3,
             .3,.3,.3,1,.4,
             .1,.2,.3,.4,1),5,5,byrow=TRUE)

data<-mvrnorm(100,c(0,0,0,0,0),X1,empirical=TRUE)
colnames(data)<-c('x1','x2','x3','m','y')
data <- as.data.frame(data)

## Example 4; 2 predictors, 1 mediator, 1 outcome, 1 covariate 
x <- c('x1','x2')
m <- 'm'
y <- 'y'
cov <-'x3'

upsilon(data,x,m,y,covs=cov)

# Generate data for example 5.
data<-mvrnorm(100,c(0,0,0,0),X,empirical=TRUE)
colnames(data)<-c('x','m1','m2','y')
data<-as.data.frame(data)

## Example 5; 1 predictor, 1 outcome, 2 parallel mediators 
x <- 'x'
m <- c('m1','m2')
y <- 'y'
meds <- c('m1+m2')

upsilon(data,x,m,y,parallel.med=meds)

# Generate data for example 6. 
data<-mvrnorm(100,c(0,0,0,0,0),X1,empirical=TRUE)
colnames(data)<-c('x','m1','m2','m3','y')
data<-as.data.frame(data)

## Example 6, 1 predictor, 1 outcome, 3 parallel mediators 
x <- 'x'
m <- c('m1','m2','m3')
y <- 'y'
meds <- c('m1+m2+m3')

upsilon(data,x,m,y,parallel.med=meds)

# Generate data for example 7.
data<-mvrnorm(100,c(0,0,0,0,0),X1,empirical=TRUE)
colnames(data)<-c('x1','x2','m1','m2','y')
data <- as.data.frame(data)

# Example 7; 2 predictors, 1 outcome, 2 parallel mediators 
x <- c('x1','x2')
m <- c('m1','m2')
y <- 'y'
parallel.med <- c('m1+m2')

upsilon(data,x,m,y,parallel.med=parallel.med)

# Generate data for example 8. 
data<-mvrnorm(100,c(0,0,0,0),X,empirical=TRUE)
colnames(data)<-c('x','m1','m2','y')
data<-as.data.frame(data)

## Example 8; 2 mediators, serial 
x <- 'x'
m <- c('m1','m2')
y <- 'y'
seq.med <- c('m2~m1')

upsilon(data,x,m,y,seq.med=seq.med)

# Generate data for example 9.
X1<-matrix(c(1,.4,.3,.3,.1,
             .4,1,-.1,.3,.2,
             .3,-.1,1,.3,.3,
             .3,.3,.3,1,.4,
             .1,.2,.3,.4,1),5,5,byrow=TRUE)

data<-mvrnorm(100,c(0,0,0,0,0),X1,empirical=TRUE)
colnames(data)<-c('x','m1','m2','m3','y')
data <- as.data.frame(data)

# Example 9; 3 sequential mediators
x <- c('x')
m <- c('m1','m2','m3')
y <- c('y')
seq.med <- list('m2~m1','m3~m2')

upsilon(data,x,m,y,seq.med=seq.med)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{design}
