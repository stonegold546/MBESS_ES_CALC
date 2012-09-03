\name{ss.aipe.reliability}
\alias{ss.aipe.reliability}

\title{ Sample Size Planning for Accuracy in Parameter Estimation for reliability coefficients. }
\description{

This function determines a necessary sample size so that the expected confidence interval width for the alpha coefficient or omega coefficient
is sufficiently narrow (when assurance=NULL) or so that the obtained confidence interval is no larger than the value specified with some desired degree of certainty (i.e., a probability that the obtained width is less than the specified width; assurance=.85). This function calculates coefficient alpha based on 
McDonald's (1999) formula for coefficient alpha, also known as Guttman-Cronbach alpha. It also uses coefficient omega from McDonald (1999). When the 'Parallel' or 'True Score' model is used, coefficient alpha is calculated. When the 'Congeneric' model is used, coefficient omega is calculated.
}
\usage{
ss.aipe.reliability(model = NULL, type = NULL, width = NULL, S = NULL, 
conf.level = 0.95, assurance = NULL, data = NULL, i = NULL, cor.est = NULL, 
lambda = NULL, psi.square = NULL, initial.iter = 500, 
final.iter = 5000, start.ss = NULL)
}

\arguments{
  \item{model}{ the type of measurement model (e.g., \code{"parallel items"}, \code{"true-score equivalent"}, or 
    \code{"congeneric model"}) for a homogeneous single common factor test }
  \item{type}{the type of method to base the formation of the confidence interval on, either the \code{"Factor Analytic"} (McDonald, 1999) or \code{"Normal Theory"} (van Zyl, Neudecker, & Nel, 2000) }
  \item{width}{ the desired full width of the confidence interval}
  \item{S}{ a symmetric covariance matrix}
  \item{conf.level}{ the desired confidence interval coverage, (i.e., 1- Type I error rate) }
  \item{assurance}{ parameter to ensure that the obtained confidence interval width is narrower than the desired width with a specificied degree of certainty}
  \item{data}{the data set that the reliability coefficient is obtained from}
  \item{i}{ number of items }
  \item{cor.est}{ the estimated inter-item correlation }
  \item{lambda}{ the vector of population factor loadings}
  \item{psi.square}{ the vector of population error variances}
  \item{initial.iter}{ the number of initial iterations or generations/replications of the simulation study within the function}
  \item{final.iter}{ the number of final iterations or generations/replications of the simulation study }
  \item{start.ss}{ the initial sample size to start the simulation at }
}
%\details{
%  ~~ If necessary, more details than the description above ~~
%}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{Required.Sample.Size}{the necessary sample size}
  \item{width }{the specified full width of the confidence interval}
  \item{specified.assurance }{the specified degree of certainty}
  \item{empirical.assurance}{the empirical assurance based on the necessary sample size returned}
  \item{final.iter}{the specified number of iterations in the simulation study}
}
\references{ McDonald, R. P. (1999). \emph{Test theory: A unified approach}. Mahwah, New Jersey: Lawrence 
Erlbaum Associates, Publishers.

van Zyl, J. M., Neudecker, H., & Nel, D. G. (2000) On the distribution of the maximum likelihood estimator 
of Cronbach's alpha. \emph{Psychometrika, 65} (3), 271-280. }

\author{Leann J. Terry (Indiana University; \email{ljterry@Indiana.Edu});
Ken Kelley (University of Notre Dame; \email{KKelley@ND.Edu})}

\note{ Not all of the items can be entered into the function to represent the population values. For example, either 'data' can be used, or 'S', or  'i', 'cor.est', and 'psi.square', or 'i', 'lambda', and 'psi.square'. With a
large number of interations ('final.iter') this function may take a while to execute. Please be patient.}

\section{Warning }{In some conditions you may receive a warning, such as "In sem.default(ram = ram, S = S, N = N, param.names = pars, var.names = vars,;Could not compute QR decomposition of Hessian. Optimization probably did not converge." 
This indicates that the model likely did not converge. In certain conditions this may occur because the model is not being fit well due to small sample size, a low number of iterations, or a poorly behaved covariance matrix.  } 

\seealso{ \code{\link{CFA.1}}; \code{\link[sem]{sem}}; \code{\link{ci.reliability}};  }
\examples{

# Pop.Mat<-rbind(c(1.0000000, 0.3813850, 0.4216370, 0.3651484, 0.4472136), 
# c(0.3813850, 1.0000000, 0.4020151, 0.3481553, 0.4264014), c(0.4216370, 
# 0.4020151, 1.0000000, 0.3849002, 0.4714045), c(0.3651484, 0.3481553, 
# 0.3849002, 1.0000000, 0.4082483), c(0.4472136, 0.4264014, 0.4714045, 
# 0.4082483, 1.0000000))

# ss.aipe.reliability (model='Parallel', type='Normal Theory', width=.1, i=6, 
# cor.est=.3, psi.square=.2, conf.level=.95, assurance=.85, initial.iter=500, 
# final.iter=5000, start.ss=NULL)

# ss.aipe.reliability (model='True Score', type='Normal Theory', width=.15, i=5, 
# cor.est=.3, psi.square=c(.2, .3, .3, .2, .3), conf.level=.95, assurance=.85, 
# initial.iter=500, final.iter=5000, start.ss=110)

#ss.aipe.reliability (model='True Score', type='Normal Theory', width=.15, 
# S=Pop.Mat, conf.level=.95, assurance=.85, initial.iter=500, final.iter=5000, 
# start.ss=NULL)

#ss.aipe.reliability (model='Congeneric', type='Factor Analytic', width=.1, i=5, 
# lambda=c(.4, .4, .3, .3, .5), psi.square=c(.2, .4, .3, .3, .2), conf.level=.95, 
# assurance=.85, initial.iter=500, final.iter=5000, start.ss=340)

}

\keyword{ design }

