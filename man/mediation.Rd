\name{mediation}
\alias{mediation}

\title{ Effect sizes and confidence intervals in a mediation model }
\description{

Automate the process of simple mediation analysis (one independent ) and effect size estimation for mediation models, as discussed in Preacher and Kelley (2010). 
}
\usage{
mediation(x, mediator, dv, S = NULL, N = NULL, x.location.S = NULL, mediator.location.S = NULL, dv.location.S = NULL, mean.x = NULL, mean.m = NULL, mean.dv = NULL, conf.level = 0.95, bootstrap = FALSE, B = 1000)
}

\arguments{
  \item{x}{vector of the predictor/independent variable}
  \item{mediator}{vector of the mediator variable}
  \item{dv}{vector of the dependent/outcome variable}
  \item{S}{Covariance matrix}
  \item{N}{Sample size, necessary when a covariance matrix (\code{S}) is used}
  \item{x.location.S}{ location of the predictor/independent variable in the covariance matrix (\code{S})}
  \item{mediator.location.S}{location of the mediator variable in the covariance matrix (\code{S})}
  \item{dv.location.S}{location of the dependent/outcome variable in the covariance matrix (\code{S})}
  \item{mean.x}{mean of the \code{x} (independent/predictor) variable when a covariance matrix (\code{S}) is used}
  \item{mean.m}{mean of the \code{m} (mediator) variable when a covariance matrix (\code{S}) is used}
  \item{mean.dv}{mean of the \code{y/dv} (dependent/outcome) variable when a covariance matrix (\code{S}) is used}
  \item{conf.level}{desired level of confidence (e.g., .90, .95, .99, etc.)}
  \item{bootstrap}{TRUE or FALSE, based on whether or not a bootstrap procedure is performed to obtain confidence intervals for the various effect sizes}
  \item{B}{number of bootstrap replications when \code{bootstrap=TRUE} (e.g., 10000)}
}
\details{
Based on the work of Preacher and Kelley (2010), this function implements (simple) mediation analysis in a way that automates much of the results that are generally of interest. More specifically, three regression outputs are automated as is the calculation of effect sizes that are thought to be useful or potentially useful in the context of mediation. Much work on mediation models exists in the literature, which should be consulted for proper interpretation of the effect sizes, models, and meaning of results. 

}
\value{
  \item{Y.on.X$Regression.Table}{Regression table of \code{Y} conditional on \code{X}}
  \item{Y.on.X$Model.Fit}{Summary of model fit for the regression of \code{Y} conditional on \code{X}}
  \item{M.on.X$Regression.Table}{Regression table of \code{X} conditional on \code{M}}
  \item{M.on.X$Model.Fit}{Summary of model fit for the regression of \code{X} conditional on \code{M}}
  \item{Y.on.X.and.M$Regression.Table}{Regression table of \code{Y} conditional on \code{X} and \code{M}}
  \item{Y.on.X.and.M$Model.Fit}{Summary of model fit for the regression of \code{Y} conditional on \code{X} and \code{M}}

\item{Indirect.Effect}{the product of \eqn{\hat{a} \times \hat{b}}, where \eqn{\hat{a}} and \eqn{\hat{b}} are the estimated coefficients of the path from the independent variable to the mediator and the path from the mediator to the dependent variable}

\item{Indirect.Effect.Partially.Standardized}{It is the indirect effect (see \code{Indirect.Effect} above) divided by the estimated standard deviation of \code{Y} (MacKinnon, 2008)}

\item{Index.of.Mediation}{Index of mediation (indirect effect multiplied by the ratio of the standard deviation of x to the standard deviation of y) (Preacher and Hayes, 2008)}

\item{R2_4.5}{An index of explained variance see MacKinnon (2008) for details}

\item{R2_4.6}{An index of explained variance see MacKinnon (2008) for details}

\item{R2_4.7}{An index of explained variance see MacKinnon (2008) for details}

\item{Maximum.Possible.Mediation.Effect}{the maximum attainable value of the mediation effect (i.e., the indirect effect), in the direction of the observed indirect effect, could have been observed, conditional on the sample variances and on the magnitudes of relationships among some of the variables}

\item{ab.to.Maximum.Possible.Mediation.Effect_kappa.squared}{the proportion of the maximum possible indirect effect; Uses the indirect effect in the numerator with the maximum possible mediation effect in the denominator (Preacher & Kelley, 2010)}

\item{Mediation.Ratio}{ratio of the indirect effect to the total effect (Freedman, 2001); also known as mediation ratio (Ditlevsen, Christensen, Lynch, Damsgaard, & Keiding, 2005); in epidemiological research and as the relative indirect effect (Huang, Sivaganesan, Succop, & Goodman, 2004); often loosely interpreted as the relative indirect effect}

\item{Indirect.to.Direct}{ratio of the indirect effect to the direct effect (Sobel, 1982)}

\item{Success.of.Surrogate.Endpoint}{Success of a surrogate endpoint (Buyse and Molenberghs, 1998) }

\item{SOS}{shared over simple effects (SOS) index, which is the ratio of the variance in Y explained by both \code{X} and \code{M} divided by the variance in \code{Y} explained by \code{X} (Lindenberger and Potter, 1998)}

\item{Residual.Based_Gamma}{A residual based index (Preacher & Kelley, 2010)}

\item{Residual.Based.Standardized_gamma}{A residual based index that is standardized, where the scales of M and Y are removed by using standardized values of M and Y (Preacher & Kelley, 2010)}

\item{ES.for.two.groups}{When X is 0 and 1 representing a two group structure, 
Hansen and McNeal's (1996) Effect Size Index for Two Groups}


}
\references{ 
Buyse, M., & Molenberghs, G. (1998). Criteria for the validation of surrogate endpoints in randomized experiments. \emph{Biometrics, 54}, 1014-1029.

Ditlevsen, S., Christensen, U., Lynch, J., Damsgaard, M. T., & Keiding, N. (2005). The mediation proportion: A structural equation approach for estimating the proportion of exposure effect on outcome explained by an intermediate variable. \emph{Epidemiology, 16}, 114-120.

Freedman, L. S. (2001). Confidence intervals and statistical power of the 'Validation' ratio for surrogate or intermediate endpoints. \emph{Journal of Statistical Planning and Inference, 96}, 143-153.

Hansen, W. B., & McNeal, R. B. (1996). The law of maximum expected potential effect: Constraints placed on program effectiveness by mediator relationships. \emph{Health Education Research, 11}, 501-507.

Huang, B., Sivaganesan, S., Succop, P., & Goodman, E. (2004). Statistical assessment of mediational effects for logistic mediational models. \emph{Statistics in Medicine, 23}, 2713-2728.

Lindenberger, U., & Potter, U. (1998). The complex nature of unique and shared effects in hierarchical linear regression: Implications for developmental psychology. \emph{Psychological Methods, 3}, 218-230.

MacKinnon, D. P. (2008). \emph{Introduction to statistical mediation analysis}. Mahwah, NJ: Erlbaum.

Preacher, K. J., & Hayes, A. F. (2008b). Asymptotic and resampling strategies for assessing and comparing indirect effects in multiple mediator models. \emph{Behavior Research Methods, 40}, 879-891.

Preacher, K. J., & Kelley, K. (2010). Effect size measures for mediation models: Quantitative and graphical strategies for communicating indirect effects. \emph{Manuscript in preparation}.

Sobel, M. E. (1982). Asymptotic confidence intervals for indirect effects in structural equation models. In S. Leinhardt (Ed.), \emph{Sociological Methodology 1982} (pp. 290-312). Washington DC: American Sociological Association.
}

\author{Ken Kelley (University of Notre Dame; KKelley@nd.edu)}
%\note{ }
\seealso{\code{\link{mediation.effect.plot}}, \code{\link{mediation.effect.bar.plot}}}
\examples{
# Using the Jessor data discussed in Preacher and Kelley (2010), to illustrate the methods based on summary statistics. 
 
# mediation(S=rbind(c(2.26831107,  0.6615415, -0.08691755), 
#  c(0.66154147,  2.2763549, -0.22593820), 
#  c(-0.08691755, -0.2259382,  0.09218055)), 
# N=432, x.location.S=1, mediator.location.S=2, dv.location.S=3, mean.x=7.157645, mean.m=5.892785, mean.dv=1.649316, conf.level=.95)

}

\keyword{multivariate}
