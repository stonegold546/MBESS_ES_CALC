\name{HS.data}
\alias{HS.data}
\docType{data}
\title{ Complete Data Set of Holzinger and Swineford's (1939) Study}
\description{
The \emph{complete} data set of scores of 301 subjects in 26 tests in Holzinger and Swineford's (1939) study.
}

\usage{data(HS.data)}
\format{
  A data frame with 301 observations on the following 32 variables.
  \describe{
    \item{\code{id}}{subject's ID number}
    \item{\code{Gender}}{subject's gender }
    \item{\code{grade}}{the grade the subject is on}
    \item{\code{agey}}{the year part of the subject's age}
    \item{\code{agem}}{the month part of the subject's age}
    \item{\code{school}}{the school the subject is from}
    \item{\code{visual}}{scores on visual perception test, test 1}
    \item{\code{cubes}}{scores on cubes test, test 2}
    \item{\code{paper}}{scores on paper form board test, test 3}
    \item{\code{flags}}{scores on lozenges test, test 4}
    \item{\code{general}}{scores on general information test, test 5}
    \item{\code{paragrap}}{scores on paragraph comprehension test, test 6}
    \item{\code{sentence}}{scores on sentence completion test, test 7}
    \item{\code{wordc}}{scores on word classification test, test 8}
    \item{\code{wordm}}{scores on word meaning test, test 9}
    \item{\code{addition}}{scores on add test, test 10}
    \item{\code{code}}{scores on code test, test 11}
    \item{\code{counting}}{scores on counting groups of dots test, test 12}
    \item{\code{straight}}{scores on straight and curved capitals test, test 13}
    \item{\code{wordr}}{scores on word recognition test, test 14}
    \item{\code{numberr}}{scores on number recognition test, test 15}
    \item{\code{figurer}}{scores on figure recognition test, test 16}
    \item{\code{object}}{scores on object-number test, test 17}
    \item{\code{numberf}}{scores on number-figure test, test 18}
    \item{\code{figurew}}{scores on figure-word test, test 19}
    \item{\code{deduct}}{scores on deduction test, test 20}
    \item{\code{numeric}}{scores on numerical puzzles test, test 21}
    \item{\code{problemr}}{scores on problem reasoning test, test 22}
    \item{\code{series}}{scores on series comopletion test, test 23}
    \item{\code{arithmet}}{scores on Woody-McCall mixed fundamentals, form I test, test 24}
    \item{\code{paperrev}}{scores on additional paper form board test, test 25}
    \item{\code{flagssub}}{scores on flags test, test 26}
  }
}

\details{
Holzinger and Swineford (1939) data is a widely-cited one, but almost always only a part of the dataset,
 the scores for Grant-White School, were used. The present dataset contains the complete data of Holzinger and 
 Swineford (1939).
 
 A total number of 301 pupils from Paster School and Grant-White School participated in Holzinger and 
Swineford's (1939) study. This study consists of 26 tests, which are used to measure 
the subjects' spatial, verbal, mental speed, memory, and mathematical ability. 

The spatial tests consist of \code{visual}, \code{cubes}, \code{paper}, \code{flags}, \code{paperrev}, 
and \code{flagssub}. The test 25, paper form board test (\code{paperrev}), can be used as a subsitute for 
test 3, paper form board test (\code{paper}). The test 26, flags test (\code{flagssub}), is a possible
 substitute for test 4, lozenges test (\code{flags}).
 
The verbal tests consist of \code{general}, \code{paragrap}, \code{sentence}, \code{wordc}, and \code{wordm}.

The speed tests consist of \code{addition}, \code{code}, \code{counting}, and \code{straight}.

The memory tests consist of \code{wordr}, \code{numberr}, \code{figurer}, \code{object}, \code{numberf}, 
and \code{figurew}.

The mathematical-ability tests consis of \code{deduct}, \code{numeric}, \code{problemr}, 
\code{series}, and \code{arithmet}.

}

\source{
Holzinger, K. J. and Swineford, F. A. (1939). A study in factor
analysis: The stability of a bi-factor solution. \emph{Supplementary Education
Monographs, 48}. University of Chicago.
}

\references{
Holzinger, K. J. and Swineford, F. A. (1939). A study in factor
analysis: The stability of a bi-factor solution. \emph{Supplementary Education
Monographs, 48}. University of Chicago.
}

%\examples{

}
\keyword{datasets}
