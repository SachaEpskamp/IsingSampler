\name{EstimateIsing}
\alias{EstimateIsing}
\alias{EstimateIsingPL}
\alias{EstimateIsingUni}
\alias{EstimateIsingBi}
\alias{EstimateIsingLL}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
non-regularized estimation methods for the Ising Model
}
\description{
This function can be used for several non-regularized estimation methods of the Ising Model. See details.
}
\usage{
EstimateIsing(data, responses, beta = 1, method = c("pl", "uni", "bi", "ll"), ...)
EstimateIsingUni(data, responses, beta = 1, ...)
EstimateIsingBi(data, responses, beta = 1, ...)
EstimateIsingPL(data, responses, beta = 1, ...)
EstimateIsingLL(data, responses, beta = 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
 Data frame with binary responses to estimate the Ising model over
}
  \item{responses}{
 Vector of length two indicating the response coding (usually \code{c(0L, 1L)} pr \code{c(-1L, 1L)})
}
  \item{beta}{
 Inverse temperature parameter
}
  \item{method}{
  The method to be used. \code{pl} uses pseudolikelihood estimation, \code{uni} sequential univariate regressions, \code{bi} bivariate regressions and \code{ll} estimates the Ising model as a loglinear model.
}
  \item{\dots}{
 Arguments sent to estimator functions
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list containing the estimation results:
\item{graph}{The estimated network}
\item{thresholds}{The estimated thresholds}
\item{results}{The results object used in the analysis}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (data, responses, beta = 1, method = c("pl", "uni", 
    "bi", "ll"), ...) 
{
    switch(method[[1]], pl = EstimateIsingPL(data, responses, 
        beta, ...), uni = EstimateIsingUni(data, responses, beta, 
        ...), bi = EstimateIsingBi(data, responses, beta, ...), 
        ll = EstimateIsingLL(data, responses, beta, ...))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
