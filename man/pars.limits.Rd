\name{pars.limits}
\alias{pars.limits}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Set limits for the parameter values}
\description{
  The functions \code{\link{likfit}} and \code{\link{variofit}} in the
  package \pkg{geoR}
}
\usage{
pars.limits(phi = c(lower = 0, upper = +Inf),
            sigmasq = c(lower = 0, upper = +Inf),
            nugget.rel = c(lower = 0, upper = +Inf),
            kappa = c(lower = 0, upper = +Inf),
            kappa2 = c(lower = 0, upper = +Inf), 
            lambda = c(lower = -3, upper = 3),
            psiR = c(lower = 1, upper = +Inf),
            psiA = c(lower = 0, upper = 2 * pi),
            tausq.rel = nugget.rel)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phi}{a two elements vector with limits for the parameter phi. Defaults to [0, +Inf]}
  \item{sigmasq}{idem for sigmasq. Defaults to [0, +Inf]}
  \item{nugget.rel}{idem for nugget.rel. Defaults to [0, +Inf]}
  \item{kappa,kappa2}{idem. Defaults to [0, +Inf]}
  \item{lambda}{idem for lambda. Defaults to [-3, +3]. Only used in
    \code{\link{likfit}}. }
  \item{psiR}{idem for psiR. Defaults to [1, +Inf]. Only used in
    \code{\link{likfit}}. }
  \item{psiA}{idem for psiA. Defaults to [0, 2 pi]. Only used in
    \code{\link{likfit}}. }
  \item{tausq.rel}{idem for tausq.rel. Defaults to [0, +Inf]}
}
\details{
  Lower and upper limits for parameter values can be
  individually specified.
  For example, including the following in the function call in
  \code{likfit} or \code{variofit}:\cr
  \code{limits = pars.limits(phi=c(0, 10), lambda=c(-2.5, 2.5))}, \cr
  will change the limits for the parameters \eqn{\phi}{phi} and \eqn{\lambda}{lambda}.
  Default values are used if the argument \code{limits} is not provided.  
  }
\value{
  A list of a 2 elements vector with limits for each parameters
}
%\references{ ~put references to the literature/web site here ~ }
%\author{ ~~who you are~~ }
%\note{ ~~further notes~~ }

% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{likfit}}, \code{\link{variofit}}}

\examples{
pars.limits(phi=c(0,10))
pars.limits(phi=c(0,10), sigmasq=c(0, 100))
}
\keyword{spatial}% at least one, from doc/KEYWORDS
\keyword{models}% __ONLY ONE__ keyword per line
