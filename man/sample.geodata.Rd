\name{sample.geodata}
\alias{sample.geodata}

\title{Sampling from geodata objects}
\description{
  This functions facilitates extracting samples from geodata objects.
}
\usage{
sample.geodata(x, size, replace = FALSE, prob = NULL, coef.logCox,
               external)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an object of the class \code{geodata}. }
  \item{size}{non-negative integer giving the number of items to choose.}
  \item{replace}{Should sampling be with replacement?}
  \item{prob}{A vector of probability weights for obtaining the elements of
          the data points being sampled.}
  \item{coef.logCox}{optional. A scalar with the coeficient for the
    log-Cox process. See DETAILS below.  }
  \item{external}{numeric values of a random field to be used
    in the log-Cox inhomogeneous poisson process. }
}
\details{
  If \code{prob=NULL} and 
  the argument \code{coef.logCox}, is provided,
  sampling follows a log-Cox proccess, i.e.
  the probability of each point being sampled is proportional to:
  \deqn{exp(b Y(x))}{exp(b * Y(x))}
  with \eqn{b}{b} given by the value passed to the argument
  \code{coef.logCox} and \eqn{Y(x)}{Y(x)} taking values passed to
  the argument \code{external} or, if this is missing,
  the element \code{data} of the \code{geodata} object.
  Therefore, the latter generates a preferential sampling. 
}

\value{
  a list which is an object of the class \code{geodata}.
}

%\references{ ~put references to the literature/web site here ~ }
%\author{ ~~who you are~~ }
%\note{ ~~further notes~~ 
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}

\seealso{\code{\link{as.geodata}}, \code{\link{sample}}.}
\examples{
par(mfrow=c(1,2))
S1 <- grf(2500,  grid="reg", cov.pars=c(1, .23))
image(S1, col=gray(seq(0.9,0.1,l=100)))
y1 <- sample.geodata(S1, 80)
points(y1$coords, pch=19)
## Now a preferential sampling
y2 <- sample.geodata(S1, 80, coef=1.3)
## which is equivalent topps
## y2 <- sample.geodata(S1, 80, prob=exp(1.3*S1$data))
points(y2$coords, pch=19, col=2)
## and now a clustered (but not preferential)
S2 <- grf(2500,  grid="reg", cov.pars=c(1, .23))
y3 <- sample.geodata(S1, 80, prob=exp(1.3*S2$data))
## which is equivalent to
## points(y3$coords, pch=19, col=4)
image(S2, col=gray(seq(0.9,0.1,l=100)))
points(y3$coords, pch=19, col=4)
}
\keyword{spatial}
\keyword{manip}
