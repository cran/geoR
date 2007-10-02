\name{names.geodata}
\alias{names.geodata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Lists names of the key elements of a geodata object}
\description{
  Produces a list with the names of the main elements of
  \code{geodata} object: coords, data, units.m, covariate and
  realisation.
  Can be useful to list names before using \code{{subset.geodata}}.
}
\usage{
\method{names}{geodata}(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an object of the class \code{geodata}.  }
}
%\details{
%  ~~ If necessary, more details than the description above ~~
%}
\value{
  A list with
  \item{coords}{names of the coordinates in the geodata object. }
  \item{data}{name(s) of the data elements in the geodata object.  }
  \item{units.m}{returns the string \code{units.m}. }
  \item{covariates}{return the covariate(s) name(s) if present
     in the \code{geodata object}}
  \item{realisations}{returns the string \code{units.m}  if present
    in the \code{geodata object}. }
  \item{other}{other elements in the \code{geodata}  object.}
}
%\references{ ~put references to the literature/web site here ~ }
%\author{ ~~who you are~~ }
%\note{ ~~further notes~~ 
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{\code{\link{names}}, \code{\link{subset.geodata}}, \code{\link{as.geodata}}. }
\examples{
names(ca20)
}
\keyword{spatial}% at least one, from doc/KEYWORDS
\keyword{manip}% __ONLY ONE__ keyword per line
