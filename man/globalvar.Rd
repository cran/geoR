\name{globalvar}
\alias{globalvar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
globalvar(geodata, locations, coords = geodata$coords, krige)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{geodata}{ ~~Describe \code{geodata} here~~ }
  \item{locations}{ }
  \item{coords}{ ~~Describe \code{coords} here~~ }
  \item{krige}{a list defining the model components and the type of
    kriging. It can take an output to a call to \code{krige.control} or
    a list with elements as for the arguments in \code{krige.control}.
  More details in the documention for  }
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{Isaaks, E.S and Srivastava, R.M. (1989) An Introduction to
  Applied Geostatistics, pag. 508, eq. 20.7. Oxford University Press. }
\author{ ~~who you are~~ }
%\note{ ~~further notes~~ 
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{\code{\link{krige.conv}} for the kriging algorithm. }
\examples{
data(s100)
}
\keyword{spatial}% at least one, from doc/KEYWORDS
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
