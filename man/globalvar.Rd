\name{globalvar}
\alias{globalvar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Computes global variance }
\description{
  Global variance computation for a set of locations using the covarianve
  model
}
\usage{
globalvar(geodata, locations, coords = geodata$coords, krige)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{geodata}{an object of the class \code{geodata} }
  \item{locations}{n by 2 matrix with a set of locations, typically a
    prediction grid }
  \item{coords}{data coordinates }
  \item{krige}{a list defining the model components and the type of
    kriging. It can take an output to a call to \code{krige.control} or
    a list with elements as for the arguments in \code{krige.control}. }
}
%\details{
%  ~~ If necessary, more details than the description above ~~
%}
\value{
  An scalar with the value of the global variance
%  If it is a LIST, use
%  \item{comp1 }{Description of 'comp1'}
%  \item{comp2 }{Description of 'comp2'}
%  ...
}
\references{Isaaks, E.S and Srivastava, R.M. (1989) An Introduction to
  Applied Geostatistics, pag. 508, eq. 20.7. Oxford University Press. }

\author{
  Paulo Justiniano Ribeiro Jr. \email{paulojus@leg.ufpr.br}, \cr
  Peter J. Diggle \email{p.diggle@lancaster.ac.uk}.
}

%\note{ ~~further notes~~ 
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{\code{\link{krige.conv}} for the kriging algorithm. }
%\examples{
%}
\keyword{spatial}% at least one, from doc/KEYWORDS
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
