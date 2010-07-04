\name{head}
\alias{head}
\docType{data}
\title{Head observations in a regional confined aquifer}
\description{
  Measurements of potentiometric head at 29 locations in a regional
  confined sandstone aquifer.  Extract from Kitanidis' book.
}
\usage{data(head)}
\format{
  An object of the class \code{geodata} which is a list with the elements:
  \describe{
    \item{coords}{coordinates of the data location. }
    \item{data}{the data vector with head measurements (feet). }
  }
}
%\details{
%  ~~ If necessary, more details than the __description__ above ~~
%}
\source{
  Kitanidis, P.K. Introduction to geostatistics - applications in
  hidrology (1997). Cambridge University Press.
}
%\references{
%  ~~ possibly secondary sources and usages ~~
%}
\examples{
summary(head)
plot(head)
}
\keyword{datasets}
