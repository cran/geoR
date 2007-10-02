\name{tce}
\alias{tce}
\docType{data}
\title{TCE concentrations in groundwater in a vertical cross section}
\description{
  Measurements at 56 locations of concentration of trichloroethylene
  (TCE) in groundwater on a transect in a fine-sand superficial aquifer.
  Extract from Kitanidis' book.
}

\usage{data(tce)}
\format{
  An object of the class \code{geodata} which is a list with the elements:
  \describe{
    \item{coords}{coordinates of the data location (feet). }
    \item{data}{the data vector with measurements of the TCE
      concentration (ppb). }
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
summary(tce)
summary(tce, lambda=0)
plot(tce)
points(tce)
points(tce, lambda=0)
}
\keyword{datasets}
