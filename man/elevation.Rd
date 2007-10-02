\name{elevation}
\docType{data}
\alias{elevation}
\title{Surface Elevations}
\description{
  Surface elevation data taken from Davis (1972).
  An onject of the class \code{geodata} with elevation values at 52 locations.
}
\usage{
data(elevation)
}
\format{
  An object of the class \code{geodata} which is a list with the
  following elements:
  \describe{
    \item{\code{coords}}{x-y coordinates (multiples of 50 feet).}
    \item{\code{data}}{elevations (multiples of 10 feet).}
  }
}
\source{
  Davis, J.C. (1973) \emph{Statistics and Data Analysis in Geology.} Wiley.
}
\examples{
summary(elevation)
plot(elevation)
}
\keyword{datasets}
