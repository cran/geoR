\name{SIC}
\docType{data}

\alias{sic}
\alias{SIC}
\alias{sic.100}
\alias{sic.367}
\alias{sic.all}
\alias{sic.some}
\alias{sic.borders}

\title{Spatial Interpolation Comparison data}

\usage{data(SIC)}

\description{
  Data from the SIC-97 project: Spatial Interpolation Comparison.
}

\format{
  Four objects of the \code{\link{class}}  \code{"geodata"}:
  \code{sic.all}, \code{sic.100}, \code{sic.367}, \code{sic.some}.
  Each is a list with two components:
  
  \describe{
    \item{\code{coords}}{the coordinates of the data locations. The
      distance are given in kilometers.  }
    \item{\code{data}}{rainfall values. The unit is milimeters.  }
    \item{\code{altitude}}{elevation values. The unit is milimeters.  }
  }
  
  Additionally an matrix \code{sic.borders} with the borders of the country.
}

\source{
  Data from the project \emph{Spatial Interpolation Comparison 97};
see \url{https://soft.mines-paristech.fr/gstlearn/latest/data/benchmark/SIC97_Readme.pdf}.
  %  see \url{https://wiki.52north.org/bin/view/AI_GEOSTATS/EventsSIC97}.
}

\references{
  
  Christensen, O.F., Diggle, P.J. and Ribeiro Jr, P.J. (2001)  Analysing positive-valued spatial data: the transformed Gaussian model.  In Monestiez, P., Allard, D. and Froidevaux (eds), GeoENV III - Geostatistics for environmental applications. Quantitative Geology and Geostatistics, Kluwer Series, 11, 287--298.

}

\examples{
points(sic.100, borders=sic.borders)
points(sic.all, borders=sic.borders)
}

\keyword{datasets}
\keyword{spatial}
