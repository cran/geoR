\name{SIC}

\alias{sic}
\alias{sic.100}
\alias{sic.367}
\alias{sic.all}
\alias{sic.some}
\alias{sic.borders}

\non_function{}

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
  see \url{ftp://ftp.geog.uwo.ca/SIC97/}
}

%\references{
%  Cressie, N.A.C (1993) \emph{Statistics for Spatial Data}. New York: Wiley.
  
%  Papritz, A. and Moyeed, R. (2001) Parameter uncertainty in spatial
%  prediction: checking its importance by cross-validating the Wolfcamp
%  and Rongelap data sets. \emph{GeoENV 2000: Geostatistical for
%    Environmental Applications. Ed. P. Monestiez, D. Allard,
%    R. Froidevaux}. Kluwer. 
%}

\examples{
data(SIC)
points(sic.100, borders=sic.borders)
}

\keyword{datasets}
\keywords{spatial}
 