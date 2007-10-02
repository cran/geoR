\name{wolfcamp}
\docType{data}

\alias{wolfcamp}
\alias{wolf}

\title{Wolfcamp Aquifer Data}
\usage{data(wolfcamp)}

\description{
  Piezometric head measurements taken at the Wolfcamp Aquifer, Texas,
  USA.
  See Cressie (1993, p.212--214) for description of the scientific
  problem and the data.
  Original data were converted to SI units: coordinates are given in
  kilometers and pressure heads to meters.
}

\format{
  An object of the \code{\link{class}}  \code{"geodata"}, which is list with two components:
  
  \describe{
    \item{\code{coords}}{the coordinates of the data locations. The
      distance are given in kilometers.  }
    \item{\code{data}}{values of the piezometric head. The unit is heads
      to meters.  }
  }
}

\source{
  Harper, W.V and Furr, J.M. (1986) Geostatistical analysis of
  potentiometric data in the Wolfcamp Aquifer of the Palo Duro Basin,
  Texas. \emph{Technical Report BMI/ONWI-587, Bettelle Memorial Institute,
    Columbus, OH}.   
}

\references{
  Cressie, N.A.C (1993) \emph{Statistics for Spatial Data}. New York: Wiley.
  
  Papritz, A. and Moyeed, R. (2001) Parameter uncertainty in spatial
  prediction: checking its importance by cross-validating the Wolfcamp
  and Rongelap data sets. \emph{GeoENV 2000: Geostatistical for
    Environmental Applications. Ed. P. Monestiez, D. Allard,
    R. Froidevaux}. Kluwer. 
}

\examples{
summary(wolfcamp)
plot(wolfcamp)
}

\keyword{datasets}


