\name{kattegat}
\alias{kattegat}

\docType{data}
\title{Kattegat basin salinity data}
\description{
  Salinity measurements at the Kattegat basin, Denmark. 
}
\usage{data(kattegat)}
\format{
  An object of the \code{\link{class}}  \code{"geodata"},
  which is list with three components:  
  \describe{
    \item{\code{coords}}{the coordinates of the data locations. The
      distance are given in kilometers.  }
    \item{\code{data}}{values of the piezometric head. The unit is heads
      to meters.  }
    \item{\code{dk}}{a list with cooordinates of lines defining borders
      and islands across the study area.  }    
  }
}
%\details{
%  ~~ If necessary, more details than the __description__ above ~~
%}
\source{National Environmental Research Institute, 
  Arhus University, Denmark and the Swedish Meteorological and Hydrological Institute.
}

\references{
  Diggle, P. J. and Lophaven, S. (2006). Bayesian geostatistical
  design. \emph{Scandinavian Journal of Statistics}, 33: 55-64.
}

\examples{
plot(c(550,770),c(6150,6420),type="n",xlab="X UTM",ylab="Y UTM")
points(kattegat, add=TRUE)
lapply(kattegat$dk, lines, lwd=2)
}
\keyword{datasets}
