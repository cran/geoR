\name{ca20}
\alias{ca20}

\docType{data}

\title{Calcium content in soil samples}
\usage{data(ca20)}
\description{
  This data set contains the calcium content measured in soil samples
  taken from the 0-20cm layer at 178 locations within a certain study area divided in three
  sub-areas. The elevation at each location was also recorded.

  The first region is tipically flooded during the rain season and not
  used as an experimental area. The calcium levels would represent
  the natural content in the region. The second region has received
  fertilizers a while ago and is tipically occupied by rice fields.
  The third region has recieved fertilizers recently and is frequently
  used as an experimental area.
}

\format{
  The object \code{ca20} belongs to the class \code{geodata}
  and is a list
  with the following elements:
  
  \describe{    
    \item{\code{coords}}{a matrix with the coordinates of the soil samples.  }
    \item{\code{data}}{calcium content measured in \eqn{mmol_c/dm^3}.  }
    \item{\code{covariate}}{a data-frame with the covariates
      \describe{
	\item{\code{altitude}}{a vector with the elevation of each
	  sampling location, in meters (\eqn{m}).}
	\item{\code{area}}{a factor indicating the sub area to which the
	  locations belongs.}
      }
    }
    \item{\code{borders}}{a matrix with the coordinates defining the
      borders of the area.  }
    \item{\code{reg1}}{a matrix with the coordinates of the limits of
      the sub-area 1. }
    \item{\code{reg1}}{a matrix with the coordinates of the limits of
      the sub-area 2. }
    \item{\code{reg1}}{a matrix with the coordinates of the limits of
      the sub-area 3. }
  }
}

\source{
  The data was collected by researchers from PESAGRO and EMBRAPA-Solos, Rio
  de Janeiro, Brasil
  and provided by Dra. Maria Cristina Neves de Oliveira.

  Capeche, C.L. \emph{et al.} (1997) Caracterização pedológica da
  fazenda Angra - PESAGRO/RIO - Estação experimental de Campos (RJ).
  CD-ROM In: \emph{Congresso BRASILEIRO de Ciência do Solo} 26,
  Informação, globalização, uso do solo; Rio de Janeiro,
  1997. trabalhos. Rio de Janeiro: Embrapa/SBCS, 1997.
}

\references{
  Further information on the package \pkg{geoR} can be found at:\cr
  \url{http://www.est.ufpr.br/geoR}.
}

\keyword{datasets}
