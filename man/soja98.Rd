\name{soja98}
\alias{soja98}
\docType{data}
\encoding{latin1}

\title{Soya bean production and other variables in a uniformity trial}
\description{
  Data on soya bean production in a uniformity trial measured at
  plots of size 5x5 meters and other soil properties measured
  in points given by the data coordinates.
}
\usage{data(soja98)}
\format{
  A data frame with 256 observations on the following 10 variables.
  \describe{
    \item{\code{X}}{a numeric vector with X-coordinates of the plot centres.  }
    \item{\code{Y}}{a numeric vector with X-coordinates of the plot centres.  }
    \item{\code{P}}{a numeric vector, phosphorous content.  }
    \item{\code{PH}}{a numeric vector, soil pH.  }
    \item{\code{K}}{a numeric vector, potassium content. (Cmol/dm^3)  }
    \item{\code{MO}}{a numeric vector, organic matter. (percentage) }
    \item{\code{SB}}{a numeric vector, basis saturation.  }
    \item{\code{iCone}}{a numeric vector, cone index, measuring mechanic
    resistence of the soil. (kg/cm^2)  }
    \item{\code{Rend}}{a numeric vector, total soya production within
      the plot (kg). }
    \item{\code{PROD}}{a numeric vector, production converted to
      productivity by a unit of area - hectare (ton/ha). }
  }
}

\source{
Souza, E.G.; Jojann, J. A.; Rocha, J. V.; Ribeiro, S. R. A.; Silva, M. S., Upazo, M. A. U.; Molin, J. P.; Oliveira, E. F.; Nóbrega, L. H. P. (1999) Variabilidade espacial dos atributos químicos do solo em um latossolo roxo distrófico da região de Cascavel-PR. \emph{Engenharia Agrícola}. Jaboticabal, volume 18, nr. 3, p.80-92.
}

%\references{
%  ~~ possibly secondary sources and usages ~~
%}

\examples{
data(soja98)
plot(soja98)
require(geoR)
prod98 <- as.geodata(soja98, coords.col=1:2, data.col=)
plot(prod98, low=TRUE)
}

\keyword{datasets}


