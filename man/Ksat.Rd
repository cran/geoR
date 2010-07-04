\name{Ksat}
\alias{Ksat}
\docType{data}
\encoding{latin1}
\title{Saturated Hydraulic Conductivity}
\description{
  The data consists of 32 measurements of the saturated hydraulic
  conductivity of a soil. 
}
\usage{data(Ksat)}
\format{
  The object \code{Ksat} is a list of the class \code{geodata}
  with the following elements:
  \describe{    
    \item{\code{coords}}{a matrix with the coordinates of the soil samples.  }
    \item{\code{data}}{measurements of the saturated hidraulic conductivity.  }
    \item{\code{borders}}{a data-frame with the coordinates of a polygon defining the
      borders of the area.  }
    }
}
%\details{
%  ~~ If necessary, more details than the __description__ above ~~
%}
\source{
  Data provided by Dr. Décio Cruciani, ESALQ/USP, Brasil.
}
%\references{
%  ~~ possibly secondary sources and usages ~~
%}
\examples{
summary(Ksat)
plot(Ksat, border=borders)
}
\keyword{datasets}
