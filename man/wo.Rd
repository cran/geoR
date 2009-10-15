\name{wo}
\alias{wo}
\docType{data}
\title{Kriging example data from Webster and Oliver}
\description{
  Data used in Chapter 8,
  page 156 of Webster and Oliver (2001) to illustrate
  properties of the kriging predictor.  
}
\usage{data(wo)}
\format{
  An object of the class \code{geodata} which is a list with the elements:
  \describe{
    \item{coords}{coordinates of the data location. }
    \item{data}{the data vector. }
    \item{x1}{coordinate of the centrally located prediction point. }
    \item{x2}{coordinate of the off-centre prediction point. }
    }
}
%\details{
%  ~~ If necessary, more details than the __description__ above ~~
%}
\source{
  Webster, R. and Oliver, M.A. (2001). Geostatistics for Environmental
  Scientists. Wiley.
}
%\references{
%  ~~ possibly secondary sources and usages ~~
%}
\examples{
attach(wo)
par(mfrow=c(1,2))
plot(c(-10,130), c(-10,130), ty="n", asp=1)
points(rbind(coords, x1))
text(coords[,1], 5+coords[,2], format(data))
text(x1[1]+5, x1[2]+5, "?", col=2)
plot(c(-10,130), c(-10,130), ty="n", asp=1)
points(rbind(coords, x2))
text(coords[,1], 5+coords[,2], format(data))
text(x2[1]+5, x2[2]+5, "?", col=2)
}
\keyword{datasets}


