\name{isaaks}
\docType{data}

\alias{isaaks}

\title{Data from Isaaks and Srisvastava's book}
\description{
  Toy example used in the book \emph{An Introduction to Geostatistics} to
  illustrate the effects of different models
  and parameters in the kriging results when predicting at a given point.
}
\usage{data(isaaks)}
\format{
  An object of the class \code{geodata} which is a list with the elements:
  \describe{
    \item{coords}{coordinates of the data location. }
    \item{data}{the data vector. }
    \item{x0}{coordinate of the prediction point. }
  }
}
%\details{
%  ~~ If necessary, more details than the __description__ above ~~
%}
\source{
  Isaaks, E.H. & Srisvastava, R.M. (1989)
  An Introduction to Applied Geostatistics.
  Oxford University Press.
}

%\references{
%  ~~ possibly secondary sources and usages ~~
%}
\examples{
isaaks
summary(isaaks)
plot(isaaks$coords, asp=1, type="n")
text(isaaks$coords, as.character(isaaks$data))
points(isaaks$x0, pch="?", cex=2, col=2)
}
\keyword{datasets}
