\name{landim1}
\encoding{latin1}
\alias{landim1}
\docType{data}
\title{Data from Landim's book}
\description{
  Artificial or non-specified data from Paulo Landim's book 
}
\usage{data(landim1)}
\format{
  A data frame with 38 observations on the following 4 variables.
  \describe{
    \item{\code{EW}}{a numeric vector with the east-west coordinates. }
    \item{\code{NS}}{a numeric vector with the north-south coordinates. }
    \item{\code{A}}{a numeric vector with data on a first variable. }
    \item{\code{B}}{a numeric vector with data on a second variable. }
  }
}
%\details{
%  ~~ If necessary, more details than the __description__ above ~~
%}
\source{
 Landim, P. M. B. (2004) \emph{Análise estatística de dados geológicos}. Editora Unesp. 
 Data from Table~1, pg.12.
}
%\references{
%  ~~ possibly secondary sources and usages ~~
%}
\examples{
data(landim)
plot(as.geodata(landim1, data.col=3))
plot(as.geodata(landim1, data.col=4))
}
\keyword{datasets}
