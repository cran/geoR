\name{soil}
\alias{soil}
\docType{data}
\title{Soil chemistry properties data set}
\description{
  Several soil chemistry properties measured on a 10x25 grid of points.
}
\usage{data(soil)}
\format{
  A data frame with 250 observations on the following 22 variables.
  \describe{
    \item{Linha}{x-coordinate}
    \item{Coluna}{y-coordinate}
    \item{Cota}{elevation}
    \item{AGrossa}{a numeric vector}
    \item{Silte}{a numeric vector}
    \item{Argila}{a numeric vector}
    \item{pHAgua}{a numeric vector}
    \item{pHKCl}{a numeric vector}
    \item{Ca}{a numeric vector, calcium content}
    \item{Mg}{a numeric vector, magnesium content}
    \item{K}{a numeric vector, potassio content}
    \item{Al}{a numeric vector, aluminium content}
    \item{H}{a numeric vector, hidrogen content}
    \item{C}{a numeric vector, carbon content}
    \item{N}{a numeric vector, nitrogenm content}
    \item{CTC}{a numeric vector, catium exchange capability}
    \item{S}{a numeric vector}
    \item{V}{a numeric vector}
    \item{M}{a numeric vector}
    \item{NC}{a numeric vector}
    \item{CEC}{a numeric vector}
    \item{CN}{a numeric vector}
  }
}
%\details{
%  ~~ If necessary, more details than the __description__ above ~~
%}
\source{
  Bassoi and Moraes
}
\references{
  Bassoi thesis...
}
\examples{
data(soil)
ctc <- as.geodata(soil, data.col=16)
plot(ctc)
}
\keyword{datasets}
\keyword{spatial}
