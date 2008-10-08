\name{soil250}
\alias{soil250}
\docType{data}
\title{Soil chemistry properties data set}
\description{
  Several soil chemistry properties measured on a regular grid with 10x25 points
  spaced by 5 meters.
}
\usage{data(soil250)}
\format{
  A data frame with 250 observations on the following 22 variables.
  \describe{
    \item{Linha}{x-coordinate}
    \item{Coluna}{y-coordinate}
    \item{Cota}{elevation}
    \item{AGrossa}{a numeric vecto, sand portion of the sample.}
    \item{Silte}{a numeric vector, silt portion of the sample.}
    \item{Argila}{a numeric vector, sand portion of the sample.}
    \item{pHAgua}{a numeric vector, soil pH at water }
    \item{pHKCl}{a numeric vector, soil pH by KCl }
    \item{Ca}{a numeric vector, calcium content }
    \item{Mg}{a numeric vector, magnesium content }
    \item{K}{a numeric vector, potassio content }
    \item{Al}{a numeric vector, aluminium content }
    \item{H}{a numeric vector, hidrogen content }
    \item{C}{a numeric vector, carbon content }
    \item{N}{a numeric vector, nitrogen content }
    \item{CTC}{a numeric vector, catium exchange capability }
    \item{S}{a numeric vector, enxofrar content }
    \item{V}{a numeric vector}
    \item{M}{a numeric vector}
    \item{NC}{a numeric vector}
    \item{CEC}{a numeric vector}
    \item{CN}{a numeric vector, carbon/nitrogen relation }
  }
}
\details{
  Uniformity trial with 250 undisturbed soil samples collected at 25cm
  soil depth of spacing of 5 meters, resulting on a regular grid of
  \eqn{25 \times 10}{25 x 10} points.

  See also the data-set \emph{wrc} with other variables colected at the same points.
}
\source{
  Bassoi thesis
}
\references{
  Bassoi papers
}
\examples{
data(soil250)
ctc <- as.geodata(soil250, data.col=16)
plot(ctc)
}
\keyword{datasets}
\keyword{spatial}
