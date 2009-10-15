\name{wrc}
\alias{wrc}
\docType{data}
\encoding{latin1}
\title{Points of a water retention curve data set}
\description{
  Soil density and measures of the water retention curve obtained at
  different pressures on a regular grid with 10x25 points
  spaced by 5 meters.
}
\usage{data(wrc)}
\format{
  A data frame with 250 observations on the following 11 variables.
  \describe{
    \item{\code{CoordX}}{a numeric vector with the X coordinates of the
      samples. }
    \item{\code{CoordY}}{a numeric vector with the Y coordinate of the
      samples. }
    \item{\code{Densidade}}{a numeric vector, soil density \eqn{(g/cm^3)}}
    \item{\code{Pr5}}{a numeric vector, water content at a pressure of 5
      mca -- \eqn{5 \times 10^2}{5 x 10^2} Pa (atm)}
    \item{\code{Pr10}}{a numeric vector, water content at a pressure of
      10 mca -- \eqn{1 \times 10^3}{1x10^3} Pa (atm)}
    \item{\code{Pr60}}{a numeric vector, water content at a pressure of 60 mca -- \eqn{6 \times 10^3}{6x10^3} Pa (atm)}
    \item{\code{Pr100}}{a numeric vector, water content at a pressure of 100 mca -- \eqn{1 \times 10^4}{1x10^4} Pa (atm)}
    \item{\code{Pr306}}{a numeric vector, water content at a pressure of 306 mca -- \eqn{3 \times 10^4}{3x10^4} Pa (atm)}
    \item{\code{Pr816}}{a numeric vector, water content at a pressure of 816 mca -- \eqn{8 \times 10^4}{8x10^4} Pa (atm)}
    \item{\code{Pr3060}}{a numeric vector, water content at a pressure of 3060 mca -- \eqn{3 \times 10^5}{8x10^5} Pa (atm)}
    \item{\code{Pr15300}}{a numeric vector, water content at a pressure of 15300 mca -- \eqn{1.5 \times 10^6}{1.5x10^6} Pa (atm)}
  }
}

\details{
  Uniformity trial with 250 undisturbed soil samples collected at 25cm
  soil depth of spacing of 5 meters, resulting on a regular grid of
  \eqn{25 \times 10}{25 x 10} sampling points.

  For each sampling point there are measurents of the soil density and
  water content obtained at eight pressures: 5, 10, 60, 100, 306, 816,
  3060 and 15300 meters of column of water (mca), corresponding to
  \eqn{5 \times 10^2}{5 x 10^2}, 
  \eqn{1 \times 10^3}{1 x 10^3}, \eqn{6 \times 10^3}{6 x 10^3}, \eqn{1 \times 10^4}{1 x 10^4},
  \eqn{3 \times 10^4}{3 x 10^}, \eqn{8 \times 10^4}{8 x 10^4}, \eqn{3 \times 10^5}{3 x 10^5},
  \eqn{1.5 \times10^6}{1.5 x 10^6} Pa.

  The experiment aimed to use the water contents of the samples to
  estimate the water retention curve at the 250 data points.

  See also the data-set \code{\link{soil250}} with soil chemistry properties measured at the same points.
}

\source{
  MORAES, S.O. (1991) Heterogeneidade hidráulica de uma terra roxa
  estruturada. PhD Thesis. ESALQ/USP.
}

\references{
  MORAES, S. O. ; LIBARDI, P. L. ; REICHARDT, K. (1993) Problemas
  metodológicos na obtenção da curva de retenção de água pelo
  solo. Scientia Agricola, Piracicaba, v. 50, n. 3, p. 383-392.
  
  MORAES, S. O. ; LIBARDI, P. L. ; REICHARDT, K. ; BACCHI,
  O. O. S. (1993) Heterogeneidade dos pontos experimentais
  de curvas de retenção da água do solo.. Scientia Agricola, Piracicaba, v. 50, n. 3, p. 393-402.

  MORAES, S. O. ; LIBARDI, P. L. (1993) Variabilidade da água disponível
  em uma terra roxa estruturada latossólica.
  Scientia Agricola, Piracicaba, v. 50, n. 3, p. 393-402, 1993. 
}

\examples{
pr100 <- as.geodata(wrc, data.col=7)
summary(pr100)
plot(pr100)
}
\keyword{datasets}
\keyword{spatial}

