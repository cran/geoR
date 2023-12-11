\name{hoef}
\docType{data}

\alias{hoef}

\title{Data for spatial analysis of experiments}
\usage{data(hoef)}
\description{
  The \code{hoef} data frame has 25 rows and 5 columns.\cr
  The data consists of a uniformity trial for which \emph{artificial}
  treatment effects were assign to the plots.
}

\format{
  This data frame contains the following columns:
  \describe{
    \item{x1}{x-coordinate of the plot.}
    \item{x2}{y-coordinate of the plot.}
    \item{dat}{the \emph{artificial} data.}
    \item{trat}{the treatment number.}
    \item{ut}{the data from the uniformity trial, without the treatment effect.}
  }
}

\details{
  The treatment effects assign to the plots are:
  \itemize{
    \item Treatment 1: \eqn{\tau_1 = 0}{tau1 = 0}
    \item Treatment 2: \eqn{\tau_2 = -3}{tau2 = -3}
    \item Treatment 3: \eqn{\tau_3 = -5}{tau3 = -5}
    \item Treatment 4: \eqn{\tau_4 = 6}{tau4 = 6}
    \item Treatment 5: \eqn{\tau_5 = 6}{tau5 = 6}
  }    
}

\source{
  Ver Hoef, J.M. & Cressie, N. (1993) Spatial statistics: analysis of field experiments.
  In Scheiner S.M. and Gurevitch, J. (Eds) \emph{Design and Analysis of
  Ecological Experiments}. Chapman and Hall. 
}
%\references{
%~~ possibly secondary sources and usages ~~
%}
\examples{
hoef.geo <- as.geodata(hoef, covar.col=4)
summary(hoef)
summary(hoef.geo)
points(hoef.geo, cex.min=2, cex.max=2, pt.div="quintiles")
}
\keyword{datasets}
