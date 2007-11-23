\name{s100 and s121}

\docType{data}

\alias{s100}
\alias{s121}


\title{Simulated Data-Sets which Illustrate the Usage of the Package
  geoR  }

\usage{
data(s100)

data(s121)
}

\description{
  These two simulated data sets are the ones used in the Technical
  Report which describes the package \pkg{geoR}
  (see reference below). 
  These data-sets are used in several examples throughout the package
  documentation.
}


\format{
  Two objects of the \code{\link{class}} \code{geodata}. Both are lists
  with the following components:
  \describe{
    \item{\code{coords}}{the coordinates of data locations.  }
    \item{\code{data}}{the simulated data. Notice that for \code{s121}
      this a \eqn{121 \times 10}{101x10} matrix with 10 simulations.  }
    \item{\code{cov.model}}{the correlation model.  }
    \item{\code{nugget}}{the values of the nugget parameter.  }
    \item{\code{cov.pars}}{the covariance parameters.  }
    \item{\code{kappa}}{the value of the parameter \emph{kappa}.  }
    \item{\code{lambda}}{the value of the parameter \emph{lambda}.  }
  }
}

\references{
  Ribeiro Jr, P.J. and Diggle, P.J. (1999) geoS: A geostatistical
  library for S-PLUS. \emph{Technical report ST-99-09, Dept of Maths and
  Stats, Lancaster University}. 
  
  Further information on the package \pkg{geoR} can be found at:\cr
  \url{http://www.leg.ufpr.br/geoR}.
}

\examples{
plot(s100)
plot(s121, type="l")
}
\keyword{datasets}

