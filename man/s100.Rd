\name{s100 and s121}

\alias{s100}
\alias{s121}

\non_function{s100, s121}

\title{Simulated Data-Sets which Illustrate the Usage of the Package
  geoR  }

\description{
  These two simulated data sets are the ones used in the Technical
  Report which describes the package \bold{geoR}
  (see reference below). 
  These data-sets are used in several examples throughout the package
  documentation.
}

\usage{
data(s100)
}

\format{
  Two objects of the \code{\link{class}} \code{geodata}. Both are lists
  with the following components:
  \describe{
    \item{\code{coords}}{the coordinates of data locations.  }
    \item{\code{data}}{the simulated data.  }
    \item{\code{cov.model}}{the correlation model.  }
    \item{\code{nugget}}{the values of the nugget parameter.  }
    \item{\code{cov.pars}}{the covariance parameters.  }
    \item{\code{kappa}}{the value of the parameter \emph{kappa}.  }
    \item{\code{lambda}}{the value of the parameter \emph{lambda}.  }
  }
}

\references{
  Ribeiro Jr, P.J. and Diggle, P.J. (1999). geoS: A geostatistical
  library for S-PLUS. Technical report ST-99-09, Dept of Maths and
  Stats, Lancaster University. 
  
  Further information about \bold{geoR} can be found at:\cr
  \url{http://www.maths.lancs.ac.uk/~ribeiro/geoR.html}.
}

\keyword{datasets}
\keyword{spatial}
