\name{s256i}

\alias{s256i}

\non_function{s256i}

\title{Simulated Data-Set which Illustrate the Usage of krige.bayes  }

\description{
  This is the simulated data-set used in the Technical
  Report which describes the implementation of the function \bold{krige.bayes}
  (see reference below). 
}

\usage{
data(s256i)
}

\format{
  Two objects of the \code{\link{class}} \code{geodata}. Both are lists
  with the following components:
  \describe{
    \item{\code{coords}}{the coordinates of data locations.  }
    \item{\code{data}}{the simulated data.  }
  }
}

\references{
  Ribeiro Jr, P.J. and Diggle, P.J. (1999)
  Bayesian inference in Gaussian model-based geostatistics.
  \emph{Technical report ST-99-08, Dept of Maths and
    Stats, Lancaster University}. 
  
  Further information about \bold{geoR} can be found at:\cr
  \url{http://www.maths.lancs.ac.uk/~ribeiro/geoR}.
}

\keyword{datasets}
\keyword{spatial}