\name{read.geodata}

\alias{read.geodata}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Reads and Converts Data to geoR Format}

\description{
  Reads data from a \emph{ASCII} file and converts it to an object of the
  \code{\link{class}} \code{geodata}, the standard data format for the
  \bold{geoR} package.
}

\usage{
read.geodata(file, header = FALSE, coords.col = 1:2, data.col = 3,
             data.names = NULL, covar.col = NULL,
             covar.names = "header", ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{file}{a string with the name of the \emph{ASCII} file.  }
  \item{header}{logical.
    Indicates whether the variables names should be
    read from the first line of the input file.   }
  \item{coords.col}{a vector with the numbers of the columns containing the
    coordinates.  }
  \item{data.col}{a scalar or vector with the number of the column(s)
    containing the data.  }
  \item{data.names}{a string or vector of strings with names for
    the data columns. Only valid if there is more than one column of
    data.
    By default the names in the original object are used.  }
  \item{covar.col}{optional. A scalar or vector with the number
    of the column(s) with the values of the covariate(s).  }
  \item{covar.names}{a string or vector of strings with the name(s) of
    the covariates. By default the names in the original object are used.  }
  \item{\dots}{further arguments to be passed to the function \code{\link{read.table}}.  }
}

\details{
  The function \code{\link{read.table}} is used to read the data from the
\emph{ASCII} file and then \code{\link{as.geodata}} is used to convert
to an object of the \code{\link{class}} \code{geodata}.
}

\value{
  An object of the \code{\link{class}} \code{geodata}.
  See documentation for the function \code{\link{as.geodata}} for
  further details.
}

\references{
  Further information about \bold{geoR} can be found at:\cr
  \url{http://www.maths.lancs.ac.uk/~ribeiro/geoR.html}.
}

\author{
  Paulo Justiniano Ribeiro Jr. \email{Paulo.Ribeiro@est.ufpr.br}, \cr
  Peter J. Diggle \email{p.diggle@lancaster.ac.uk}.
}

%\note{ ~~further notes~~ }
%
% ~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{\code{\link{as.geodata}} to convert existing \R objects,
  \code{\link{read.table}}, the basic \R function used to  read \emph{ASCII} files,
  and \code{\link{list}} for detailed information about lists.}

%\examples{}
\keyword{spatial}
\keyword{manip}%-- one or more ...
