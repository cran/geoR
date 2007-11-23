\name{.nlmP}

\alias{.nlmP}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Adapts nlm for Constraints in the Parameter Values}

\description{
  This function adapts the \R function \code{\link{nlm}} to allow for
constraints (upper and/or lower bounds) in the values of the parameters. 
}
\usage{
.nlmP(objfunc, params, lower=rep(-Inf, length(params)),
      upper=rep(+Inf, length(params)), \dots)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{objfunc}{the function to be minimized.  }
  \item{params}{starting values for the parameters.  }
  \item{lower}{lower bounds for the variables. Defaults to \eqn{-Inf}.  }
  \item{upper}{upper bounds for the variables. Defaults to \eqn{-Inf}.  }
  \item{\dots}{further arguments to be passed to the function
    \code{\link{nlm}}.  }
}

\details{
    Constraints on the parameter values are internally imposed by using
  exponential, logarithmic, and logit transformation of the parameter
values.
}

\value{
  The output is the same as for the function \code{\link{nlm}}.
}

\references{
  Further information on the package \pkg{geoR} can be found at:\cr
  \url{http://www.leg.ufpr.br/geoR}.
}


\author{
  Patrick E. Brown \email{p.brown@lancaster.ac.uk}.\cr
  Adapted and included in \pkg{geoR} by \cr
  Paulo Justiniano Ribeiro Jr. \email{paulojus@leg.ufpr.br}

}


%\note{ ~~further notes~~ }

% ~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{\code{\link{nlm}}, \code{\link{optim}}.  }

%\examples{
%}

\keyword{spatial}
\keyword{optimize}%-- one or more ...
