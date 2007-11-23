\name{eyefit}

\alias{eyefit}
\alias{lines.eyefit}
\alias{summary.eyefit}
\alias{print.summary.eyefit}
\alias{print.eyefit}
\alias{plot.eyefit}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Interactive Variogram Estimation}
\description{
  Function to fit an empirical variogram "by eye" using an interactive
  Tcl-Tk interface.
}
\usage{
eyefit(vario, silent = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vario}{An empirical variogram object as returned by the function
  \code{\link{variog}}.  }
  \item{silent}{logical indicating wheather or not the fitted variogram
    must be returned.  }
}
%\details{
%  ~~ If necessary, more details than the __description__  above ~~
%}
\value{
  Returns a list of list
  with the model parameters for each of the saved fit(s).
}

%\references{ ~put references to the literature/web site here ~ }
\author{Andreas Kiefer \email{andreas@inf.ufpr.br} \cr
        Paulo Justiniano Ribeiro Junior \email{paulojus@leg.ufpr.br}. }
%\note{ ~~further notes~~ }

% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{variofit}} for least squares variogram fit,
  \code{\link{likfit}} for likelihood based parameter estimation
  and \code{\link{krige.bayes}} to obtain the posterior distribution for the model parameters. }
%\examples{
%}
\keyword{spatial}% at least one, from doc/KEYWORDS
\keyword{models}% __ONLY ONE__ keyword per line
\keyword{dynamic}% __ONLY ONE__ keyword per line
