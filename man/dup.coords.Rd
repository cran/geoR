\name{dup.coords}
\alias{dup.coords}
\alias{dup.coords.default}
\alias{dup.coords.geodata}

\title{Locates duplicated coordinates}
\description{
  This funtions takes an object with 2-D coordinates and returns the
  positions of the duplicated coordinates.
}
\usage{
dup.coords(x)
dup.coords.default(x)
dup.coords.geodata(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a two column numeric matrix or data frame.}
}
%\details{
%  ~~ If necessary, more details than the __description__  above ~~
%}
\value{
  A list where each component is a vector with the positions  or the
  rownames, if different from \code{NULL}, of the duplicates coordinates. 
}
%\references{ ~put references to the literature/web site here ~ }
\author{
  Paulo Justiniano Ribeiro Jr. \email{Paulo.Ribeiro@est.ufpr.br} \cr
  Peter J. Diggle \email{p.diggle@lancaster.ac.uk}.
}

%\note{ ~~further notes~~ }

% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{as.geodata}}}
\examples{
## simulating data
foo <- grf(30, cov.p=c(1, .3)) 
## forcing some duplicated locations
foo$coords[4,] <- foo$coords[14,] <- foo$coords[24,] <- foo$coords[2,]
foo$coords[17,] <- foo$coords[23,] <- foo$coords[8,]
## locations the duplicated coordinates
dup.coords(foo)
}
\keyword{spatial}
\keyword{manip}
%\keyword{logic}
