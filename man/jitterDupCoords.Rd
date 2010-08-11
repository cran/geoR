\name{jitterDupCoords}

\alias{jitterDupCoords}
\alias{jitterDupCoords.default}
\alias{jitterDupCoords.geodata}
\alias{jitter2d}

\title{
  Jitters (duplicated) coordinates.
}
\description{
  Jitters 2D coordinates uniformily on a region around (duplicated) points.
}
\usage{
jitter2d(coords, max, min = 0.2 * max, fix.one = TRUE,
         which.fix = c("random", "first", "last")) 

jitterDupCoords(x, ...)

\method{jitterDupCoords}{default}(x, ...)

\method{jitterDupCoords}{geodata}(x, ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, coords}{a matrix or data frame with 2D coordinates or geodata
    object.  }
  \item{max}{numeric scalar defining maximum jittering distance. }
  \item{min}{numeric scalar defining minimum jittering distance. }
  \item{fix.one}{logical. Whether or not one of the coordinates should
    not be jittered.  }
  \item{which.fix}{single element vector of integer or character,
    defining which coordinate won't be jittered.  Only used if \code{fix.one=TRUE}.  }
  \item{\dots}{arguments passed to \code{jitter2d}.  }
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
  \code{jitter2d} returns an object of the same type fo the input with
  jittered values

  \code{jitterDupCoords} returns an object of the same type fo the input with
  jittered coordinate values only at the duplicated locations
  
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}

%\references{
%% ~put references to the literature/web site here ~
%}

\author{
  Paulo Justiniano Ribeiro Jr. \email{paulojus@leg.ufpr.br}, \cr
  Peter J. Diggle \email{p.diggle@lancaster.ac.uk}.
}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{dup.coords}}, \code{\link{duplicated.geodata}} for
  functions identifying duplicated locations. 
}

\examples{
## simulating data
dt <- grf(30, cov.p=c(1, .3)) 
dt$coords <- round(dt$coords, dig=2)
## "forcing" some duplicated locations
dt$coords[4,] <- dt$coords[14,] <- dt$coords[24,] <- dt$coords[2,]
dt$coords[17,] <- dt$coords[23,] <- dt$coords[8,]

## jittering a matrix of duplicated coordinates
dt$coords[c(2,4,14,24),]
jitter2d(dt$coords[c(2,4,14,24),], max=0.01)

## jittering only the duplicated locations and comparing with original
cbind(dt$coords, jitterDupCoords(dt$coords, max=0.01))

## creating a now geodata object jittering the duplicated locations of the original one:
dup.coords(dt)
dt1 <- jitterDupCoords(dt, max=0.01)
dup.coords(dt1)
}
\keyword{spatial}
\keyword{manip}

