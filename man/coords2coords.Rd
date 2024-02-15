\name{coords2coords}
\alias{coords2coords}
\alias{zoom.coords}
\alias{zoom.coords.default}
\alias{zoom.coords.geodata}
\alias{rect.coords}

\title{Operations on Coordinates}
\description{
  Functions for shifting, zooming and envolving rectangle of a set of coordinates.
}
\usage{
coords2coords(coords, xlim, ylim, xlim.ori, ylim.ori)

zoom.coords(x, \ldots)

\method{zoom.coords}{default}(x, xzoom, yzoom, xlim.ori, ylim.ori, xoff=0, yoff=0, \ldots)

\method{zoom.coords}{geodata}(x, \ldots)

rect.coords(coords, xzoom = 1, yzoom=xzoom, add.to.plot=TRUE,
            quiet = FALSE, \ldots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{coords, x}{two column matrix or data-frame with coordinates.  }
  \item{xlim}{range of the new x-coordinates. }
  \item{ylim}{range of the new y-coordinates. }
  \item{xlim.ori}{optional. Range of the original x-coordinates, by
    default the range of the original x-coordinates.  }
  \item{ylim.ori}{optional. Range of the original y-coordinates, by
    default the range of the original y-coordinates.  }
  \item{xzoom}{scalar, expanding factor in the x-direction.  }
  \item{yzoom}{scalar, expanding factor in the y-direction.  }
  \item{xoff}{scalar, shift in the x-direction.  }
  \item{yoff}{scalar, shift in the y-direction.  }
  \item{add.to.plot}{logical, if \code{TRUE} the retangle is added to
    the current plot.  }
  \item{quiet}{logical, none is returned.  }
  \item{\ldots}{further arguments to be passed to \code{\link{rect}}.  }
}
%\details{
%  ~~ If necessary, more details than the __description__  above ~~
%}
\value{
  \item{coords2coords and zoom.coords}{return an object of the same type as given in the argument
  \code{coords} with the transformed coordinates.  }
  \item{rect.coords}{returns a matrix with the 4 coordinates of
  the rectangle defined by the coordinates.  }
}
%\references{ ~put references to the literature/web site here ~ }
%\author{ ~~who you are~~ }
\author{
  Paulo Justiniano Ribeiro Jr. \email{paulojus@leg.ufpr.br}, \cr
  Peter J. Diggle \email{p.diggle@lancaster.ac.uk}.
}
%\note{ ~~further notes~~ }

% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{subarea}}, \code{\link{rect}}}

\examples{
foo <- matrix(c(4,6,6,4,2,2,4,4), nc=2)
foo1 <- zoom.coords(foo, 2)
foo1
foo2 <- coords2coords(foo, c(6,10), c(6,10))
foo2
plot(1:10, 1:10, type="n")
polygon(foo)
polygon(foo1, lty=2)
polygon(foo2, lwd=2)
arrows(foo[,1], foo[,2],foo1[,1],foo1[,2], lty=2)
arrows(foo[,1], foo[,2],foo2[,1],foo2[,2])
legend("topleft", 
       c("foo", "foo1 (zoom.coords)", "foo2 (coords2coords)"), lty=c(1,2,1), lwd=c(1,1,2))

## "zooming" part of The Gambia map
gb <- gambia.borders/1000
gd <- gambia[,1:2]/1000
plot(gb, ty="l", asp=1, xlab="W-E (kilometres)", ylab="N-S (kilometres)")
points(gd, pch=19, cex=0.5)
r1b <- gb[gb[,1] < 420,]
rc1 <- rect.coords(r1b, lty=2)

r1bn <- zoom.coords(r1b, 1.8, xoff=90, yoff=-90)
rc2 <- rect.coords(r1bn, xz=1.05)
segments(rc1[c(1,3),1],rc1[c(1,3),2],rc2[c(1,3),1],rc2[c(1,3),2], lty=3)

lines(r1bn)
r1d <- gd[gd[,1] < 420,]
r1dn <- zoom.coords(r1d, 1.7, xlim.o=range(r1b[,1],na.rm=TRUE), ylim.o=range(r1b[,2], na.rm=TRUE), 
                    xoff=90, yoff=-90)
points(r1dn, pch=19, cex=0.5)
text(450,1340, "Western Region", cex=1.5)
}
\keyword{spatial}% at least one, from doc/KEYWORDS
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
