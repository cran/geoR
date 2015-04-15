\name{subarea}
\alias{subarea}

\title{Selects a Subarea from a Geodata Object}
\description{
  Selects elements of a \code{geodata} object wich are within a
  rectangular (sub)area
}
\usage{
subarea(geodata, xlim, ylim, \ldots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{geodata}{an object of the class \code{geodata} as defined in \code{\link{as.geodata}}.  }
  \item{xlim}{optional, a vector with selected range of the x-coordinates. }
  \item{ylim}{optional, a vector with selected range of the y-coordinates. }
  \item{\ldots}{further arguments to be passed to \code{\link{zoom.coords}}.  }
}
\details{
  The function copies the original \code{geodata} object and
  selects values of \code{$coords}, \code{$data}, \code{$borders},
  \code{$covariate} and \code{$units.m} which lies within the selected
  sub-area.
  Remaining components of the geodata objects are untouched.
  
  If \code{xlim} and/or \code{ylim} are not provided the function
  prompts the user to click 2 points defining an retangle defining the 
  subarea on a existing plot.

}
\value{
  Returns an \code{geodata} object, subsetting of the original one provided.
}
%\references{ ~put references to the literature/web site here ~ }
\author{
  Paulo Justiniano Ribeiro Jr. \email{paulojus@leg.ufpr.br}, \cr
  Peter J. Diggle \email{p.diggle@lancaster.ac.uk}.
}
%\note{ ~~further notes~~ }

% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{zoom.coords}}, \code{\link{locator}}}

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
legend("topleft", c("foo", "foo1 (zoom.coords)", "foo2 (coords2coords)"), 
       lty=c(1,2,1), lwd=c(1,1,2))

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

if(require(geoRglm)){
data(rongelap)
points(rongelap, bor=borders)
## zooming the western area
rongwest <- subarea(rongelap, xlim=c(-6300, -4800))
points(rongwest, bor=borders)
## now zooming in the same plot
points(rongelap, bor=borders)
rongwest.z <- zoom.coords(rongwest, xzoom=3, xoff=2000, yoff=3000)
points(rongwest.z, bor=borders, add=TRUE)
rect.coords(rongwest$sub, quiet=TRUE)
rect.coords(rongwest.z$sub, quiet=TRUE)
}
}
\keyword{spatial}% at least one, from doc/KEYWORDS
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
