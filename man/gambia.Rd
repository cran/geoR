\name{gambia}
\alias{gambia}
\alias{gambia.borders}
\alias{gambia.map}
\docType{data}
\title{Gambia Malaria Data}
\description{
  Malaria prevalence in children recorded at villages in The Gambia, Africa.
}
\usage{data(gambia)}
\format{
  Two objects are made available:
  \enumerate{
    \item \code{gambia}\cr
    A data frame with 2035 observations on the following 8 variables.
    \describe{
      \item{x}{x-coordinate of the village (UTM).  }
      \item{y}{y-coordinate of the village (UTM).  }
      \item{pos}{presence (1) or absence (0) of malaria in a blood
	sample taken from the child.}
      \item{age}{age of the child, in days}
      \item{netuse}{indicator variable denoting whether (1) or not (0)
	the child regularly sleeps under a bed-net.  }
      \item{treated}{indicator variable denoting
	whether (1) or not (0) the bed-net is
	treated (coded 0 if netuse=0).}
      \item{green}{satellite-derived measure of the green-ness of vegetation
	in the immediate vicinity of the village (arbitrary units).  }
      \item{phc}{indicator variable denoting the presence (1) or
	absence (0) of a health center in the
	village.  }
    }
    \item \code{gambia.borders}\cr
    A data frame with 2 variables:
    \describe{
      \item{x}{x-coordinate of the country borders.  }
      \item{y}{y-coordinate of the country borders.  }
    }
  }
}
%\details{
%  ~~ If necessary, more details than the __description__ above ~~
%}
%\source{
%  ~~ reference to a publication or URL from which the data were obtained ~~
%}
\references{
  Thomson, M., Connor, S., D Alessandro, U., Rowlingson, B., Diggle, P., Cresswell, M. & Greenwood, B.
(1999). Predicting malaria infection in Gambian children from satellite data and bednet use surveys: the
importance of spatial correlation in the interpretation of results. \emph{American Journal of Tropical
Medicine and Hygiene} 61: 2--8.

  Diggle, P., Moyeed, R., Rowlingson, B. & Thomson, M. (2002). Childhood
  malaria in The Gambia: a case-study in model-based geostatistics,
  \emph{Applied Statistics}.
}
\examples{
plot(gambia.borders, type="l", asp=1)
points(gambia[,1:2], pch=19)
# a built-in function for a zoomed map
gambia.map()
# Building a "village-level" data frame
ind <- paste("x",gambia[,1], "y", gambia[,2], sep="")
village <- gambia[!duplicated(ind),c(1:2,7:8)]
village$prev <- as.vector(tapply(gambia$pos, ind, mean))
plot(village$green, village$prev)
}
\keyword{datasets}

