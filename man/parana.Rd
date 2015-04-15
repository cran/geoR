\name{parana}
\alias{parana}
\alias{maijun}
\docType{data}
\encoding{latin1}
\title{Rainfall Data from Parana State, Brasil}
\usage{data(parana)}

\description{
  This data-set was used by Diggle and Ribeiro (2001) to illustrate the methods
  discussed in the paper. The data reported analysis was
  carried out using the package \pkg{geoR}. 
  
  The data refers to average rainfall over different years for the period May-June
  (dry-season). It was collected at 143 recording stations throughout \eqn{\mbox{Paran\'{a}}}{Parana} State,
  Brasil.
}

\format{
  The object \code{parana} of the class \code{geodata}, which is a list
  containing the following components:
  
  \describe{    
    \item{\code{coords}}{a matrix with the coordinates of the recording
      stations.  }
    \item{\code{data}}{a vector with the average recorded rainfall for the May-June period.  }
    \item{\code{borders}}{a matrix with the coordinates defining the borders of
      \eqn{\mbox{Paran\'{a}}}{Parana} state.  }
    \item{\code{loci.paper}}{a matrix with the coordinates of the four
      prediction locations discussed in the paper. }
  }
}

\source{
  The data were collected at several recording stations at
  \eqn{\mbox{Paran\'{a}}}{Parana} State, Brasil, belonging to the following companies:
  COPEL, IAPAR, DNAEE, SUREHMA  and INEMET.
  
  The data base was organized by Laura Regina Bernardes Kiihl (IAPAR,
  Instituto \eqn{\mbox{Agron\^{o}mico}}{Agronomico} do \eqn{\mbox{Paran\'{a}}}{Parana}, Londrina, Brasil)
  and the fraction of the data included in this data-set was
  provided by Jacinta Loudovico Zamboti (Universidade Estadual de
  Londrina, Brasil).
  The coordinates of the borders of \eqn{\mbox{Paran\'{a}}}{Parana} State were provided
  by \eqn{\mbox{Jo\~{a}o}}{Joao} Henrique Caviglione (IAPAR).
}

\references{
  Diggle, P.J. \& Ribeiro Jr, P.J. (2002) Bayesian inference in
  Gaussian model-based geostatistics. Geographical and Environmental
  Modelling, Vol. 6, No. 2, 129-146.
}

\examples{
summary(parana)
plot(parana, bor=borders)
}
\keyword{datasets}

