\name{parana}

\alias{maijun}
\alias{loci.paper}
\alias{borders}

\non_function{mainjun, loci.paper, borders}

\title{Rainfall Data from Parana State, Brasil}

\description{
  This data-set was used by Diggle and Ribeiro (2001) to illustrate the methods
  discussed in the paper. The data reported analysis was
  carried out using the package \bold{geoR}. 
  
  The data refers to average rainfall over different years for the period May-June
  (dry-season). It was collected at 143 recording stations throughout \eqn{\mbox{Paran\'{a}}}{Parana} State,
  Brasil.
}

\usage{data(parana)}

\format{
  Three objects are provided:
  
  \describe{
    \item{\code{maijun}}{a list with coordinates of the recording stations and
      the May-June average data.  }
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
  Diggle, P.J. and Ribeiro Jr, P.J. (2001). Bayesian inference in
  Gaussian model-based geostatistics. \emph{Geographical and Environmental Modeling} (to appear) 
  
  Further information about \bold{geoR} can be found at:\cr
  \url{http://www.maths.lancs.ac.uk/~ribeiro/geoR.html}.
}

\keyword{datasets}
\keyword{spatial}
