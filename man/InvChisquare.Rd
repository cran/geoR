\name{InvChisquare}

\alias{dinvchisq}
\alias{rinvchisq}

\title{The (Scaled) Inverse Chi-Squared Distribution}

\description{
  Density and random generation for the scaled inverse chi-squared
  (\eqn{\chi^2_{ScI}}{chi^2_{ScI}}) distribution with
  \code{df} degrees of freedom and optional non-centrality parameter
  \code{scale}.
}
\usage{
dinvchisq(x, df, scale, log = FALSE)
rinvchisq(n, df, scale = 1/df)
}

\arguments{
  \item{x}{vector of quantiles.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length
    is taken to be the number required.}
  \item{df}{degrees of freedom.}
  \item{scale}{scale parameter.}
  \item{log}{logical; if TRUE, densities d are given as log(d).}
}

\value{
  \code{dinvchisq} gives the density and \code{rinvchisq}
  generates random deviates.
}

\details{
  The inverse chi-squared distribution with \code{df}\eqn{= n}
  degrees of freedom has density
  \deqn{f(x) = \frac{1}{{2}^{n/2} \Gamma (n/2)} {(1/x)}^{n/2+1} {e}^{-1/(2x)}}{%
    f(x) = 1 / (2^(n/2) Gamma(n/2))  (1/x)^(n/2-1) e^(-1/(2x))}
  for \eqn{x > 0}.
  The mean and variance are \eqn{\frac{1}{(n-2)}}{1/(n-2)} and
  \eqn{\frac{2}{(n-4)(n-2)^2}}{2/((n-4)(n-2)^2)}.
  
  The non-central chi-squared distribution with \code{df}\eqn{= n}
  degrees of freedom and non-centrality parameter \code{scale}
  \eqn{= S^2} has density
  \deqn{
    f(x) = \frac{{n/2}^{n/2}}{\Gamma (n/2)} S^n {(1/x)}^{n/2+1}
    {e}^{-(n S^2)/(2x)}}{%
    f(x) = ((n/2)^(n/2))/(\Gamma (n/2)) S^n (1/x)^((n/2)+1) e^(-(n
    S^2)/(2x))},
  for \eqn{x \ge 0}.
  The first is a particular case of the latter for \eqn{\lambda = n/2}. 
}

\seealso{
  \code{\link{rchisq}} for the chi-squared distribution which
  is the basis for this function.
}
\examples{
set.seed(1234); rinvchisq(5, df=2)
set.seed(1234); 1/rchisq(5, df=2)

set.seed(1234); rinvchisq(5, df=2, scale=5)
set.seed(1234); 5*2/rchisq(5, df=2)

## inverse Chi-squared is a particular case
x <- 1:10
all.equal(dinvchisq(x, df=2), dinvchisq(x, df=2, scale=1/2))
}
\keyword{distribution}
