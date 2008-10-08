\name{soja98}
\alias{soja98}
\docType{data}
\encoding{latin1}
\title{Soya bean production and other variables in a uniformity trial}
\description{
  Data on soya bean production in a uniformity trial measured at
  plots of size 5x5 meters and other soil properties measured
  in points given by the data coordinates.
}

\usage{data(soja98)}

\format{
  A data frame with 256 observations on the following 10 variables.
  \describe{
    \item{\code{X}}{a numeric vector with X-coordinates of the plot centres. }
    \item{\code{Y}}{a numeric vector with X-coordinates of the plot centres. }
    \item{\code{P}}{a numeric vector, phosphorous content. }
    \item{\code{PH}}{a numeric vector, soil pH. }
    \item{\code{K}}{a numeric vector, potassion content. }
    \item{\code{MO}}{a numeric vector, organic matter. }
    \item{\code{SB}}{a numeric vector, basis saturation. }
    \item{\code{iCone}}{a numeric vector, cone index, measuring mechanic
    resistence of the soil.  }
    \item{\code{Rend}}{a numeric vector, total soya production within
      the plot (kg). }
    \item{\code{PROD}}{a numeric vector, production converted to
      productivity by a unit of area - hectare (ton/ha). }
  }
}

%\details{
%  O dados foram produzidos em pesquisa da Universidade Estadual do Oeste do Paran\'a -- Unioeste em a\'area comercial de 17.418,24 m$^2$ de Latossolo vermelho distr\'ofico, com declividade m\'edia de 0,19\%, localizada localizado no Centro de Pesquisa Eloy Gomes, da Cooperativa Central Agropecu\'aria de Desenvolvimento Tecnol\'ogico e Econômico Ltda. (COODETEC), situada na BR 467, km 98, em Cascavel-PR. Nessa \'area, no final do ano de 1997, cultivou-se soja em sistema de semeadura direta. Em abril de 1998, foram demarcadas 256 parcelas de 7,20 $\times$ 7,20 m, medidos os atributos qu\'imicos: pH, Mat\'eria Org\^anica (\%), Pot\'assio ($Cmol_c \;  dm^{-3}$), F\'osforo ($mg \; dm^{-3}$) e Índice de Satura\c{c}\~ao de Bases (\%) e o atributo f\'isico \'Indice de Cone (kg cm$^{-2}$). As amostras para atributos qu\'imicos foram obtidas com 7 cm de di\^ametro e 15 cm de profundidade dentro de cada uma das 256 parcelas, usando-se o sistema desalinhado, sistem\'atico estratificado. Para a Produtividade, foram colhidas e identificadas as parcelas de 5,0 $\times$ 5,0 m, exclu\'idos bordaduras e carreador.
%}

\source{
Souza, E.G.; Jojann, J. A.; Rocha, J. V.; Ribeiro, S. R. A.; Silva, M. S., Upazo, M. A. U.; Molin, J. P.; Oliveira, E. F.; N\'obrega, L. H. P. (1999) Variabilidade espacial dos atributos qu\'imicos do solo em um latossolo roxo distr\'ofico da regi\~ao de Cascavel-PR. \emph{Engenharia Agr\'icola}. Jaboticabal, volume 18, nr. 3, p.80-92.

}
%\references{
%  ~~ possibly secondary sources and usages ~~
%}
\examples{
data(soja98)
plot(soja98)
require(geoR)
prod98 <- as.geodata(soja98, coords.col=1:2, data.col=)
plot(prod98, low=TRUE)
}
\keyword{datasets}
