##
## Simulando do Modelo Linear Generalizado
##
## Exemplo 1: Modelo Binomial (Bernoulli)
##            com ligação "logit"
## valores da covariável:
x <- c(23, 28, 25, 27, 32, 45, 51)
n <- length(x)
## Definindo o valor do parâmetro $\beta$
beta <- .03

## calculando o vetor de médias $\mu$
mu <- exp(x * beta)/(1+ exp(x*beta))
mu
## Simulando dados
y <- rbinom(n, size=1, prob=mu)
y

## Exemplo 2: Modelo Poisson
##            com ligação "log"
##            (usando a mesma covariável)
beta <- 0.05

## calculando o vetor de médias $\mu$
mu <- exp(x * beta)
mu

## Simulando dados
y <- rpois(n, lam=mu)
y

## Exemplo 3: Modelo Poisson com efeito aleatório
##            com ligação "log"
beta <- 0.05
x * beta

## simulando o efeito aleatório
u <- rnorm(n, sd=1)
u
## termo $X\beta + U$
x * beta + u

## calculando o vetor de médias $\mu$
mu <- exp(x * beta + u)

## Simulando dados com efeito aleatório
set.seed(123)
y <- rpois(n, lam=mu)
y

## Mesma simulação sem o efeito aleatório
set.seed(123)
y <- rpois(n, lam=exp(x*beta))
y

## Exemplo 4: Simulando do modelo Poisson
##            com efeito espacial
##            com função de ligação "log"

## definindo as coordenadas dos pontos
cp <- expand.grid(1:10, 1:10)
cp
plot(cp, asp=1)

## simulando valores nestes pontos
require(geoR)
s <- grf(grid=cp, cov.pars=c(2, 3))$data

## calculando a média
mu <- exp(s)

## simulando da Poisson
y <- rpois(100, lam=mu)
plot(cp, asp=1, type="n")
text(cp[,1], cp[,2], y, cex=2)

## visualizando S e Y
par(mfrow=c(1,2), mar=c(3,3,0.5,0.5), mgp=c(2,1,0))
plot(cp, asp=1, type="n")
text(cp[,1], cp[,2], round(s, dig=2), cex=1.5)
plot(cp, asp=1, type="n")
text(cp[,1], cp[,2], y, cex=1.5)

