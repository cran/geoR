##
## Ilustrando os cálculos dos passos básicos de uma análise geoestatística
##

## 1) Lendo os dados
##    --------------
xy <- read.table("xy.txt")
xy
x <- xy[,1:2]
x
y <- xy[,3]
y
n <- length(y)
n
## especificando ponto a ser predito
x0 <- c(4,4)
x0

## 2) Cálculos da estimação do variograma
##    -----------------------------------
## Calculando matrix de distâncias (euclidianas) entre dados
require(mva)
u <- dist(x)
u

## Calculando a semivariância para cada par de pontos
dist(y)
v <-  as.vector(0.5 * (dist(y)^2))
v

## gráfico da nuvem variográfica
## (note que o gráfico deve incluir a origem) 
plot(u, v, xlim=c(0, max(u)), ylim=c(0, max(v)))

## calculando variograma em classes de distância
## limites das classes de distância: (0, 2, 4, 6, 8)
## pontos médios : (1, 3, 5, 7)
ind <- cut(u, br=c(0,2,4,6,8))
ind
vm <- tapply(v, ind, mean)
um <- c(1,3,5,7)
plot(um, vm, xlim=c(0, max(um)), ylim=c(0, max(vm, na.rm=T)))

## ajustando o modelo gaussiano ("a olho")
plot(function(x) 0.5 + 5 * (1 - exp(-(x/3)^2)), 0, 7, add=T)

## considerando que as estimativas dos parâmetros são:
tau2 = 0.5
sigma2 = 5
phi = 3

## 3) Estimando a média geral
##    -----------------------
## Calculando a matrix R de correlação entre dados
R <- exp(-(u/phi)^2)
R
R <- as.matrix(R)
R
diag(R) <- 1
R
## Matrix de variâncias e covariâncias
Sigma <- tau2 * diag(5) + sigma2 * R

## Estimando a média
um <- rep(1,5)
mu <- solve(t(um) %*% solve(Sigma) %*% um) %*% (t(um) %*% solve(Sigma) %*% y) 
mu

## 4) Predição em um certo ponto (krigagem)
##    ------------------------------------
## Estimando ponto no posição (4,4)
## Calculando distâncias entre cada dado e o ponto a ser predito
u0 <- dist(rbind(x0,x))[1:n]
u0
## Calculando a matrix r de correlação entre dados e ponto a ser predito
r <- exp(-(u0/phi)^2)
r
## Encontrando o valor da estimativa
y0 <- mu +
  (sigma2*t(r))%*%solve(tau2*diag(n)+sigma2*R)%*%(y - mu) 
y0

## e a variância da estimativa
v.y0 <- sigma2 -
   (sigma2*t(r)) %*% solve(tau2*diag(n)+sigma2*R) %*% (sigma2*r) 
v.y0



