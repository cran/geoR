## Definindo parâmetros
mu <- 5
sigma2 <- 3
phi <- 2 # função de correlação exponencial
tau2 <- 0.5

## Dados
y <- c(10,7,8)
y
n <- length(y)
n
## Coordenadas dos dados
x <- matrix(c(1,4,3,4,3,1), nc=2)
x
## Coordenada do ponto a ser predito
x0 <- c(2,2)
x0
## Calculando distâncias entre dados
d <- dist(x)
d
## Calculando a matrix R de correlação entre dados
R <- exp(-d/phi)
R
R <- as.matrix(R)
R
diag(R) <- 1
R
## Calculando distâncias entre cada dado e o ponto a ser predito
d0 <- dist(rbind(x0,x))[1:n]
d0
## Calculando a matrix r de correlação entre dados e ponto a ser predito
r <- exp(-d0/phi)
r
## Encontrando o valor da estimativa
y0 <- mu +
  (sigma2*t(r))%*%solve(tau2*diag(n)+sigma2*R)%*%(y - mu) 
y0

## e a variância da estimativa
v.y0 <- sigma2 -
   (sigma2*t(r)) %*% solve(tau2*diag(n)+sigma2*R) %*% (sigma2*r) 
v.y0
