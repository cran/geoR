## Exemplo 1: Processo de 0 e 1's
n <- 100
p <- 7/8

u <- numeric(n)
u[1] <- rbinom(1,1,.5)
#for(i in 2:n)
#  u[i] <- ifelse(sum(rbinom(3,1,.5))==3,1-u[i-1],u[i-1])
#for(i in 2:n)
#  u[i] <- ifelse(rbinom(1,1,1-p),1-u[i-1],u[i-1])
for(i in 2:n)
  u[i] <- ifelse(runif(1)<1-p,1-u[i-1],u[i-1])
u
plot(u, pch=19)

pe01 <- function(p, n=100){
  u <- numeric(n)
  u[1] <- rbinom(1,1,.5)
  for(i in 2:n)
    u[i] <- ifelse(runif(1)<1-p,1-u[i-1],u[i-1])
  return(u)
}
sim1 <- pe01(7/8)
plot(sim1, pch=19)

x11()
par(mfrow=c(3,1))
  plot(pe01(7/8), pch=19)
  plot(pe01(3/4), pch=19)
  plot(pe01(1/2), pch=19)

## Exemplo 2: Processo gerando uma variável contínua
pe02 <- function(ro, n=100){
  u <- numeric(n)
  u[1] <- rnorm(1)
  for(i in 2:n)
    u[i] <- ro * u[i-1] + rnorm(1)
  return(u)
}
sim1 <- pe02(0.5)
plot(sim1, pch=19, ty="o")

par(mar=c(3.5,3.5,0,0), mgp=c(2,1,0))

par(mfrow=c(3,1))
for(i in 1:3)
  plot(pe02(0.5), pch=19, ty="o")

x11()
par(mar=c(3.5,3.5,0,0), mgp=c(2,1,0))
par(mfrow=c(3,1))
for(i in 1:3)
  plot(pe02(0.9), pch=19, ty="o", ylim=c(-4,4))

x11()
par(mar=c(3.5,3.5,0,0), mgp=c(2,1,0))
par(mfrow=c(3,1))
for(i in 1:3)
  plot(pe02(0.2), pch=19, ty="o", ylim=c(-4,4))

x11()
par(mar=c(3.5,3.5,0,0), mgp=c(2,1,0))
par(mfrow=c(3,1))
for(i in 1:3)
  plot(pe02(0), pch=19, ty="o", ylim=c(-4,4))

x11()
par(mar=c(3.5,3.5,0,0), mgp=c(2,1,0))
par(mfrow=c(3,1))
  plot(pe02(0.1), pch=19, ty="o", ylim=c(-4,4))
  plot(pe02(0.5), pch=19, ty="o", ylim=c(-4,4))
  plot(pe02(0.9), pch=19, ty="o", ylim=c(-4,4))

## Exemplo 3: processo com mais de 1 parâmetro
pe03 <- function(ro, variancia=1, n=100){
  u <- numeric(n)
  dp <- sqrt(variancia)
  u[1] <- rnorm(1, sd=dp)
  for(i in 2:n)
    u[i] <- ro * u[i-1] + rnorm(1, sd=dp)
  return(u)
}

sim1 <- pe03(0.5)
sim2 <- pe03(0.5, var=2)
sim3 <- pe03(0.5, var=4)
limites <- range(c(sim1,sim2,sim3))

par(mar=c(3.5,3.5,0,0), mgp=c(2,1,0))
par(mfrow=c(3,1))

plot(sim1, ty="o")
plot(sim2, ty="o")
plot(sim3, ty="o")

x11()
par(mar=c(3.5,3.5,0,0), mgp=c(2,1,0))
par(mfrow=c(3,1))
plot(sim1, ty="o", ylim=limites)
plot(sim2, ty="o", ylim=limites)
plot(sim3, ty="o", ylim=limites)

# controlando a semente do gerador de num. aleat.
x11()
par(mar=c(3.5,3.5,0,0), mgp=c(2,1,0))
par(mfrow=c(3,1))
semente <- .Random.seed
plot(pe02(0.1), pch=19, ty="o", ylim=c(-4,4))
.Random.seed <- semente
plot(pe02(0.5), pch=19, ty="o", ylim=c(-4,4))
.Random.seed <- semente
plot(pe02(0.9), pch=19, ty="o", ylim=c(-4,4))
