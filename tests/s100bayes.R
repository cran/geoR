## Exemplificando uma análise Bayesiana
require(geoR)
set.seed(345)

data(s100)
summary(s100)
plot(s100)
points(s100)
plot(variog(s100, max.dist=1.1))

## Usando a função krige.bayes()
args(krige.bayes)
args(model.control)
args(prior.control)
args(output.control)

## Definindo o modelo
MC <- model.control()
MC
## Definindo a priori
PC <- prior.control(phi.discrete=seq(0,1,l=11))
PC
## Definindo a resultados a serem retornados
## - quantis (quartis)  .25, .50 e .75
## - Probabilidade da variável ser menor que 1.5
OC <- output.control(n.pos=500, n.pred=100,quantile=c(.25,.5, .75), thres=1.5)
OC
## Definindo o grid de predição
gp <- expand.grid(seq(0,1,l=30), seq(0,1,l=30))
## Rodando a função
s100.kb <- krige.bayes(s100, loc=gp, model=MC, prior=PC, out=OC)
## examinando os resultados
names(s100.kb)
## examinando a posteriori:
names(s100.kb$posterior)
s100.kb$posterior$beta
s100.kb$posterior$sigmasq
s100.kb$posterior$phi
s100.kb$posterior$samples

## examinando a preditiva
names(s100.kb$predictive)
## fazendo um mapa da média
image(s100.kb)
## fazendo um mapa da média
image(s100.kb)
## fazendo um mapa dos erros padrão de predição
image(s100.kb, val=sqrt(s100.kb$pred$variance))
image(s100.kb, val=sqrt(s100.kb$pred$variance), coords.data=s100$coords)

## examinando e mapeando os quartis
s100.kb$pred$quant[1:5,]
## mapa do 1o quartil
image(s100.kb, val=s100.kb$pred$quant[,1])
## mapa do 2o quartil (mediana)
image(s100.kb, val=s100.kb$pred$quant[,2])
## mapa do 3o quartil 
image(s100.kb, val=s100.kb$pred$quant[,3])

## Mapeando a probabilidade P[S < 1.5]
image(s100.kb, val=s100.kb$pred$prob)
## Mapeando a probabilidade P[S > 1.5]
image(s100.kb, val=1-s100.kb$pred$prob)

## inspecionando as simulações da preditiva
dim(s100.kb$pred$sim)
## fazendo mapas de 4 simulações de [S|y]
par.ori <- par(no.readonly=T)
par(mfrow=c(2,2), mar=c(1,1,0,0), mgp=c(1.5,.8,0))
image(s100.kb, val=s100.kb$pred$sim[,1])
image(s100.kb, val=s100.kb$pred$sim[,2])
image(s100.kb, val=s100.kb$pred$sim[,3])
image(s100.kb, val=s100.kb$pred$sim[,4])
par(par.ori)

## Obtendo a distribuição a posterior de [T|y]
## onde T é a porcentagem da área que na qual o
## valor do atributo está acima de 2

## calculando para 1 simulação apenas 
sum(s100.kb$pred$sim[,3] > 2)/900

## agora definindo uma função para aplicar em todas
## as simulações
farea2 <- function(x) sum(x>2)/length(x)
## Aplicando a função para obter a preditiva de [T|y]
t.y <- apply(s100.kb$pred$sim,2,farea2)
hist(t.y)
## obtendo um intervalo de credibilidade a 95%
quantile(t.y, prob=c(0.025, 0.975))

