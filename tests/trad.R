# carregando o pacote geoR
require(geoR)
par.ori <- par(no.readonly=TRUE)

# Lendo os dados de um arquivo texto no R
cru.df <- read.table("Cruciani.dat", head=T)
cru.df
cru.df$Ksat

## Lendo os dados no formato "geodata" (para uso no pacote geoR)
args(read.geodata)
cru <- read.geodata("Cruciani.dat", head=T, coords=2:3, data.col=4)
cru

## Lendo arquivo com bordas da região
cru.b <- as.matrix(read.table("Cruciani.border", head=T)[,2:3])

## O objeto criado é uma lista
is.list(cru)
names(cru)
cru$coords

## inspecionando e manipulando os elementos da lista.
## Alguns exemplos:
min(cru$coords[,1])
max(cru$coords[,1])
min(cru$coords[,2])
max(cru$coords[,2])

apply(cru$coords,2,range)
max(dist(cru$coords))

apply(cru$coords,2,mean)
apply(cru$coords,2,summary)

summary(cru$data)
hist(cru$data)
boxplot(cru$data)

# um primeiro gráfico com descrição espacial dos dados
plot(cru)
args(plot.geodata)
plot(cru, bor=cru.b)

##
## Verificando tendências nos dados
##

## inpecionando outro conjunto de dados
data(parana)
## vendo informações sobre os dados
#help(parana)

## visualizando os dados
plot(parana)
plot(parana, bord=parana$bor)
## removendo a tendência e visualizando os resíduos
plot(parana, trend="1st", bord=parana$bor)

##
## verificando necessidade de transformação
##
## dados cru são muito assimétricos
## tentando transformações
hist(log(cru$data))
length(cru$data)
boxplot(log(cru$data))
plot(cru, lam=0)

## transformacao BoxCox para dados independentes
require(MASS)
boxcox(cru$data~1)
boxcox(cru$data~1, lam=seq(-1,1,l=21))

## outra função para visualização dos dados
points(cru, lambda=0, bor=cru.b)
points(cru, lambda=0, pt.div="equal", bor=cru.b)
points(cru, lambda=0, pt.div="equal", col=gray(seq(1,0, l=11)), bor=cru.b)
args(points.geodata)

##
## Variogramas
##

## dados originais
cru.vario <- variog(cru)
plot(cru.vario)

## dados transformados (log)
cru.vario <- variog(cru, lam=0)
plot(cru.vario)

## mudando outros argumentos
args(variog)

cru.vario <- variog(cru, lam=0, max.dist=10)
plot(cru.vario)

cru.vario <- variog(cru, lam=0, uvec=seq(0,10,l=11))
plot(cru.vario)

cru.vario <- variog(cru, lam=0, uvec=c(0, 1, 2, 4, 7, 10))
plot(cru.vario)

cru.vario <- variog(cru, lam=0, max.dist=10)
plot(cru.vario)

cru.vario <- variog(cru, lam=0, max.dist=10, trend="1st")
plot(cru.vario)

## variogramas removendo tendência nos dados do Paraná
parana.vario <- variog(parana, max.dist=400)
plot(parana.vario)

parana.vario <- variog(parana, max.dist=400, trend="1st")
plot(parana.vario)

## envelope de variograma para testar a presenca de dependencia espacial
parana.vario.env <- variog.mc.env(parana, obj.v=parana.vario)
plot(parana.vario, env=parana.vario.env)

## agora nos dados de condutividade saturada
cru.vario <- variog(cru, lam=0, max.dist=9)
plot(cru.vario)
cru.vario.env <- variog.mc.env(cru, obj.v=cru.vario)
plot(cru.vario, env=cru.vario.env)


## variogramas direcionais para dados do Parana'
var2 <- variog(parana, max.dist=400, dir=pi/2, trend="1st", uvec=seq(0,400,l=8))
plot(var2, col="blue", ty="l")
var1 <- variog(parana, max.dist=400, dir=0, trend="1st", uvec=seq(0,400,l=8))
lines(var1)
## fazendo várias direções de uma só vez
plot(variog4(parana, max.d=500, trend="1st", uvec=seq(0, 350, l=10)))

##
## "Estimando" parâmetros usando o variograma
##
## variograma escolhido 
cru.vario <- variog(cru, lam=0, max.dist=9)
plot(cru.vario)

## ajustando (estimando parâmetros) "pelo olho"  ("a sentimento")
## Alguns exemplos
args(lines.variomodel.default)
lines.variomodel(seq(0,10,l=100), nug=0.5, cov.pars=c(1.5, 2.5), cov.model="exp", max.dist=10)
lines.variomodel(seq(0,10,l=100), nug=0.5, cov.pars=c(2, 2.5), cov.model="exp", max.dist=10, lty=2)
lines.variomodel(seq(0,10,l=100), nug=0, cov.pars=c(2.3, 2.5), cov.model="exp", max.dist=10, lwd=2)
lines.variomodel(seq(0,10,l=100), nug=0, cov.pars=c(2.3, 7.5), cov.model="sph", max.dist=10, lwd=2, lty=2)

## ajustando (estimando parâmetros) usando mínimos quadrados
args(variofit)
cru.wls <- variofit(cru.vario, ini=c(2.3, 2.5))
cru.wls
lines(cru.wls, col="red", lwd=2)

args(xvalid)
cru.xvm1 <- xvalid(cru, model=cru.wls)
names(cru.xvm1)

cru.ols <- variofit(cru.vario, ini=c(2.3, 2.5), wei="equal", min="optim")
cru.ols
lines(cru.ols, lty=2, col="blue", lwd=3)

## envelope para modelo ajustado
args(variog.model.env)
cru.vario.env.mod <- variog.model.env(cru, model.pars=cru.wls, obj.variog=cru.vario)
plot(cru.vario, env=cru.vario.env.mod)

## ajustando com o modelo Matérn com kappa=2
cru.wls.mat <- variofit(cru.vario, ini=c(2.3, 2.5), cov.model="matern", kappa=2)
cru.wls.mat
lines(cru.wls.mat, col="darkgreen", lwd=2)
##
## Predição espacial 
##
## definindo 2 pontos a serem preditos
loci <- matrix(c(10,15,5,5), ncol=2)
loci

points(cru, lambda=0, pt.div="equal", bor=cru.b)
text(loci)

## usando a função de krigagem (interpolação espacial)
args(krige.conv)
#help(krige.conv)

cru.kc <- krige.conv(cru, loc=loci, krige=krige.control(obj=cru.wls))
cru.kc[1:2]


## Probabilidade de K ser maior que 1.5 em cada um dos pontos
OC <- output.control(thres=1.5, quan=c(.25,.5,.75))
cru.kc <- krige.conv(cru, loc=loci, krige=krige.control(obj=cru.wls), output=OC)
cru.kc$pred
cru.kc$prob
cru.kc$quant
names(cru.kc)
1-cru.kc$prob
dim(cru.kc[[5]])

## predizendo em uma malha de pontos
## definindo um grid regular
points(cru, lambda=0, pt.div="equal", bor=cru.b)

loci0 <- expand.grid(seq(0,24,l=40), seq(-0.5,13, l=30))
points(loci0, pch="+")

loci1 <- polygrid(seq(0,24,l=40), seq(-0.5,13, l=30), cru.b)
points(loci1, col="blue", pch="+")

cru.loci1 <- krige.conv(cru, loc=loci0, bor=cru.b, krige=krige.control(obj=cru.wls))
cru.loci11 <- krige.conv(cru, loc=loci1, krige=krige.control(obj=cru.wls))

names(cru.loci1)
cru.loci1$pred
cru.loci1$krige.var

## mapeando valores preditos
args(image.kriging)
image(cru.loci1)
image(cru.loci11, loc=loci0, bor=cru.b)
## adicionando legenda
image(cru.loci1, x.leg=c(2,21), y.leg=c(-4,-2))
## trocando padrão de cores
image(cru.loci1, col=gray(seq(1,0,l=21)), x.leg=c(2,21), y.leg=c(-4,-2))
## e agora incluindo a localização dos dados
image(cru.loci1, col=gray(seq(1,0,l=21)), x.leg=c(2,21), y.leg=c(-4,-2), coords.dat=cru$coords)

## Outros tipo de gráfico para visualização dos resultados 
contour(cru.loci1)
contour(cru.loci1, fill=T)
persp(cru.loci1)
persp(cru.loci1, theta=20)
persp(cru.loci1, theta=20, phi=30)


# mapeando os erros padrão de predição
image(cru.loci1, col=gray(seq(1,0,l=21)), val=sqrt(cru.loci1$krige.var), x.leg=c(2,21), y.leg=c(-4,-2))
## incluindo as localizações dos dados
image(cru.loci1, col=gray(seq(1,0,l=21)), val=sqrt(cru.loci1$krige.var), x.leg=c(2,21), y.leg=c(-4,-2), coords.dat=cru$coords)


## re-fazendo a krigagem agora com o modelo de Matérn
cru.loci1.mat <- krige.conv(cru, loc=loci0, bor=cru.b, krige=krige.control(obj=cru.wls.mat))
image(cru.loci1.mat, col=gray(seq(1,0,l=21)), coords.dat=cru$coords)

## comparando as predições dos dois modelos
par(mfrow=c(1,2), mar=c(2.5,2.5,0,0), mgp=c(1.2, .5, 0))
## limites comuns para comparar resultados
limis <- range(c(cru.loci1.mat$pred, cru.loci1$pred))
image(cru.loci1, border=cru.b, coords.dat=cru$coords, x.leg=c(2,21), y.leg=c(-4,-2), zlim=limis)
title("modelo exponencial", line=0.5)
image(cru.loci1.mat, col=gray(seq(1,0,l=21)), coords.dat=cru$coords, x.leg=c(2,21), y.leg=c(-4,-2), zlim=limis)
title("modelo Matern", line=0.5)

par(mfrow=c(1,1))
plot(cru.loci1$pred, cru.loci1.mat$pred)
abline(0,1, lwd=2)
##
## Outras opções no output
##
## Gerando simulações em 2 pontos a serem preditos
cru.kc <- krige.conv(cru, loc=loci, krige=krige.control(obj=cru.wls), out=output.control(n.pred=1000))

names(cru.kc)

## visualizando os valores simulados
cru.kc$sim
dim(cru.kc$sim)
hist(cru.kc$sim[1,],prob=T)

## distribuição baseada nas simulações em um ponto de predição
plot(density(cru.kc$sim[1,], bw=0.75), main="Ponto 1")
hist(cru.kc$sim[1,],prob=T,add=T)

## probabilidade de exceder certo valor em cada um dos pontos 
sum(cru.kc$sim[1,] > 3)/1000
sum(cru.kc$sim[2,] > 3)/1000

## opção para calcular prob. de superar valor de referência (threshold)
cru.kc <- krige.conv(cru, loc=loci, krige=krige.control(obj=cru.wls), out=output.control(n.pred=1000, thre=2))
names(cru.kc)
cru.kc$prob

## calculando probabilidades e quantis no grid de predição
cru.loci1 <- krige.conv(cru, loc=loci0, bor=cru.b, krige=krige.control(obj=cru.wls), out=output.control(n.pred=1000, thre=1.5, quant=c(0.10, .5, .9)))
names(cru.loci1)
dim(cru.loci1$sim)
## representando probabilidades em um mapa
par(mfrow=c(1,1))
image(cru.loci1, col=gray(seq(1,0,l=21)), coords.dat=cru$coords, val=1-cru.loci1$prob, x.leg=c(2,21), y.leg=c(-4,-2))

## mapeando  quantis
dim(cru.loci1$quant)
## mapeando o quantil 0.1
image(cru.loci1, col=gray(seq(1,0,l=21)), coords.dat=cru$coords, val=cru.loci1$quan[,1], x.leg=c(2,21), y.leg=c(-4,-2))
## mapeando a mediana
image(cru.loci1, col=gray(seq(1,0,l=21)), coords.dat=cru$coords, val=cru.loci1$quan[,2], x.leg=c(2,21), y.leg=c(-4,-2))
## mapeando o quantil 0.9
image(cru.loci1, col=gray(seq(1,0,l=21)), coords.dat=cru$coords, val=cru.loci1$quan[,3], x.leg=c(2,21), y.leg=c(-4,-2))

##
## Outra forma de visualizar os resultados:
## usando funções do pacote lattice
require(lattice)
## uma imagem "ruim" (sem preservar escala) 
levelplot(cru.loci1$pred ~ loci1$x * loci1$y, col.r=gray(seq(1,0,l=20)))
## fazendo na escala correta
levelplot(cru.loci1$pred ~ loci1$x * loci1$y, col.r=gray(seq(1,0,l=20)),aspect="xy")

## colocando fundo branco
lset()
levelplot(cru.loci1$pred ~ loci1$x * loci1$y, col.r=gray(seq(1,0,l=20)),aspect="xy")

## mudando posição da legenda
levelplot(cru.loci1$pred ~ loci1$x * loci1$y, col.r=gray(seq(1,0,l=20)),aspect="xy", colorkey=list(space="top"))

## vendo outras opções a função
args(levelplot)
#help(levelplot)

## Nota: exportando resultados como arquivo texto
## (para usar resultados em outro software)
## exportando as coordenadas dos pontos de predição
write(t(loci1), file="pp.txt", ncol=2)
## exportando os valores preditos
write(cru.loci1$pred, file="pppred.txt", ncol=1)
## exportando as variâncias de krigagem
write(cru.loci1$krige.var, file="ppvar.txt", ncol=1)

##
## Alternativas para análise:
##     Estimando parâmetros por máxima verossimilhanca
##
cru.ml <- likfit(cru, ini=c(2,1.5), lambda=0)
cru.ml
summary(cru.ml)
#lines(cru.ml, lty=2, col="red", lwd=2)

## estimando com outros valore iniciais para checar
## estabilidade da optimização numérica
cru.ml <- likfit(cru, ini=c(3,0.5), lambda=0)
cru.ml
cru.ml <- likfit(cru, ini=c(1,2.5), lambda=0)
cru.ml
cru.ml <- likfit(cru, ini=c(1,0.5), lambda=0)
cru.ml

## fornecendo diversos valores iniciais
ini.m <- expand.grid(c(1, 1.5, 2, 2.5), c(0,1,2))
ini.m
cru.ml <- likfit(cru, ini=c(2,1.5), lambda=0, nug=c(0, 0.5, 1))
cru.ml

## estimando o parâmetro de transformação
cru.ml <- likfit(cru, ini=c(2,1.5), lambda=0, fix.lam=F)
cru.ml

## estimando com modelo de Matérn com kappa=1
cru.ml.k1 <- likfit(cru, ini=c(2,1.5), lambda=0, cov.model="matern", kappa=1)
cru.ml.k1

## estimando com modelo de Matérn com kappa=2
cru.ml.k2 <- likfit(cru, ini=c(2,1.5), lambda=0, cov.model="matern", kappa=2)
cru.ml.k2

## estimando com modelo esférico
cru.ml.sp <- likfit(cru, ini=c(2, 1.5), lambda=0, cov.model="sph")
cru.ml.sp
lines(cru.ml.sp)

## comparando estimativas dos parâmetros
cru.ml
cru.ml.k1
cru.ml.k2
cru.ml.sp

plot(cru.vario)
lines(cru.ml)
lines(cru.ml.k1, lty=2)
lines(cru.ml.k2, lty=2, lwd=2)
lines(cru.ml.sp, lwd=2)

## estimando um modelo com covariáveis (tendência)
parana.ml <- likfit(parana, ini=c(1200, 100), trend="1st")
parana.ml
summary(parana.ml)

## agora ajustando sem tendência
parana.ml0 <- likfit(parana, ini=c(1200, 100))
summary(parana.ml0)

## comparando os ajustes
summary(parana.ml)
summary(parana.ml0)

## agora ajustando com tendência de 2o grau
parana.ml2 <- likfit(parana, ini=c(1200, 100), trend="2nd")
parana.ml2
summary(parana.ml2)

## comparando os ajustes
summary(parana.ml)
summary(parana.ml2)

## validação cruzada
args(xvalid)
cru.xv <- xvalid(cru, model=cru.ml)
par(mfcol=c(5,2), mar=c(1,1,.5,.5))
plot(cru.xv)
par(par.ori)

## Inferência Bayesiana

cru.mod <- model.control(kappa=1, lambda=0)
cru.mod
cru.prior <- prior.control(phi.dis=seq(0,15, l=16), tausq.rel.prior="unif", tausq.rel.disc=seq(0,1,l=11))
cru.prior

cru.bayes <- krige.bayes(cru, loc=loci, model=cru.mod, prior=cru.prior)
cru.bayes

names(cru.bayes)
names(cru.bayes$post)
names(cru.bayes$post$beta)
cru.bayes$post$beta
cru.bayes$post$sigmasq

cru.bayes$post$joint
dim(cru.bayes$post$joint)

cru.bayes$post$samp

par(mfrow=c(1,2))
plot(cru.bayes)

apply(cru.bayes$post$sampl, 2, mean)

names(cru.bayes$pred)
cru.bayes$pred$mean
cru.bayes$pred$variance
cru.bayes$pred$dist
cru.bayes$pred$med
cru.bayes$pred$variabil
cru.bayes$pred$sim
dim(cru.bayes$pred$sim)

