## Ilustrando:
##  - a montagem de uma covariável definida como a distância
##      da cada ponto a uma certa localidade
##  - a obtenção de IC para coeficiente de uma covariável
require(geoR)
## carregando um conjunto de dados
cru <- read.geodata("Cruciani.dat", head=T, coords.col=2:3, data.col=4)
cru.borda <- read.table("Cruciani.border", head=T)[,2:3]
cru.borda <- rbind(cru.borda, cru.borda[1,])
points(cru)

## marcando a localidade no mapa
#pt <- locator()
#pt
#pt <- unlist(pt)
#pt
pt <- c(14,9)

## montando a covariável:
##   calculando a distância da localidade a cada ponto
d1 <- dist(rbind(pt, cru$coords))[1:32]

## variograma após "remoção" da tendência
variog(cru, trend=~d1)

## ajustando o modelo com a covariável
ml <- likfit(cru, ini=c(2, 3), trend=~d1)
ml

## vendo as estimativas dos parâmetros de médias, seus erros padrão
summary(ml)
names(ml)
ml$beta
ml$beta.var
## montando o IC (95%) para o coeficiente da covariável 
ml$beta[2] +  qnorm(c(0.025, 0.975)) * sqrt(ml$beta.var[2,2])

