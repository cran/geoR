p20 <- read.table("profund020.txt", head=T)
p20
summary(p20)
names(p20)
names(p20) <- tolower(names(p20))
names(p20)
is.factor(p20$tipo)
p20$linha <- as.factor(p20$linha)
p20$ponto <- as.factor(p20$ponto)

p20$x <- (p20$x - min(p20$x))/1000
p20$y <- (p20$y - min(p20$y))/1000

boxplot(p20$k ~ p20$linha)

require(geoR)
## montando o geodata para potássio (K)
k20 <- as.geodata(p20, data.col=10, covar.col=c(3,4,5,6))
summary(k20)
plot(k20)

x11()
plot(k20, trend=~tipodesolo)

x11()
plot(k20, trend=~tipodesolo+linha)

