##
## Aspectos iniciais de análise de dados utilizando pacote geoR 
## ============================================================
##
##
## 0. Comandos iniciais (ajustando o ambiente do R)
## ------------------------------------------------
##
## Mostrar ajuda no "browser" (opcional)
#  options(htmlhelp = TRUE)
## Salvando parâmetros gráficos originais (opcional) 
par.ori <- par(no.readonly = TRUE)
## Carregando o pacote geoR
require(geoR)
## library(geoR)
## library(geoR, lib="coloque_diretorio_onde_geoR_foi_instalado_aqui")
##
##
## 1. Entrada de dados
## -------------------
##
## 1.1. Carregando dados incluídos no pacote 
##
ls()
data(s100)
ls()
#help(s100)

data(SIC)
ls()
#help(SIC)
##
## 1.2. Convertendo data-frames para objetos no formato geodata
##
require(MASS)
data(topo)
topo
gtopo <- as.geodata(topo)
gtopo

args(as.geodata)
#help(as.geodata)
##
## 1.3. Lendo dados de um arquivo ASCII (texto)
##
## exemplo: copiar os arquivos Cruciani.dat e Cruciani.border
##          para o diretório dados sob o de trabalho.
##          (veja o link "dados" na sessão de tutoriais da página do curso)
##
## lendo dados
cru <- read.geodata("Cruciani.dat", head=T, coords.col=2:3, data.col=4)
cru
## lendo arquivo com bordas da área
cru.borda <- read.table("Cruciani.border", head=T)[,2:3]
cru.borda <- rbind(cru.borda, cru.borda[1,])
##
##
## 2. Análise exploratária I - visualisando os dados
##    (usando funções/métodos do pacote geoR)
##
plot(cru, bord = cru.borda)

## Inspecting the options and documentation for the plot options
args(plot.geodata)
#help(plot.geodata)

plot(cru, bord=cru.borda, lambda=0)

hist(cru$data, main="", xlab="K")
hist(log(cru$data), main="", xlab="log(K)")

points(cru, bord=cru.borda)

args(points.geodata)
#help(points.geodata)

points(cru, bord=cru.borda, cex.min=1, cex.max=1, pt.div="quartile")
points(cru, bord=cru.borda, cex.min=1, cex.max=1, col=gray(seq(0.9,0,l=length(cru$data))))
points(cru, bord=cru.borda, cex.min=1, cex.max=1, col=gray(seq(0.9,0,l=length(cru$data))), xla="Coord X", ylab="Coord Y")

points(cru, lambda=0, bord=cru.borda)
points(cru, lambda=0, bord=cru.borda, cex.min=1, cex.max=1, pt.div="quartile")
points(cru, lambda=0, bord=cru.borda, cex.min=1, cex.max=1, col=gray(seq(0.9,0,l=length(cru$data))))

##
##
## 3. Análise exploratória 2 - variogramas
##
cru.v1 <-  variog(cru)
plot(cru.v1)
cru.v1
names(cru.v1)

args(variog)
#help(variog)

cru.cl1 <- variog(cru, option = "cloud")
plot(cru.cl1)

cru.v2 <- variog(cru, lambda = 0)
plot(cru.v2)
cru.cl2 <- variog(cru, lambda = 0, option = "cloud")
plot(cru.cl2)

cru.v3 <- variog(cru, lambda = 0, max.dist=9)
plot(cru.v3)

cru.v4 <- variog(cru, lambda = 0, uvec=seq(0,9,l=8))
plot(cru.v4)

cru.v5 <- variog(cru, lambda = 0, uvec=seq(0,12,l=8))
plot(cru.v5)

cru.v6 <- variog(cru, lambda = 0, uvec=seq(0,15,l=8))
plot(cru.v6)


cru.v <- variog(cru, lambda = 0, uvec=seq(0,9,l=8), bin.cl=T)
plot(cru.v)
plot(cru.v, bin.cloud = T)

cru.mv <- variog(cru, lambda = 0, uvec=seq(0,9,l=8), est="mod")
plot(cru.mv)

##
## Sugestões:
##  - verifique outros argumentos para serem utilizados/explorados na função variog()
##  - veja tambem a função variog4()

cru.v <- variog(cru, lambda = 0, uvec=seq(0,9,l=8), bin.cl=T)
cru.env <- variog.mc.env(cru, obj.variog = cru.v)
plot(cru.v, env = cru.env) 

##
## Sugestão:
##  - explore outros conjuntos de dados incluídos no pacote
##
#data(package=geoR)

data(wolfcamp)
#help(wolfcamp)

data(parana)
#help(parana)


