##
## Comandos no arquivo Rintro.html
##
n1 <- 15
n1

messa <- "Hello you"
messa

n2 <- 1:4
n2

country <- c("France", "Spain", "Italy")
country

ls()
ls(pattern="n")

ls.str()

rm(n1)
ls()

rm(list=ls())
ls()
##
##
##
args(ls)
args(anova)

#help(ls)

#options(htmlhelp = TRUE)  # windows only
#help.start()

#help(ls)
#help(anova)
##
##
##
x1 <- 10
x1

x2 <- c(1, 3, 6)
x2
x2[1]
x2[2]
length(x2)
is.vector(x2)
is.matrix(x2)
is.numeric(x2)
is.character(x2)

x3 <- 1:10
x3

x4 <- seq(0,1, by=0.1)
x4
x4[x4 > 0.5]
x4 > 0.5

x5 <- seq(0,1,  l=11)
x5

x6 <- rep(1, 5)
x6

x7 <- rep(c(1, 2), c(3, 5))
x7

x8 <- rep(1:3, rep(5,3))
x8

x9 <- rnorm(10, mean=70, sd=10)
x9
sum(x9)
mean(x9)
var(x9)
min(x9)
max(x9)
summary(1:10)

x10 <- x9[x9 > 72]

args(seq)
args(rep)
##
##
##
m1 <- matrix(1:12, ncol=3)
m1
length(m1)
dim(m1)
nrow(m1)
ncol(m1)
m1[1,2]
m1[2,2]
m1[,2]
m1[3,]
dimnames(m1)
dimnames(m1) <- list(c("L1", "L2", "L3","L4"), c("C1","C2","C3"))
dimnames(m1)
m1[c("L1","L3"),]
m1[c(1,3),]

m2 <- cbind(1:5, 6:10)
m2

m3 <- cbind(1:5, 6)
m3

ar1 <- array(1:24, dim=c(3,4,2))
ar1[,2:3,]
ar1[2,,1]
sum(ar1[,,1])
sum(ar1[1:2,,1])
##
##
##
d1 <- data.frame(X = 1:10, Y = c(51, 54, 61, 67, 68, 75, 77, 75, 80, 82))
d1
names(d1)
d1$X
d1$Y
plot(d1)
plot(d1$X, d1$Y)

d2 <- data.frame(Y= c(10+rnorm(5, sd=2), 16+rnorm(5, sd=2), 14+rnorm(5, sd=2)))
d2$lev <- gl(3,5)
d2
by(d2$Y, d2$lev, summary)

d3 <- expand.grid(1:3, 4:5)
##
##
##
l1 <- list(A=1:10, B="THIS IS A MESSAGE", C=matrix(1:9, ncol=3))
l1

l2 <- lm(Y ~ X, data=d1)
l2
is.list(l2)
class(l2)
summary(l2)
anova(l2)
names(l2)
l2$pred
l2$res
par(mfrow=c(2,2))
plot(l2)
par(mfrow=c(1,1))
l3 <- aov(Y ~ lev, data=d2)
l3
summary(l3)
##
##
##
lm
plot
plot.default

min
max
lines
##
##
##    
cru <- read.table("Cruciani.dat")
cru
##
##
##
library(geoR)
#library(geoR, lib="PATH_TO_YOUR_geoR_DIRECTORY")

cru <- read.geodata("Cruciani.dat", head=T, coords.col=2:3, data.col=4)
cru.borda <- read.table("Cruciani.border", head=T)[,2:3]
cru.borda <- rbind(cru.borda, cru.borda[1,])

cru
is.list(cru)
class(cru)

plot(cru, bord=cru.borda)
 
args(plot.geodata)
#help.start()
#help(plot.geodata)

plot(cru, bord=cru.borda, lambda=0)

points(cru, bord=cru.borda)
args(points.geodata)
points(cru, bord=cru.borda, cex.min=1, cex.max=1, pt.div="quartile")
points(cru, bord=cru.borda, cex.min=1, cex.max=1,col=gray(seq(0.9,0,l=length(cru$data))))
points(cru, bord=cru.borda, cex.min=1, cex.max=1, col=gray(seq(0.9,0,l=length(cru$data))), xla="Coord X", ylab="Coord Y")

points(cru, lambda=0, bord=cru.borda)
points(cru, lambda=0, bord=cru.borda, cex.min=1, cex.max=1, pt.div="quartile")
points(cru, lambda=0, bord=cru.borda, cex.min=1, cex.max=1, col=gray(seq(0.9,0,l=length(cru$data))))

hist(cru$data, main="", xlab="K")
hist(log(cru$data), main="", xlab="log(K)")
##
##
##
data(wolfcamp)
wolfcamp
#help(wolfcamp)

data(s100)

data(parana)
