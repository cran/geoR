require(geoR)
options(digits = 3, width = 80)
set.seed(123)
##
## Os comandos abaixo para dados duplicados foram incorporados em as.geodata()
## em geoR_1.2-9
##
ap <- matrix(sample(1:10, 32, rep=T), ncol=2)
ap[2,] <- ap[13,]
ap[4,] <- ap[1,]
ap[8,] <- ap[6,]
ap[14,] <- ap[6,]
ap[10,] <- ap[3,]
ap <- cbind(ap, (1:16), (1:16)*10, (1:16)*100, (1:16)*1000)
ap
as.geodata(ap)
as.geodata(ap, rep.data.ac=mean)
as.geodata(ap, rep.data.ac=median)
as.geodata(ap, rep.data.ac="first")

ap
as.geodata(ap, data.col=3:4)
as.geodata(ap, rep.data.ac=mean, data.col=3:4)
as.geodata(ap, rep.data.ac=median, data.col=3:4)
as.geodata(ap, rep.data.ac="first", data.col=3:4)

ap
as.geodata(ap, data.col=3:4, covar.col=4)
as.geodata(ap, rep.data.ac=mean, data.col=3:4, covar.col=4)
as.geodata(ap, rep.data.ac=median, data.col=3:4, covar.col=4)
as.geodata(ap, rep.data.ac="first", data.col=3:4, covar.col=4)

ap
as.geodata(ap, data.col=3:4, covar.col=4:5)
as.geodata(ap, rep.data.ac=mean, data.col=3:4, covar.col=4:5)
as.geodata(ap, rep.data.ac=median, data.col=3:4, covar.col=4:5)
as.geodata(ap, rep.data.ac="first", data.col=3:4, covar.col=4:5)

ap[2,3] <- NA
ap[6,3] <- NA
ap[8,3] <- NA
ap[9,3] <- NA
ap[16,3] <- NA
ap[10,4] <- NA
ap[11,5] <- NA
ap
as.geodata(ap, data.col=3:4, covar.col=4:5)
