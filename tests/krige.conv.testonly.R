require(geoR)
set.seed(123)
options(digits = 3, width = 80)
data(s100)
# Defining a prediction grid
loci <- expand.grid(seq(0,1,l=21), seq(0,1,l=21))
# predicting by ordinary kriging
kc <- krige.conv(s100, loc=loci,
                 krige=krige.control(cov.pars=c(1, .25)))
kc
# mapping point estimates and variances
par.ori <- par(no.readonly = TRUE)
par(mfrow=c(1,2), mar=c(3.5,3.5,1,0), mgp=c(1.5,.5,0))
image(kc, main="kriging estimates")
image(kc, val=sqrt(kc$krige.var), main="kriging std. errors")
# Now setting the output to simulate from the predictive
# (obtaining conditional simulations),
# and to compute quantile and probability estimators
s.out <- output.control(n.predictive = 1000, quant=0.9, thres=2)
kc <- krige.conv(s100, loc = loci,
         krige = krige.control(cov.pars = c(1,.25)),
         output = s.out)
kc
par(mfrow=c(2,2))
image(kc, val=kc$sim[,1], main="a cond. simul.")
image(kc, val=kc$sim[,1], main="another cond. simul.")
image(kc, val=(1 - kc$prob), main="Map of P(Y > 2)")
image(kc, val=kc$quant, main="Map of y s.t. P(Y < y) = 0.9")
par(par.ori)

loci <- expand.grid(seq(0,1,l=6), seq(0,1,l=6))
kc <- krige.conv(s100, loc=loci,
                 krige=list(cov.pars=c(1, .25), kappa = 1))
kc
