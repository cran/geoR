require(geoR)
options(digits = 3, width = 80)
set.seed(15)
## =======================================================
##
##  EXAMPLES OF
##        - KTE (KRIGING WITH EXTERNAL TREND)
##          AND
##        - UK (UNIVERSAL KRIGING)
##  USING geoR
##
## =======================================================
##
## 1) KTE
##
## 1.1 Creating an artificial data set with trend on covariates
##
ap <- grf(100, cov.pars=c(1, .25))
ap1.d <- rnorm(100, mean=2)
ap2.d <- runif(100) * 100
ap$data <- 20 - 10 * ap1.d + 2 * ap2.d + ap$data
##
## 1.2 Estimating parameters by ML
##
ap.ml <- likfit(ap, trend = ~ap1.d+ap2.d, ini=c(.5,.5))
ap.ml
##
## 1.3 Defining locations for spatial interpolation
##
pred.grid <- expand.grid((0:20)/20, (0:20)/20)
##
## 1.4 Notice that the external trend information (covariates) must be
## also available at the prediction locations.
## Since this is an artificial example lets create these covariates
##
ap1.l <- rnorm(441, mean=2)
ap2.l <- runif(441) * 100
##
## 1.5 Now kriging with external trend
##
kte.control <- krige.control(trend.d = ~ ap1.d + ap2.d, trend.l = ~ ap1.l + ap2.l,
                             obj.model = ap.ml)
ap.kte <- krige.conv(ap, loc=pred.grid, krige=kte.control)
ap.kte[1:2]
##
## ---------------------------------------------------------------------
##
## 2) UK
##
## Notice that there is no "tecnical" distinction between KTE and UK
## UK is that same os KTE where the covariates are the coordinates
##
## 2.1 Creating an artificial data set with 1st degree trend on the coordinates
##
aq <- grf(100, cov.pars=c(1, .25))
aq$data <- 20 - 10 * aq$coords[,1] + 2 * aq$coords[,2] + aq$data
##
## 2.2 Estimating parameters by ML
##
aq.ml <- likfit(aq, trend = "1st", ini=c(.5,.5))
aq.ml
##
## 3.3 Defining locations for spatial interpolation
##
pred.grid <- expand.grid((0:20)/20, (0:20)/20)
##
## 4.4 Notice that trend information must be
## also provided at the prediction locations.
##
## Now universal kriging with a 1st degree polynomial trend
##
uk.control <- krige.control(trend.d = "1st", trend.l = "1st",
                             obj.model = aq.ml)
aq.uk <- krige.conv(aq, loc=pred.grid, krige=uk.control)
aq.uk[1:2]

