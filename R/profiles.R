"proflik" <- 
  function(obj.likfit, geodata, coords = geodata$coords,
           data = geodata$data,
           sill.values,
           range.values, 
           nugget.values,
           nugget.rel.values,
           lambda.values,
           sillrange.values = TRUE,
           sillnugget.values = TRUE,
           rangenugget.values = TRUE, 
           sillnugget.rel.values = FALSE,
           rangenugget.rel.values = FALSE, 
           silllambda.values = FALSE,
           rangelambda.values = TRUE, 
           nuggetlambda.values = FALSE,
           nugget.rellambda.values = FALSE,
           uni.only = TRUE,
           bi.only = FALSE,
           minimisation.function = c("optim", "nlmP"), ...)
{
  ##
  ## 1. setting arguments
  ##
  require(mva)
  call.fc <- match.call()
  minimisation.function <- match.arg(minimisation.function)
  n.cov.pars <- obj.likfit$npars - length(obj.likfit$beta)
  if(obj.likfit$transform.info$fix.lambda == FALSE)
    n.cov.pars <- n.cov.pars - 1
  if(missing(sill.values)) sill.values <- FALSE
  if(missing(range.values)) range.values <- FALSE 
  if(!is.null(obj.likfit$call$fix.nugget))
    if(obj.likfit$call$fix.nugget == TRUE)
      nugget.values <-  nugget.rel.values <- FALSE
  if(missing(nugget.values)) nugget.values <- FALSE
  if(missing(nugget.rel.values)) nugget.rel.values <- FALSE
  if(missing(lambda.values) | obj.likfit$transform.info$fix.lambda == TRUE)
    lambda.values <- FALSE  
  if(uni.only == TRUE){
    sillrange.values <- sillnugget.values <- rangenugget.values <-
      sillnugget.rel.values <- rangenugget.rel.values <- 
        silllambda.values <- rangelambda.values <- 
          nugget.rellambda.values <- nuggetlambda.values <- FALSE
  }
  else{
    if(all(sillrange.values == TRUE)){
      if(all(sill.values == FALSE) | all(range.values == FALSE)){
        sillrange.values <- FALSE
        stop("if argument sillrange.values = TRUE sill.values and range.values must be provided. Alternatively a matrix can be provided in sillrange.values  or set this to FALSE")
      }
      else
        sillrange.values <- as.matrix(expand.grid(sill.values, range.values))
    }
    if(n.cov.pars == 2){
        sillnugget.values <- rangenugget.values <-
          sillnugget.rel.values <- rangenugget.rel.values <- 
            nugget.rellambda.values <- nuggetlambda.values <- FALSE
    }
    else{
      if(all(sillnugget.values == TRUE)){
        if(all(sill.values == FALSE) | all(nugget.values == FALSE)){
          sillnugget.values <- FALSE
          stop("if argument sillnugget.values = TRUE sill.values and nugget.values must be provided. Alternatively a matrix can be provided in sillnugget.values or set this to FALSE")
        }
        else
          sillnugget.values <- as.matrix(expand.grid(sill.values, nugget.values))
      }
      if(all(rangenugget.values == TRUE)){
        if(all(range.values == FALSE) | all(nugget.values == FALSE)){
          rangenugget.values <- FALSE
          stop("if argument rangenugget.values = TRUE range.values and nugget.values must be provided. Alternatively a matrix can be provided in rangenugget.values or set this to FALSE")
        }
        else
          rangenugget.values <- as.matrix(expand.grid(range.values, nugget.values))
      }
      if(all(sillnugget.rel.values == TRUE)){
        if(all(sill.values == FALSE) | all(nugget.rel.values == FALSE)){
          sillnugget.rel.values <- FALSE
          stop("if argument sillnugget.rel.values = TRUE sill.values and nugget.rel.values must be provided. Alternatively a matrix can be provided in sillnugget.rel.values or set this to FALSE")
        }
        else
          sillnugget.rel.values <- as.matrix(expand.grid(sill.values, nugget.rel.values))
      }
      if(all(rangenugget.rel.values == TRUE)){
        if(all(range.values == FALSE) | all(nugget.rel.values == FALSE)){
          rangenugget.rel.values <- FALSE
          stop("if argument rangenugget.rel.values = TRUE range.values and nugget.rel.values must be provided. Alternatively a matrix can be provided in rangenugget.rel.values or set this to FALSE")
        }
        else
          rangenugget.rel.values <- as.matrix(expand.grid(range.values, nugget.rel.values))
      }
      if(obj.likfit$transform.info$fix.lambda == TRUE){
        if(all(nuggetlambda.values == TRUE)){
          if(all(lambda.values == FALSE) | all(nugget.values == FALSE)){
            nuggetlambda.values <- FALSE
            stop("if argument nuggetlambda.values = TRUE lambda.values and nugget.values must be provided. Alternatively a matrix can be provided in nuggetlambda.values or set this to FALSE")
          }
          else
            nuggetlambda.values <- as.matrix(expand.grid(lambda.values, nugget.values))
        }
        if(all(nugget.rellambda.values == TRUE)){
          if(all(lambda.values == FALSE) | all(nugget.rel.values == FALSE)){
            nugget.rellambda.values <- FALSE
            stop("if argument nugget.rellambda.values = TRUE lambda.values and nugget.rel.values must be provided. Alternatively a matrix can be provided in nugget.rellambda.values or set this to FALSE")
          }
          else
            nugget.rellambda.values <- as.matrix(expand.grid(lambda.values, nugget.rel.values))
        }
      }
    }
    if(obj.likfit$transform.info$fix.lambda == TRUE)
      silllambda.values <- rangelambda.values <- FALSE
    else{
      if(all(silllambda.values == TRUE)){
        if(all(sill.values == FALSE) | all(lambda.values == FALSE)){
          silllambda.values <- FALSE
          stop("if argument silllambda.values = TRUE sill.values and lambda.values must be provided. Alternatively a matrix can be provided in silllambda.values or set this to FALSE")
        }
        else
          silllambda.values <- as.matrix(expand.grid(sill.values, lambda.values))
      }
      if(all(rangelambda.values == TRUE)){
        if(all(range.values == FALSE) | all(lambda.values == FALSE)){
          rangelambda.values <- FALSE
          stop("if argument rangelambda.values = TRUE range.values and lambda.values must be provided. Alternatively a matrix can be provided in rangelambda.values or set this to FALSE")
        }
        else
          rangelambda.values <- as.matrix(expand.grid(range.values, lambda.values))
      }      
    }
  }
  ##
  ## 2. data preparation
  ##
  trend <- trend.spatial(trend=obj.likfit$trend, coords=coords)
  data <- as.vector(data)
  dimnames(trend) <- list(NULL, NULL)
  if(obj.likfit$transform.info$fix.lambda == TRUE) {
    if(obj.likfit$lambda != 1) {
      if(any(data <= 0))
        stop("Data transformation not allowed when there are zeros or negative data"
             )
      if(obj.likfit$lambda == 0)
        data <- log(data)
      else data <- ((data^obj.likfit$lambda) - 1)/obj.likfit$
      lambda
    }
  }
  n <- length(data)
  dists.vec <- as.vector(dist(coords))
  d <- range(dists.vec)
  min.dist <- d[1]
  max.dist <- d[2]
  tausq <- obj.likfit$nugget
  sigmasq <- obj.likfit$cov.pars[1]
  tausq.rel <- tausq/sigmasq
  phi <- obj.likfit$cov.pars[2]
  lambda <- obj.likfit$lambda
  loglik <- obj.likfit$loglik
  sill.total <- sigmasq + tausq
  n.uni <- 0
  n.bi <- 0
  lower.phi <- 0.01 * (min.dist/max.dist) 
  upper.phi <- 1000 * max.dist
  lower.sigmasq <- 0.01 * sill.total
  result <- list()
  assign(".temp.list", list(n = n,
                            z = data,
                            beta.size = dim(trend)[2],
                            kappa = obj.likfit$kappa,
                            xmat = trend,
                            ## txmat = t(trend),
                            method = obj.likfit$method,
                            dists.lowertri = dists.vec,
                            cov.model = obj.likfit$cov.model,
                            fix.lambda = obj.likfit$transform.info$fix.lambda,
                            lambda = obj.likfit$lambda,
                            lower.phi = lower.phi,
                            upper.phi = upper.phi,
                            lower.sigmasq = lower.sigmasq, 
                            phi.est = phi,
                            tausq.rel.est = tausq.rel,
                            tausq.est = tausq,
                            sigmasq.est = sigmasq,
                            minimisation.function = minimisation.function), pos=1)
  if(obj.likfit$transform.info$fix.lambda == TRUE) {
    .temp.list$log.jacobian <<- obj.likfit$transform.info$log.jacobian
  }
  ##
  ## 3. One-dimentional profile likelihoods
  ##
  ##  
  ## 3.1 Profile for sigmasq
  ##  
  if(bi.only == FALSE) {
    if(any(sill.values != FALSE)) {
      n.uni <- n.uni + 1
      cat("proflik: computing profile likelihood for the sill\n")
      if(n.cov.pars == 2) {
        if(tausq == 0) {
          if(obj.likfit$transform.info$fix.lambda == FALSE) {
            ini.grid <- as.matrix(expand.grid(seq(min(range.values),
                                                  max(range.values),
                                                            l = 5),
                                                  seq(-1,1,l = 5)))
          }
          else {
            ini.grid <- as.matrix(
                                            seq(min(range.values),
                                                max(range.values),
                                                l = 10))
          }
          dimnames(ini.grid) <- list(NULL, NULL)
          .temp.list$ini.grid <<- ini.grid
          pl.sigmasq <- apply(matrix(sill.values,
                                     ncol = 1), 1, proflik.aux2, ...)
          .temp.list$ini.grid <<- NULL
        }
        else {
          stop("not yet implemented for fixed nugget != 0"
               )
        }
      }
      if(n.cov.pars == 3) {
        if(any(lambda.values != FALSE)) {
          ini.grid <- as.matrix(expand.grid(seq(min(range.values),
                                                max(range.values), l = 6),
                                            seq(0, 2 * tausq.rel, l = 4),
                                            seq(-1, 1, l = 5)))
        }
        else {
          ini.grid <- as.matrix(expand.grid(seq(min(range.values),
                                                max(range.values), l = 10),
                                            seq(0, 2 * tausq.rel, l = 4)))
        }
        dimnames(ini.grid) <- list(NULL, NULL)
        .temp.list$ini.grid <<- ini.grid
        pl.sigmasq <- apply(matrix(sill.values, ncol = 
                                   1), 1, proflik.aux9, ...)
        .temp.list$ini.grid <<- NULL
      }
      v.ord <- order(c(sigmasq, sill.values))
      result$sill <- list(sill = c(sigmasq, sill.values)[
                            v.ord], proflik.sill = c(loglik, pl.sigmasq)[
                                      v.ord], est.sill = c(sigmasq, loglik))
    }
    ##  
    ## 3.2 Profile for phi
    ##
    if(any(range.values != FALSE)) {
      n.uni <- n.uni + 1
      cat("proflik: computing profile likelihood for the range\n")
      if(n.cov.pars == 2) {
        if(tausq == 0) {
          .temp.list$nugget <<- 0
          pl.phi <- apply(matrix(range.values,
                                 ncol = 1), 1, proflik.aux0, ...)
          .temp.list$nugget <<- NULL
        }
        else {
          stop("not yet implemented for fixed nugget != 0"
               )
        }
      }
      if(n.cov.pars == 3) {
        pl.phi <- apply(matrix(range.values, ncol = 1),
                        1, proflik.aux7, ...)
      }
      v.ord <- order(c(phi, range.values))
      result$range <- list(range = c(phi, range.values)[
                             v.ord], proflik.range = c(loglik, pl.phi)[
                                       v.ord], est.range = c(phi, loglik))
    }
    ##  
    ## 3.3 Profile for \tau^2
    ##  
    if(n.cov.pars == 3) {
      if(any(nugget.values != FALSE)) {
        n.uni <- n.uni + 1
        cat("proflik: computing profile likelihood for the nugget\n"
              )
        pl.tausq <- apply(matrix(nugget.values, ncol = 
                                 1), 1, proflik.aux11, ...)
        v.ord <- order(c(tausq, nugget.values))
        result$nugget <- list(nugget = c(tausq, 
                                nugget.values)[v.ord], proflik.nugget
                              = c(loglik, pl.tausq)[v.ord], 
                              est.nugget = c(tausq, loglik))
      }
      ##  
      ## 3.4 Profile for relative \tau^2
      ##
      if(any(nugget.rel.values != FALSE)) {
        cat("proflik: computing profile likelihood for the relative nugget\n"
              )
        n.uni <- n.uni + 1
        pl.tausq.rel <- apply(matrix(nugget.rel.values,
                                     ncol = 1), 1, proflik.aux5, ...)
        v.ord <- order(c(tausq.rel, nugget.rel.values))
        result$nugget.rel <- list(nugget.rel = c(
                                    tausq.rel, nugget.rel.values)[v.ord],
                                  proflik.nugget = c(loglik, pl.tausq.rel
                                    )[v.ord], est.nugget.rel = c(tausq.rel,
                                                loglik))
      }
    }
    ##  
    ## 3.5 Profile for \lambda
    ##
    if(any(lambda.values != FALSE)) {
      .temp.temp.list <<- .temp.list
      .temp.temp.list$coords <<- coords
      n.uni <- n.uni + 1
      cat("proflik: computing profile likelihood for the transformation parameter\n"
            )
      if(n.cov.pars == 2) {
        if(tausq == 0) {
          .temp.temp.list$fixtau <<- T
          .temp.temp.list$ini <<- c(sigmasq,phi)
          pl.lambda <- apply(as.matrix(lambda.values), 1,
                             proflik.aux23, ...)
        }
        else {
          stop("not yet implemented for fixed nugget != 0"
               )
        }
      }
      if(n.cov.pars == 3) {
        .temp.temp.list$fixtau <<- FALSE
        .temp.temp.list$ini <<- phi
        pl.lambda <- apply(matrix(lambda.values,
                                  ncol = 1), 1, proflik.aux23, ...)
      }
      v.ord <- order(c(lambda, lambda.values))
      result$lambda <- list(lambda = c(lambda, 
                              lambda.values)[v.ord], proflik.lambda
                            = c(loglik, pl.lambda)[v.ord], 
                            est.lambda = c(lambda, loglik))
      remove(.temp.temp.list, inherits=TRUE, pos=1)
    }
  }
  ##
  ## 4. Two-dimentional profile likelihoods
  ##
  ##  
  ## 4.1 Profile for \sigma^2 and \phi
  ##
  if(uni.only == FALSE){
    if(any(sillrange.values != FALSE)) {
      n.bi <- n.bi + 1
      cat("proflik: computing 2-D profile likelihood for the sill and range parameters\n")
      if(n.cov.pars == 2) {
        if(tausq == 0) {
          .temp.list$nugget <<- 0
          if(.temp.list$fix.lambda == TRUE) {
            pl.sigmasqphi <- apply(cbind(0, sillrange.values, 1), 1, loglik.spatial, ...)
          }
          else {
            pl.sigmasqphi <- apply(sillrange.values,
                                   1, proflik.aux28, ...)
          }
          .temp.list$nugget <<- NULL
        }
        else {
          stop("not yet implemented for fixed nugget != 0"
               )
        }
      }
      if(n.cov.pars == 3) {
        pl.sigmasqphi <- apply(sillrange.values, 1, proflik.aux13, ...)
      }
      names(pl.sigmasqphi) <- NULL
      result$sillrange <- list(sill = as.numeric(levels(as.factor(sillrange.values[,1]))), range = 
                               as.numeric(levels(as.factor(sillrange.values[,2]))), proflik.sillrange = pl.sigmasqphi, 
                               est.sillrange = c(sigmasq, phi, loglik))
    }
    ##  
    ## 4.2 Profile for \sigma^2 and \tau^2
    ##  
    if(any(sillnugget.values != FALSE)) {
      n.bi <- n.bi + 1
      cat("proflik: computing 2-D profile likelihood for the sill and nugget\n")
      if(obj.likfit$transform.info$fix.lambda == FALSE)
        ini.grid <- as.matrix(expand.grid(seq(min(range.values),
                                              max(range.values), l = 
                                              10), seq(-1, 1, l = 5)))
      else
        ini.grid <- as.matrix(seq(min(range.values),
                                  max(range.values), l = 10))
      dimnames(ini.grid) <- list(NULL, NULL)
      .temp.list$ini.grid <<- ini.grid
      pl.sigmasqtausq <- apply(sillnugget.values, 1, proflik.aux15, ...)
      .temp.list$ini.grid <<- NULL
      names(pl.sigmasqtausq) <- NULL
      result$sillnugget <- list(sill = as.numeric(levels(as.factor(sillnugget.values[,1]))), nugget = as.numeric(levels(as.factor(sillnugget.values[,2]))), proflik.sillnugget = 
                                pl.sigmasqtausq, est.sillrange = c(sigmasq,
                                                   tausq, loglik))
    }
    ##  
    ## 4.3 Profile for \phi and \tau^2
    ##
    if(any(rangenugget.values != FALSE)) {
      n.bi <- n.bi + 1
      cat("proflik: computing 2-D profile likelihood for the range and nugget\n"
            )
      .temp.list$ini.grid <<- as.matrix(seq(sigmasq/4, 5 * 
                                          sigmasq, l = 15))
      pl.phitausq <- apply(rangenugget.values, 1, proflik.aux17, ...)
      .temp.list$ini.grid <<- NULL
      names(pl.phitausq) <- NULL
      result$rangenugget <- list(range = as.numeric(levels(as.factor(rangenugget.values[,1]))), nugget
                               = as.numeric(levels(as.factor(rangenugget.values[,1]))), proflik.rangenugget = 
                               pl.phitausq, est.rangenugget = c(phi, tausq,
                                              loglik))
    }
    ##  
    ## 4.4 Profile for \sigma^2 and \tau^2_{rel}
    ##
    if(any(sillnugget.rel.values != FALSE)) {
      n.bi <- n.bi + 1
      cat("proflik: computing 2-D profile likelihood for the sill and relative nugget parameters\n"
            )
      if(.temp.list$fix.lambda == FALSE)
        ini.grid <- as.matrix(expand.grid(seq(min(range.values), max(range.values), l = 
                                              10), seq(-1, 1, l = 5)))
      else
        ini.grid <- as.matrix(seq(min(range.values),
                                  max(range.values), l = 10))
      dimnames(ini.grid) <- list(NULL, NULL)
      .temp.list$ini.grid <<- ini.grid
      pl.sigmasqtausq.rel <- apply(sillnugget.rel.values, 1, 
                                   proflik.aux19, ...)
      .temp.list$ini.grid <<- NULL
      names(pl.sigmasqtausq.rel) <- NULL
      result$sillnugget.rel <- list(sill = as.numeric(levels(as.factor(sillnugget.rel.values[,1]))), 
                                    nugget.rel = as.numeric(levels(as.factor(sillnugget.rel.values[,2]))), 
                                    proflik.sillnugget.rel = pl.sigmasqtausq.rel,
                                    est.sillrange.rel = c(sigmasq, tausq.rel, 
                                      loglik))
    }
    ##  
    ## 4.5 Profile for \phi and \tau^2_{rel}
    ##
    if(any(rangenugget.rel.values != FALSE)) {
      n.bi <- n.bi + 1
      cat("proflik: computing 2-D profile likelihood for the range and relative nugget parameters\n"
            )
      pl.phitausq.rel <- apply(rangenugget.rel.values, 1, proflik.aux30, ...)
      names(pl.phitausq.rel) <- NULL
      result$rangenugget.rel <- list(range = as.numeric(levels(as.factor(rangenugget.rel.values[,1]))),
                                   nugget.rel = as.numeric(levels(as.factor(rangenugget.rel.values[,2]))), 
                                   proflik.rangenugget.rel = pl.phitausq.rel,
                                   "est.rangenugget .rel" = c(phi, tausq.rel,
                                     loglik))
    }
  }
  ##  
  ## 4.6 Profile for \sigma^2 and \lambda
  ##
  if(any(silllambda.values != FALSE)) {
    n.bi <- n.bi + 1
    cat("proflik: computing 2-D profile likelihood for the sill and transformation parameters\n"
          )
    if(n.cov.pars == 2) {
      ini.grid <- as.matrix(seq(min(range.values), max(
                                                       range.values), l = 10))
      dimnames(ini.grid) <- list(NULL, NULL)
      .temp.list$ini.grid <<- ini.grid
      if(tausq == 0) {
        .temp.list$nugget <<- 0
        pl.sigmasqlambda <- apply(silllambda.values, 1,
                                  proflik.aux24, ...)
        .temp.list$ini.grid <<- .temp.list$nugget <<- NULL
      }
      else {
        stop("not yet implemented for fixed nugget != 0"
             )
      }
    }
    if(n.cov.pars == 3) {
      ini.grid <- as.matrix(expand.grid(seq(min(range.values),
                                            max(range.values), l = 10), seq(0, 1, l = 5)))
      dimnames(ini.grid) <- list(NULL, NULL)
      .temp.list$ini.grid <<- ini.grid
      pl.sigmasqlambda <- apply(sigmasqlambda.values, 1, 
                                proflik.aux27, ...)
      .temp.list$ini.grid <<- NULL
    }
    names(pl.sigmasqlambda) <- NULL
    result$silllambda <- list(sill = as.numeric(levels(as.factor(silllambda.values[,1]))), lambda = as.numeric(levels(as.factor(silllambda.values[,2]))), proflik.silllambda = pl.sigmasqlambda,
                              est.silllambda = c(sigmasq, lambda, loglik))
  }
  ##  
  ## 4.7 Profile for \phi and \lambda
  ##
  if(any(rangelambda.values != FALSE)) {
    .temp.list$data <<- .temp.list$z
    n.bi <- n.bi + 1
    cat("proflik: computing 2-D profile likelihood for the range and transformation parameters\n"
              )
    if(n.cov.pars == 2) {
      if(tausq == 0) {
        .temp.list$nugget <<- 0
        pl.philambda <- apply(rangelambda.values, 1, 
                              proflik.aux1, ...)
        .temp.list$nugget <<- NULL
      }
      else {
        stop("not yet implemented for fixed nugget != 0"
             )
      }
    }
    if(n.cov.pars == 3) {
      pl.philambda <- apply(rangelambda.values, 1, proflik.aux31, ...)
    }
    names(pl.philambda) <- NULL
    result$rangelambda <- list(range = as.numeric(levels(as.factor(rangelambda.values[,1]))), lambda = as.numeric(levels(as.factor(rangelambda.values[,2]))), proflik.rangelambda = pl.philambda,
                               est.rangelambda = c(phi, lambda, loglik))
  }
  ##  
  ## 4.8 Profile for \tau^2 and \lambda
  ##                                        
  if(any(nuggetlambda.values != FALSE)) {
    n.bi <- n.bi + 1
    cat("proflik: computing 2-D profile likelihood for the nugget and transformation parameters\n"
          )
    pl.nuggetlambda <- apply(nuggetlambda.values, 1, proflik.aux32, ...)
      names(pl.nuggetlambda) <- NULL
    result$nuggetlambda <- list(nugget = as.numeric(levels(as.factor(nuggetlambda.values[,1]))), lambda = as.numeric(levels(as.factor(nuggetlambda.values[,2])))
                                , proflik.nuggetlambda = pl.nuggetlambda,
                                est.nuggetlambda = c(tausq, lambda, loglik))
  }
  ##  
  ## 4.9 2-D Profile for \tau^2_{rel} and \lambda
  ##
  if(any(nugget.rellambda.values != FALSE)) {
    n.bi <- n.bi + 1
    pl.nugget.rellambda <- apply(nugget.rellambda.values, 1, proflik.aux33, ...)
    names(pl.nugget.rellambda) <- NULL
    result$nugget.rellambda <- list(nugget.rel = as.numeric(levels(as.factor(nugget.rellambda.values[,1]))),
                                    lambda = as.numeric(levels(as.factor(nugget.rellambda.values[,2]))), proflik.nugget.rellambda = 
                                    pl.nugget.rellambda, est.nugget.rellambda = c(tausq.rel,
                                                           lambda, loglik))
  }
  result$n.uni <- n.uni
  result$n.bi <- n.bi
  result$method <- obj.likfit$method
  result$call <- call.fc
  class(result) <- "proflik"
  return(result)
}

"proflik.aux0" <-
  function(phi, ...)
{
  ## This function computes the value of the profile likelihood for the correlation parameter \phi when nugget effect is not included in the model.
  ## It requires the minimisation of the function wrt \phi for each value of \lambda, if this transformation parameter is included in the model
  ## This is an auxiliary function called by likfit.proflik
  ##
  if(.temp.list$fix.lambda == TRUE)
    proflik <- proflik.aux1(phi = phi)
  else {
    .temp.list$phi <<- phi
    if(.temp.list$minimisation.function == "optim")
      proflik <-  - (optim(.temp.list$lambda, proflik.aux1.1, method="L-BFGS-B", lower
                           = -2, upper = 2, ...)$value)
    else
      proflik <-  - (nlmP(proflik.aux1.1, .temp.list$lambda, lower
                          = -2, upper = 2, ...)$minimum)
    .temp.list$phi <<- NULL
  }
  return(proflik)
}
"proflik.aux1" <-
  function(philambda, ...)
{
  ## This function computes the value of the profile likelihood for the correlation function scale parameter \phi when nugget effect = 0
  if(length(philambda) == 2) lambda <- philambda[2]
  else lambda <- 1
  n <- .temp.list$n
  main <- proflik.main(tausq=.temp.list$nugget, sigmasq=1, phi=philambda[1], lambda = lambda)
  if(.temp.list$method == "ML") {
    proflik <-  - (n/2) * log(2 * pi) - main$log.det.to.half -
      (n/2) * log(main$ssresmat/n) - (n/2) + main$
    log.jacobian
  }
  if(.temp.list$method == "RML") {
    eigentrem <- eigen(main$ixix, symmetric = TRUE, only.values = TRUE)
    proflik <-  - ((n - beta.size)/2) * log(2 * pi) - main$
    log.det.to.half - ((n - beta.size)/2) * log(main$ssresmat/
                                                n) - (n/2) + 0.5 * sum(log(eigentrem$values)) + 
                                                  main$log.jacobian
  }
  return(proflik)
}
"proflik.aux10" <-
  function(phitausq.rel.lambda, ...)
{
  if(length(phitausq.rel.lambda) == 3) lambda <- phitausq.rel.lambda[3]
  else lambda <- 1
  phitausq.rel.lambda <- as.vector(phitausq.rel.lambda)
  n <- .temp.list$n
  phi <- phitausq.rel.lambda[1]
  tausq <- phitausq.rel.lambda[2]
  sigmasq <- .temp.list$sigmasq
  main <- proflik.main(tausq=tausq, sigmasq=1, phi=phi, lambda = lambda)
  if(.temp.list$method == "ML") {
    neglik <- (n/2) * log(2 * pi) + main$log.det.to.half +
      (n/2) * log(sigmasq) + (0.5/sigmasq) * main$ssresmat -
        main$log.jacobian
  }
  if(.temp.list$method == "RML") {
    eigentrem <- eigen(main$ixix, symmetric = TRUE, only.values = TRUE)
    neglik <- ((n - beta.size)/2) * log(2 * pi) + main$
    log.det.to.half + ((n - beta.size)/2) * log(sigmasq) +
      (0.5/sigmasq) * main$ssresmat - 0.5 * sum(log(eigentrem$
                                               values)) - main$log.jacobian
  }
  return(as.vector(round(neglik, dig=8)))
}
"proflik.aux11" <-
  function(tausq, ...)
{
  ## This function computes the value of the profile likelihood for the parameter \tau^2.
  ## It requires the minimisation of the function wrt \sigma^2, \phi and \lambda (if the case)  for each value of \tau^2.
  ## This is an auxiliary function called by proflik.
  .temp.list$nugget <<- as.vector(tausq)
  if(.temp.list$fix.lambda == TRUE) {
    if(.temp.list$minimisation.function == "optim"){
      sigmasqphi.res <- optim(c(.temp.list$sigmasq.est, .temp.list$phi.est),
                              proflik.aux12,method="L-BFGS-B",
                              lower = c(.temp.list$lower.sigmasq,
                                .temp.list$lower.phi),
                              upper=c(+Inf, .temp.list$upper.phi), ...)$value
    }
    else{
      sigmasqphi.res <- nlmP(proflik.aux12,
                             c(.temp.list$sigmasq.est, .temp.list$phi.est),
                             lower = c(.temp.list$lower.sigmasq,
                               .temp.list$lower.phi),
                             upper=c(+Inf, .temp.list$upper.phi), ...)$minimum
    }
  }
  else {
    if(.temp.list$minimisation.function == "optim"){
      sigmasqphi.res <- optim(c(.temp.list$sigmasq.est, .temp.list$
                                phi.est, .temp.list$lambda), proflik.aux12,method="L-BFGS-B",  lower = c(.temp.list$lower.sigmasq, .temp.list$lower.phi, -2),
                              upper = c( + Inf, .temp.list$upper.phi, 2), ...)$value
    }
    else{
      sigmasqphi.res <- nlmP(proflik.aux12,c(.temp.list$sigmasq.est, .temp.list$
                                             phi.est, .temp.list$lambda),
                             lower = c(.temp.list$lower.sigmasq, .temp.list$lower.phi, -2),
                             upper = c( + Inf, .temp.list$upper.phi, 2), ...)$minimum
    }    
  }
  .temp.list$nugget <<- NULL
  return( - sigmasqphi.res)    
}

"proflik.aux1.1" <-
  function(lambda, ...)
{
  ## This function computes the value of the profile likelihood for the correlation function scale parameter \phi when nugget effect = 0
  phi <- .temp.list$phi
  n <- .temp.list$n
  main <- proflik.main(tausq=.temp.list$nugget, sigmasq=1, phi=phi, lambda = lambda)
  if(.temp.list$method == "ML") {
    neglik <- (n/2) * log(2 * pi) + main$log.det.to.half +
      (n/2) * log(main$ssresmat/n) + (n/2) - main$log.jacobian
  }
  if(.temp.list$method == "RML") {
    eigentrem <- eigen(main$ixix, symmetric = TRUE, only.values = TRUE)
    neglik <- ((n - beta.size)/2) * log(2 * pi) + main$
    log.det.to.half + ((n - beta.size)/2) * log(main$ssresmat/n) +
      (n/2) - 0.5 * sum(log(eigentrem$values)) - 
        main$log.jacobian
  }
  return(as.vector(round(neglik, dig=8)))
}

"proflik.aux12" <-
  function(sigmasqphi.lambda, ...)
{
  ## This function computes the value of the profile likelihood for the nugget parameter \tau^2, minimizing the likelihood wrt correlation function scale parameter \phi (range), the random field scale parameter \sigma^2 (sill) and the transformation parameter \lambda. 
  if(length(sigmasqphi.lambda) == 3) lambda <-  sigmasqphi.lambda[3]
  else lambda <- 1
  sigmasqphi.lambda <- as.vector(sigmasqphi.lambda)
  n <- .temp.list$n
  sigmasq <- sigmasqphi.lambda[1]
  phi <- sigmasqphi.lambda[2]
  main <- proflik.main(tausq=.temp.list$nugget, sigmasq=sigmasq, phi=phi, lambda = lambda)
  if(.temp.list$method == "ML") {
    neglik <- (n/2) * log(2 * pi) +
      main$log.det.to.half +
        0.5 * (main$ssresmat) -
          main$log.jacobian
  }
  if(.temp.list$method == "RML") {
    eigentrem <- eigen(main$ixix, symmetric = TRUE, only.values = TRUE)
    neglik <- ((n - beta.size)/2) * log(2 * pi) +
      main$log.det.to.half +
        0.5 * (main$ssresmat) -
          0.5 * sum(log(eigentrem$values)) -
            main$log.jacobian
  }
  return(as.vector(round(neglik, dig=8)))
}
"proflik.aux13" <-
  function(sigmasqphi, ...)
{
  ## This function computes the value of the 2-D profile likelihood for the parameters \sigma^2 and \phi when the nugget is included.
  ## It requires the minimisation of the function wrt \tau^2 and \lambda (if the case) for each value of (\sigma^2, \phi)
  ## This is an auxiliary function called by likfit.proflik
  .temp.list$sigmasqphi <<- as.vector(sigmasqphi)
  if(.temp.list$fix.lambda == TRUE) {
    if(.temp.list$minimisation.function == "optim") 
      tausq.res <- optim(.temp.list$tausq.est, proflik.aux14, method="L-BFGS-B", lower
			 = 0, ...)$value
    else
      tausq.res <- nlmP(proflik.aux14, .temp.list$tausq.est, lower
                        = 0, ...)$minimum
  }
  else {
    if(.temp.list$minimisation.function == "optim") 
      tausq.res <- optim(
                         c(.temp.list$tausq.est, .temp.list$lambda), proflik.aux14, method="L-BFGS-B",lower = c(0, -2
                                                                                                        ), upper = c( +Inf, 2), ...)$value
    else
      tausq.res <- nlmP(
                        proflik.aux14, c(.temp.list$tausq.est, .temp.list$lambda),lower = c(0, -2
                                                                                    ), upper = c( +Inf, 2), ...)$minimum
  }
  .temp.list$sigmasqphi <<- NULL
  return( - tausq.res)
}

"proflik.aux14" <-
  function(tausq.lambda, ...)
{
  ## This function computes the value of the 2-D profile likelihood for the parameters (\sigma^2, \phi), minimizing the likelihood wrt the nugget parameter \tau^2.
  ## This functions is called by the auxiliary function proflik.aux13
  if(length(tausq.lambda) == 2) lambda <- tausq.lambda[2]
  else lambda <- 1
  n <- .temp.list$n
  tausq <- tausq.lambda[1]
  main <- proflik.main(tausq=tausq, .temp.list$sigmasqphi[1], phi=.temp.list$sigmasqphi[2], lambda = lambda)
  if(.temp.list$method == "ML") {
    neglik <- (n/2) * log(2 * pi) + main$log.det.to.half + 0.5 *
      main$ssresmat - main$log.jacobian
  }
  if(.temp.list$method == "RML") {
    eigentrem <- eigen(main$ixix, symmetric = TRUE, only.values = TRUE)
    neglik <- ((n - beta.size)/2) * log(2 * pi) + main$
    log.det.to.half + 0.5 * main$ssresmat - 0.5 * sum(log(
                                                     eigentrem$values)) - main$log.jacobian
  }
  return(as.vector(round(neglik, dig=8)))
}
"proflik.aux15" <-
  function(sigmasqtausq, ...)
{
  ## This function computes the value of the 2-D profile likelihood for the parameters \sigma^2 and \tau^2
  ## It requires the minimisation of the function wrt \phi and also \lambda (if the case) for each value of (\sigma^2, \tau^2) 
  ## This is an auxiliary function called by likfit.proflik
  .temp.list$sigmasqtausq <<- as.vector(sigmasqtausq)
  if(.temp.list$fix.lambda == TRUE) {
    ini.lik <- round(100000000. * apply(.temp.list$ini.grid, 1,
                                        proflik.aux16))
    ini <- as.vector(.temp.list$ini.grid[ini.lik == min(ini.lik),
                                         ])
    if(.temp.list$minimisation.function == "optim") 
      phi.res <- optim(ini, proflik.aux16, method="L-BFGS-B", lower = 
                       .temp.list$lower.phi, upper=.temp.list$upper.phi, ...)$value
    else
      phi.res <- nlmP(proflik.aux16, ini, lower = 
                      .temp.list$lower.phi, upper=.temp.list$upper.phi, ...)$minimum
  }
  else {
    ini.lik <- round(100000000. * apply(.temp.list$ini.grid, 1,
                                        proflik.aux16))
    ini <- as.vector(.temp.list$ini.grid[ini.lik == min(ini.lik),
                                         ])
    if(.temp.list$minimisation.function == "optim") 
      phi.res <- optim(ini, proflik.aux16, method="L-BFGS-B", 
                       lower = c(.temp.list$lower.phi, -2),
                       upper = c(.temp.list$upper.phi, 2), ...)$value
    else
      phi.res <- nlmP(proflik.aux16, ini, 
                      lower = c(.temp.list$lower.phi, -2),
                      upper = c(.temp.list$upper.phi, 2), ...)$minimum
  }
  .temp.list$sigmasqtausq <<- NULL
  return( - phi.res)
}

"proflik.aux16" <-
  function(phi.lambda, ...)
{
  ## This function computes the value of the 2-D profile likelihood for the sill and nugget parameters (\sigma^2,\tau^2), minimising the profile likelihood wrt correlation function scale parameter \phi (and the transformation parameter \lambda
  ## This is an auxiliary function called by likfit.aux15  
  if(length(phi.lambda) == 2) lambda <- phi.lambda[2]
  else lambda <- 1
  n <- .temp.list$n
  phi <- phi.lambda[1]
  main <- proflik.main(tausq=.temp.list$sigmasqtausq[2], sigmasq=.temp.list$sigmasqtausq[1], phi=phi, lambda = lambda)
  if(.temp.list$method == "ML") {
    neglik <- (n/2) * log(2 * pi) + main$log.det.to.half + 0.5 *
      main$ssresmat - main$log.jacobian
  }
  if(.temp.list$method == "RML") {
    eigentrem <- eigen(main$ixix, symmetric = TRUE, only.values = TRUE)
    neglik <- ((n - beta.size)/2) * log(2 * pi) + main$
    log.det.to.half + 0.5 * main$ssresmat - 0.5 * sum(log(
                                                     eigentrem$values)) -
                                                       main$log.jacobian
  }
    return(as.vector(round(neglik, dig=8)))
}
"proflik.aux17" <-
  function(phitausq, ...)
{
  ## This function computes the value of the 2-D profile likelihood for the parameters (\phi, \tau^2)
  ## It requires the minimisation of the function wrt \sigma^2 and \lambda (if the case) for each value of (\phi, \tau^2) 
  ## This is an auxiliary function called by likfit.proflik
  .temp.list$phitausq <<- as.vector(phitausq)
  if(.temp.list$fix.lambda == TRUE) {
    ini.lik <- round(100000000. * apply(.temp.list$ini.grid, 1,
                                        proflik.aux18))
    ini <- as.vector(.temp.list$ini.grid[ini.lik == min(ini.lik),])
    if(.temp.list$minimisation.function == "optim") 
      sigmasq.res <- optim(ini, proflik.aux18, method="L-BFGS-B", 
                           lower = .temp.list$lower.sigmasq, ...)$value
    else
      sigmasq.res <- nlmP(proflik.aux18, ini, 
                          lower = .temp.list$lower.sigmasq, ...)$minimum
  }
  else {
    if(.temp.list$minimisation.function == "optim") 
      sigmasq.res <- optim(c(.temp.list$sigmasq.est, .temp.list$lambda
                             ), proflik.aux18, method="L-BFGS-B", lower = c(.temp.list$lower.sigmasq,
                                                                    -2), upper = c( + Inf, 2), ...)$value
    else
      sigmasq.res <- nlmP(proflik.aux18, c(.temp.list$sigmasq.est, .temp.list$lambda
                                           ), lower = c(.temp.list$lower.sigmasq,
                                                -2), upper = c( + Inf, 2), ...)$minimum
  }
  .temp.list$phitausq <<- NULL
  return( - sigmasq.res)
}

"proflik.aux18" <-
  function(sigmasq.lambda, ...)
{
  ## This function computes the value of the 2-D profile likelihood for the range and nugget parameters (\phi, \tau^2), minimising the likelihood wrt the random field scale parameter \sigma^2 (sill) ant the transformation parameter \lambda. 
  ## This is an auxiliary function called by likfit.aux17.
  if(length(sigmasq.lambda) == 2) lambda <- sigmasq.lambda[2]
  else lambda <- 1
  n <- .temp.list$n
  sigmasq <- sigmasq.lambda[1]
  main <- proflik.main(tausq=.temp.list$phitausq[2], sigmasq=sigmasq, phi=.temp.list$phitausq[1], lambda = lambda)
  if(.temp.list$method == "ML") {
    neglik <- (n/2) * log(2 * pi) + main$log.det.to.half + 0.5 *
      main$ssresmat - main$log.jacobian
  }
  if(.temp.list$method == "RML") {
    eigentrem <- eigen(main$ixix, symmetric = TRUE, only.values = TRUE)
    neglik <- ((n - beta.size)/2) * log(2 * pi) + main$
    log.det.to.half + 0.5 * main$ssresmat - 0.5 * sum(log(
                                                     eigentrem$values)) -
                                                       main$log.jacobian
  }
  return(as.vector(round(neglik, dig=8)))
}
"proflik.aux19" <-
  function(sigmasqtausq.rel, ...)
{
  ## This function computes the value of the 2-D profile likelihood for the parameters (\sigma^2, \tau^2_{rel})
  ## It requires the minimisation of the function wrt \phi and \lambda (if the case) for each value of (\sigma^2, \tau^2_{rel})
  ## This is an auxiliary function called by likfit.proflik
  .temp.list$sigmasqtausq.rel <<- as.vector(sigmasqtausq.rel)
  if(.temp.list$fix.lambda == TRUE) {
    if(.temp.list$minimisation.function == "optim") 
      phi.res <- optim(.temp.list$phi.est, proflik.aux20, method="L-BFGS-B", lower = 
                       .temp.list$lower.phi, upper=.temp.list$upper.phi, ...)$value
    else
      phi.res <- nlmP(proflik.aux20, .temp.list$phi.est, lower = 
                      .temp.list$lower.phi, upper=.temp.list$upper.phi, ...)$minimum
  }
  else {
    if(.temp.list$minimisation.function == "optim") 
      phi.res <- optim(c(.temp.list$phi.est, .temp.list$lambda, ...), proflik.aux20, method="L-BFGS-B", 
                       lower = c(.temp.list$lower.phi, -2),
                       upper = c(.temp.list$upper.phi, 2), ...)$value
    else
      phi.res <- nlmP(proflik.aux20, c(.temp.list$phi.est, .temp.list$lambda, ...), 
                      lower = c(.temp.list$lower.phi, -2),
                      upper = c(.temp.list$upper.phi, 2), ...)$minimum
  }
  .temp.list$sigmasqtausq.rel <<- NULL
  return( - phi.res)
}
"proflik.aux2" <-
  function(sigmasq, ...)
{
  ## This function computes the value of the profile likelihood for the random field scale (variance) parameter \sigma^2 when nugget effect is not included in the model.
  ## It requires the minimisation of the function wrt \phi for each value of \sigma^2
  ## This is an auxiliary function called by likfit.proflik
  ##
  .temp.list$sigmasq <<- as.vector(sigmasq)
  if(.temp.list$fix.lambda == TRUE) {
    ini.lik <- round(100000000. * apply(.temp.list$ini.grid, 1,
                                        proflik.aux3))
    ini <- as.vector(.temp.list$ini.grid[ini.lik == min(ini.lik),
                                         ])
    if(.temp.list$minimisation.function == "optim") 
      phi.res <- optim(ini , proflik.aux3, method="L-BFGS-B", lower = .temp.list$lower.phi,upper=.temp.list$upper.phi, ...)$value
    else
      phi.res <- nlmP(proflik.aux3, ini , lower = .temp.list$lower.phi,upper=.temp.list$upper.phi, ...)$minimum
  }
  else {
    ini.lik <- round(100000000. * apply(.temp.list$ini.grid, 1,
                                        proflik.aux3))
    ini <- as.vector(.temp.list$ini.grid[ini.lik == min(ini.lik),])
    if(.temp.list$minimisation.function == "optim") 
      phi.res <- optim(ini, proflik.aux3, method="L-BFGS-B",
                       lower = c(.temp.list$lower.phi, -2),
                       upper = c(.temp.list$upper.phi, 2), ...)$value
    else
      phi.res <- nlmP(proflik.aux3, ini,
                      lower = c(.temp.list$lower.phi, -2),
                      upper = c(.temp.list$upper.phi, 2), ...)$minimum
  }
  .temp.list$sigmasq <<- NULL
  return( - phi.res)
}

"proflik.aux20" <-
  function(phi.lambda, ...)
{
  ## This function computes the value of the 2-D profile likelihood for the sill and relative nugget parameters (\sigma^2, \tau^2_{rel}), minimising the likelihood wrt the correlation function scale parameter \phi and the transformation parameter \lambda.
  ## This is an auxiliary function called by likfit.aux19.
  phi.lambda <- as.vector(phi.lambda)
  if(length(phi.lambda) == 2) lambda <- phi.lambda[2]
  else lambda <- 1
  sigmasqtausq.rel <- as.vector(.temp.list$sigmasqtausq.rel)
  sigmasq <- sigmasqtausq.rel[1]
  tausq.rel <- sigmasqtausq.rel[2]
  phi <- phi.lambda[1]
  n <- .temp.list$n
  main <- proflik.main(tausq=tausq.rel, sigmasq=1, phi=phi, lambda = lambda)
  if(.temp.list$method == "ML") {
    neglik <- (n/2) * log(2 * pi) +
      main$log.det.to.half +
        (n/2) * log(sigmasq) +
          (0.5/sigmasq) * main$ssresmat -
            main$log.jacobian
  }
  if(.temp.list$method == "RML") {
    eigentrem <- eigen(main$ixix, symmetric = TRUE, only.values = TRUE)
    neglik <- ((n - beta.size)/2) * log(2 * pi) +
      main$log.det.to.half +
        (n/2) * log(sigmasq) +
          (0.5/sigmasq) * main$ssresmat -
            0.5 * sum(log(eigentrem$values)) -
              main$log.jacobian
  } 
  return(as.vector(round(neglik, dig=8)))
}
"proflik.aux21" <-
function(phitausq.rel, ...)
{
  ## This function computes the value of the 2-D profile likelihood for the parameters (\phi, \tau^2_{rel})
  ## This is an auxiliary function called by likfit.proflik
  phitausq.rel <- as.vector(phitausq.rel)
  phi <- phitausq.rel[1]
  tausq.rel <- phitausq.rel[2]
  n <- .temp.list$n
  main <- proflik.main(tausq=tausq.rel, sigmasq=1, phi=phi, lambda = 1)
  sigmasq.hat <- main$ssresmat/n
  if(.temp.list$method == "ML") {
    proflik <-  - (n/2) * log(2 * pi) -
      main$log.det.to.half -
      (n/2) * log(sigmasq.hat) -
        (n/2) -
          main$log.jacobian
  }
  if(.temp.list$method == "RML") {
    eigentrem <- eigen(main$ixix, symmetric = TRUE, only.values = TRUE)
    proflik <-  - ((n - beta.size)/2) * log(2 * pi) -
      main$log.det.to.half -
        (n/2) * log(sigmasq.hat) -
          (n/2) +
            0.5 * sum(log(eigentrem$values)) -
              main$log.jacobian
  }
  return(proflik)
}

"proflik.aux21.1" <-
  function(lambda, ...)
{
  ## This function computes the value of the 2-D profile likelihood for the parameters (\phi, \tau^2_{rel})
  ## This requires minimasation wrt to the transformation parameter \lambda
  ## This is an auxiliary function called by likfit.proflik
  n <- .temp.list$n
  main <- proflik.main(tausq = .temp.list$phitausq.rel[2], sigmasq = 1,
                       phi = .temp.list$phitausq.rel[1], lambda = lambda)
  sigmasq.hat <- main$ssresmat/n
  if(.temp.list$method == "ML") {
    neglik <- (n/2) * log(2 * pi) +
      main$log.det.to.half + (n/2) * log(sigmasq.hat) +
        (n/2) -
          main$log.jacobian
  }
  if(.temp.list$method == "RML") {
    eigentrem <- eigen(main$ixix, symmetric = TRUE, only.values = TRUE)
    neglik <- ((n - beta.size)/2) * log(2 * pi) +
      main$log.det.to.half +
        (n/2) * log(sigmasq.hat) +
          (n/2) -
            0.5 * sum(log(eigentrem$values)) -
              main$log.jacobian
  }
  return(as.vector(round(neglik, dig=8)))
}

"proflik.aux22" <-
  function(sigmasq, ...)
{
  ## This function computes the value of the 2-D profile likelihood for the range and nugget parameters (\phi, \tau^2), minimising the likelihood wrt the random field scale parameter \sigma^2 (sill) 
  ## This is an auxiliary function called by likfit.aux17
  n <- .temp.list$n
  main <- proflik.main(tausq=.temp.list$phitausq[2], sigmasq=sigmasq, phi= .temp.list$phitausq[1], lambda = 1)
  if(.temp.list$method == "ML") {
    neglik <- (n/2) * log(2 * pi) +
      main$log.det.to.half +
        0.5 * main$ssresmat -
          main$log.jacobian
  }
  if(.temp.list$method == "RML") {
    eigentrem <- eigen(main$ixix, symmetric = TRUE, only.values = TRUE)
    neglik <- ((n - beta.size)/2) * log(2 * pi) +
      main$log.det.to.half +
        0.5 * main$ssresmat -
          0.5 * sum(log(eigentrem$values)) -
            main$log.jacobian
  }
  return(as.vector(round(neglik, dig=8)))
}

"proflik.aux23" <-
  function(lambda, ...)
{
  ## This function computes the value of the profile likelihood for the transformation parameter \lambda
  ## It requires the minimisation of the function wrt \phi and \tau^2 and sigma^2 for each value of \lambda
  ## This is an auxiliary function called by proflik
  lambda <- as.vector(lambda)
  if(.temp.temp.list$fixtau == FALSE) {
    if(lambda == 0)
      data.l <- log(.temp.list$z)
    else data.l <- ((.temp.list$z^lambda) - 1)/lambda
    var.l <- var(data.l)
    ini.l <- c(0.10000000000000001 * var.l, var.l, .temp.temp.list$ini)
  }
  else
    ini.l <- .temp.temp.list$ini
  if(dim(.temp.list$xmat)[2] == 1 & all(.temp.list$xmat == 1))
    trend.mat <- "cte"
  else
    trend.mat <- ~ (.temp.list$xmat[,-1])
  lambda.res <- likfit(coords = .temp.temp.list$coords, data = .temp.list$z,
                       ini = ini.l, trend = trend.mat, fix.nugget = 
                       .temp.temp.list$fixtau,
                       method = .temp.list$method, cov.model = 
                       .temp.list$cov.model, kappa = .temp.list$kappa, fix.lambda = TRUE,
                       lambda = lambda, messages.screen = FALSE)$loglik
  .temp.list <<- .temp.temp.list
  return(lambda.res)
}

"proflik.aux24" <-
  function(sigmasqlambda, ...)
{
  ## This function computes the value of the 2-D profile likelihood for the parameters (\sigma^2, \lambda) when there is no nugget effect (\tau^2 = 0, fixed)
  ## It requires the minimisation of the function wrt \phi for each value of (\sigma^2, \lambda)
  ## This is an auxiliary function called by proflik
  sigmasqlambda <- as.vector(sigmasqlambda)
  .temp.list$sigmasq <<- sigmasqlambda[1]
  lambda <- sigmasqlambda[2]
  if(lambda == 1) {
    .temp.list$log.jacobian <<- 0
  }
  else {
    if(any(.temp.list$z <= 0))
      stop("Transformation option not allowed when there are zeros or negative data"
           )
    .temp.list$log.jacobian <<- sum(log(.temp.list$z^(lambda - 1)))
    if(lambda == 0)
      .temp.list$z <<- log(.temp.list$z)
    else .temp.list$z <<- ((.temp.list$z^lambda) - 1)/lambda
  }
  ini.lik <- round(100000000. * apply(.temp.list$ini.grid, 1,
                                      proflik.aux3))
  ini <- as.vector(.temp.list$ini.grid[ini.lik == min(ini.lik),
                                       ])
  if(.temp.list$minimisation.function == "optim") 
    phi.res <- optim(ini, proflik.aux3, method="L-BFGS-B", lower = .temp.list$
                     lower.phi, upper = .temp.list$upper.phi, ...)$value
  else
    phi.res <- nlmP(proflik.aux3, ini, lower = .temp.list$
                    lower.phi, upper = .temp.list$upper.phi, ...)$minimum
  .temp.list$log.jacobian <<- NULL
  .temp.list$sigmasq <<- NULL
  .temp.list$z <<- .temp.list$data
  return( - phi.res)
}

"proflik.aux27" <-
  function(sigmasqlambda, ...)
{
  ## This function computes the value of the 2-D profile likelihood for sill \sigma^2 and the transformation parameter \lambda
  ## It requires the minimisation of the function wrt \phi and \tau^2 and for each value of (\sigma^2,\lambda)
  ## This is an auxiliary function called by proflik.
  sigmasqlambda <- as.vector(sigmasqlambda)
  .temp.list$sigmasq <<- sigmasqlambda[1]
  lambda <- sigmasqlambda[2]
  if(lambda == 1) {
    .temp.list$log.jacobian <<- 0
  }
  else {
    .temp.list$fix.lambda <<- T
    if(any(.temp.list$z^(lambda - 1) <= 0))
      .temp.list$log.jacobian <<- log(prod(.temp.list$z^(lambda -
                                                         1)))
    else .temp.list$log.jacobian <<- sum(log(.temp.list$z^(lambda -
                                                           1)))
    if(lambda == 0)
      .temp.list$z <<- log(.temp.list$z)
    else .temp.list$z <<- ((.temp.list$z^lambda) - 1)/lambda
  }
  ini.lik <- round(100000000. * apply(.temp.list$ini.grid, 1,
                                      proflik.aux10))
  ini <- as.vector(.temp.list$ini.grid[ini.lik == min(ini.lik),
                                       ])        
  if(.temp.list$minimisation.function == "optim") 
    phitausq.rel.res <- optim(ini, proflik.aux10, method="L-BFGS-B",
                              lower = c(.temp.list$lower.phi,
                                0), upper=c(.temp.list$upper.phi, 100), ...)$value
  else
    phitausq.rel.res <- nlmP(proflik.aux10, ini,
                             lower = c(.temp.list$lower.phi,
                               0), upper=c(.temp.list$upper.phi, 100), ...)$minimum
  .temp.list$log.jacobian <<- NULL
  .temp.list$sigmasq <<- NULL
  .temp.list$z <<- .temp.list$data
  return( - phitausq.rel.res)
}

"proflik.aux28" <-
  function(sigmasqphi, ...)
{
  ## This function computes the value of the 2-D profile likelihood for the random field scale (variance) parameter \sigma^2  and the correlation function parameter \phi when nugget effect is not included in the model.
  ## It requires the minimisation of the function wrt \lambda for each value of (\sigma^2, \phi)
  ## This is an auxiliary function called by likfit.proflik
  ##
  ini.seq <- seq(-1.5, 1.5, l=7)
  .temp.list$sigmasqphi <<- as.vector(sigmasqphi)
  lambda.lik <- apply(as.matrix(ini.seq), 1, proflik.aux4)
  ini <- ini.seq[lambda.lik == max(lambda.lik)]
  if(.temp.list$minimisation.function == "optim") 
    lambda.res <- optim(ini, proflik.aux4, method="L-BFGS-B", lower = -2.5, upper = 2.5, ...)$value
  else
    lambda.res <- nlmP(proflik.aux4, ini, lower = -2.5, upper = 2.5, ...)$minimum
  .temp.list$sigmasqphi <<- NULL
  return( - lambda.res)
}

"proflik.aux30" <-
  function(phitausq.rel, ...)
{
  ## This function computes the value of the profile likelihood for the correlation parameter \phi when nugget effect is not included in the model.
  ## It requires the minimisation of the function wrt \phi for each value of \lambda, if this transformation parameter is included in the model
  ## This is an auxiliary function called by likfit.proflik
  ##
  if(.temp.list$fix.lambda == TRUE)
    proflik <- proflik.aux21(phitausq.rel = phitausq.rel)
  else {
    .temp.list$phitausq.rel <<- phitausq.rel
    if(.temp.list$minimisation.function == "optim") 
      proflik <-  - (optim(.temp.list$lambda, proflik.aux21.1, method="L-BFGS-B", lower =
                          -2, upper = 2, ...)$value)
    else
      proflik <-  - (nlmP(proflik.aux21.1, .temp.list$lambda, lower =
                          -2, upper = 2, ...)$minimum)
    
    .temp.list$phitausq.rel <<- NULL
  }
  return(proflik)
}
"proflik.aux3" <-
function(phi.lambda, ...)
{
  ## This function computer the negative of the likelihood function for the correlation function scale parameter \phi (and the transformation parameter \lambda) only for models with fixed nugget effect (i.e., when it is not a parameter to be estimated) 
  ## This function is used when computing the profile likelihood for \sigma^2
  ## This is an auxiliary function called by proflik.aux2
  ##  phi <- pmax(phi, .temp.list$lower.phi)
  if(length(phi.lambda) == 2) lambda <- phi.lambda[2] else lambda <- 1
  sigmasq <- .temp.list$sigmasq
  phi <- phi.lambda[1]
  n <- .temp.list$n
  main <- proflik.main(tausq=0, sigmasq=1, phi=phi, lambda = lambda)
  if(.temp.list$method == "ML") {
    neglik <- (n/2) * log(2 * pi) + main$log.det.to.half +
      (n/2) * log(sigmasq) + (0.5/sigmasq) * main$ssresmat - 
        main$log.jacobian
  }
  if(.temp.list$method == "RML") {
    eigentrem <- eigen(main$ixix, symmetric = TRUE, only.values = TRUE)
    neglik <- ((n - beta.size)/2) * log(2 * pi) + main$
    log.det.to.half + ((n - beta.size)/2) * log(sigmasq) +
      (0.5/sigmasq) * main$ssresmat - 0.5 * sum(log(eigentrem$
                                               values)) - main$log.jacobian
  }
  return(as.vector(round(neglik, dig=8)))
}

"proflik.aux31" <-
  function(philambda, ...)
{
  ## This function computes the value of the 2-D profile likelihood for range \phi and the transformation parameter \lambda.
  ## It requires the minimisation of the function wrt \tau^2_{rel} and for each value of (\phi,\lambda).
  ## This is an auxiliary function called by proflik.
  philambda <- as.vector(philambda)
  .temp.list$phi <<- philambda[1]
  .temp.list$lambda <- philambda[2]
  if(.temp.list$minimisation.function == "optim") 
    tausq.rel.res <- optim(.temp.list$tausq.rel.est, proflik.aux8, method="L-BFGS-B", lower = 
                          0, upper=100, ...)$value
  else
    tausq.rel.res <- nlmP(proflik.aux8, .temp.list$tausq.rel.est, lower = 
                          0, upper=100, ...)$minimum
  .temp.list$phi <<- NULL
  return( - tausq.rel.res)
}

"proflik.aux32" <-
  function(tausqlambda, ...)
{
  ## This function computes the value of the 2-D profile likelihood for nugget \tau^2 and the transformation parameter \lambda.
                                        # It requires the minimisation of the function wrt \phi and \sigma^2 and for each value of (\tau^2,\lambda).
                                        # This is an auxiliary function called by proflik.
  tausqlambda <- as.vector(tausqlambda)
  .temp.list$nugget <<- tausqlambda[1]
  lambda <- tausqlambda[2]
  if(lambda == 1) {
    .temp.list$log.jacobian <<- 0
  }
  else {
    if(any(.temp.list$z^(lambda - 1) <= 0))
      .temp.list$log.jacobian <<- log(prod(.temp.list$z^(lambda -
                                                         1)))
    else .temp.list$log.jacobian <<- sum(log(.temp.list$z^(lambda -
                                                           1)))
    if(lambda == 0)
      .temp.list$z <<- log(.temp.list$z)
    else .temp.list$z <<- ((.temp.list$z^lambda) - 1)/lambda
  }
  if(.temp.list$minimisation.function == "optim") 
    sigmasqphi.res <- optim(c(.temp.list$sigmasq.est, .temp.list$phi.est), proflik.aux12, method="L-BFGS-B",
                           lower = c(.temp.list$lower.sigmasq, .temp.list$
                             lower.phi), upper=c(+Inf, .temp.list$upper.phi), ...)$value
  else
    sigmasqphi.res <- nlmP(proflik.aux12, c(.temp.list$sigmasq.est, .temp.list$phi.est),
                           lower = c(.temp.list$lower.sigmasq, .temp.list$
                             lower.phi), upper=c(+Inf, .temp.list$upper.phi), ...)$minimum
  .temp.list$log.jacobian <<- NULL
  .temp.list$z <<- .temp.list$data
  .temp.list$nugget <<- NULL
  return( - sigmasqphi.res)
}

"proflik.aux33" <-
  function(tausq.rellambda, ...)
{
  ## This function computes the value of the 2-D profile likelihood for nugget \tau^2 and the transformation parameter \lambda.
  ## It requires the minimisation of the function wrt \phi for each value of (\tau^2,\lambda).
  ## This is an auxiliary function called by proflik.
  tausq.rellambda <- as.vector(tausq.rellambda)
  .temp.list$nugget.rel <<- tausq.rellambda[1]
  lambda <- tausq.rellambda[2]
  if(lambda == 1) {
    .temp.list$log.jacobian <<- 0
  }
  else {
    if(any(.temp.list$z^(lambda - 1) <= 0))
      .temp.list$log.jacobian <<- log(prod(.temp.list$z^(lambda -
                                                         1)))
    else .temp.list$log.jacobian <<- sum(log(.temp.list$z^(lambda -
                                                           1)))
    if(lambda == 0)
      .temp.list$z <<- log(.temp.list$z)
    else .temp.list$z <<- ((.temp.list$z^lambda) - 1)/lambda
  }
  if(.temp.list$minimisation.function == "optim") 
    phi.res <- optim(.temp.list$phi.est, proflik.aux6, method="L-BFGS-B", lower = .temp.list$
                    lower.phi, upper=.temp.list$upper.phi, ...)$value
  else
    phi.res <- nlmP(proflik.aux6, .temp.list$phi.est, lower = .temp.list$
                    lower.phi, upper=.temp.list$upper.phi, ...)$minimum
  .temp.list$log.jacobian <<- NULL
  .temp.list$nugget.rel <<- NULL
  .temp.list$z <<- .temp.list$data
  return( - phi.res)
}
"proflik.aux4" <-
  function(lambda, ...)
{
  ## This function computer the values of the profile likelihood function for the parameters \phi  and \sigma^2 for models with nugget effect = 0, including the tranformation parameter \lambda
  ## This is an auxiliary function called by proflik.aux28
  ##
  sigmasqphi <- as.vector(.temp.list$sigmasqphi)
  sigmasq <- sigmasqphi[1]
  phi <- sigmasqphi[2]
  n <- .temp.list$n
  if(lambda > 0.999 & lambda < 1.001)
    lambda <- 1
  main <- proflik.main(tausq=.temp.list$nugget, sigmasq = sigmasq, phi=phi, lambda = lambda)
  if(.temp.list$method == "ML") {
    neglik <- ((n/2) * log(2 * pi) +
               main$log.det.to.half +
               0.5 * main$ssresmat - 
               main$log.jacobian)
  }
  if(.temp.list$method == "RML") {
    xx.eigen <- eigen(crossprod(.temp.list$xmat), symmetric = TRUE,
                      only.values = TRUE)
    neglik <- (((n - beta.size)/2) * log(2 * pi) -
               0.5 * sum(log(xx.eigen$values)) +
               main$log.det.to.half +
               (0.5) * main$ssresmat +
               choldet +
               main$log.jacobian)
  }
  return(as.vector(round(neglik, dig=8)))
}
"proflik.aux5" <-
  function(tausq.rel, ...)
{
  ## This function computes the value of the profile likelihood for the parameter \tau^2_{rel}.
  ## It requires the minimisation of the function wrt \phi and \lambda (if the case) for each value of \tau^2_{rel}.
  ## This is an auxiliary function called by proflik.
  .temp.list$nugget.rel <<- as.vector(tausq.rel)
  if(.temp.list$fix.lambda == TRUE) {
    if(.temp.list$minimisation.function == "optim") 
      phi.res <- optim(.temp.list$phi.est, proflik.aux6, method="L-BFGS-B", lower = 
                       .temp.list$lower.phi, upper=.temp.list$upper.phi, ...)$value
    else
      phi.res <- nlmP(proflik.aux6, .temp.list$phi.est, lower = 
                      .temp.list$lower.phi, upper=.temp.list$upper.phi, ...)$minimum
  }
  else {
    if(.temp.list$minimisation.function == "optim") 
      phi.res <- optim(c(.temp.list$phi.est, .temp.list$lambda), proflik.aux6, method="L-BFGS-B", 
                       lower = c(.temp.list$lower.phi, -2),
                       upper = c(.temp.list$upper.phi, 2), ...)$value
    else
      phi.res <- nlmP(proflik.aux6, c(.temp.list$phi.est, .temp.list$lambda), 
                      lower = c(.temp.list$lower.phi, -2),
                      upper = c(.temp.list$upper.phi, 2), ...)$minimum
  }
  .temp.list$nugget.rel <<- NULL
  return( - phi.res)
}

"proflik.aux6" <-
function(phi.lambda, ...)
{
  ## This function computes the value of the profile likelihood for the relative nugget parameter \tau^2_{rel}, minimizing the likelihood wrt correlation function scale parameter \phi (range) and the transformation parameter \lambda.
  if(length(phi.lambda) == 2) lambda <- phi.lambda[2] else lambda <- 1
  phi.lambda <- as.vector(phi.lambda)
  phi <- phi.lambda[1]
  n <- .temp.list$n
  main <- proflik.main(tausq=.temp.list$nugget.rel, sigmasq=1, phi=phi, lambda = lambda)
  if(.temp.list$method == "ML") {
    neglik <- (n/2) * log(2 * pi) +
      main$log.det.to.half +
        (n/2) * log(main$ssresmat/n) +
          (n/2) -
            main$log.jacobian
  }
  if(.temp.list$method == "RML") {
    eigentrem <- eigen(main$ixix, symmetric = TRUE, only.values = TRUE)
    neglik <- ((n - beta.size)/2) * log(2 * pi) +
      main$log.det.to.half +
        ((n - beta.size)/2) * log(main$ssresmat/n) +
          (n/2) -
            0.5 * sum(log(eigentrem$values)) -
              main$log.jacobian
  }
  return(as.vector(round(neglik, dig=8)))
}
"proflik.aux7" <-
  function(phi, ...)
{
  ## This function computes the value of the profile likelihood for the parameter \phi when the nugget \tau^2 is included in the model
  ## It requires the minimisation of the function wrt relative \tau^2_{rel} for each value of \phi
  ## This is an auxiliary function called by proflik.
  .temp.list$phi <<- as.vector(phi)
  if(.temp.list$fix.lambda == TRUE) {
    .temp.list$lambda <<- 1
    if(.temp.list$minimisation.function == "optim") 
      tausq.rel.res <- optim(.temp.list$tausq.rel.est, proflik.aux8, method="L-BFGS-B", 
                             lower = 0, upper=100, ...)$value
    else
      tausq.rel.res <- nlmP(proflik.aux8,.temp.list$tausq.rel.est, 
                            lower = 0, upper=100, ...)$minimum
    .temp.list$lambda <<- NULL
  }
  else {
    if(.temp.list$minimisation.function == "optim") 
      tausq.rel.res <- optim(c(.temp.list$tausq.rel.est, .temp.list$lambda), proflik.aux8, method="L-BFGS-B", lower = c(0, -2), upper = c(100, 2), ...)$value
    else
      tausq.rel.res <- nlmP(proflik.aux8, c(.temp.list$tausq.rel.est, .temp.list$lambda), lower = c(0, -2), upper = c(100, 2), ...)$minimum
  }
  .temp.list$phi <<- NULL
  return( - tausq.rel.res)
}
"proflik.aux8" <-
  function(tausq.rel.lambda, ...)
{
  ## This function computes the value of the profile likelihood for the correlation function scale parameter \phi (and lambda), minimizing the likelihood wrt relative nugget parameter \tau^2_{rel}
  if(length(tausq.rel.lambda) == 2)
    lambda <- tausq.rel.lambda[2]
  else lambda <- .temp.list$lambda
  n <- .temp.list$n
  phi <- .temp.list$phi
  tausq.rel <- tausq.rel.lambda[1]
  main <- proflik.main(tausq=tausq.rel, sigmasq=1, phi=phi, lambda = lambda)
  if(.temp.list$method == "ML") {
    neglik <- (n/2) * log(2 * pi) + main$log.det.to.half + (
                                                              n/2) * log(main$ssresmat/n) + (n/2) - main$log.jacobian
  }
  if(.temp.list$method == "RML") {
    eigentrem <- eigen(main$ixix, symmetric = TRUE, only.values = TRUE)
    neglik <- (((n - beta.size)/2) * log(2 * pi) + main$
               log.det.to.half + ((n - beta.size)/2) * log(main$ssresmat/
                                                           n) + (n/2) - 0.5 * sum(log(eigentrem$values))) - 
                                                             main$log.jacobian
  }
  return(as.vector(round(neglik, dig=8)))
}
"proflik.aux9" <-
  function(sigmasq, ...)
{
  ## This function computes the value of the profile likelihood for the parameter \sigma^2 when \tau^2 is included in the model
  ## It requires the minimisation of the function wrt \phi and \tau^2 for each value of \sigma^2
  ## This is an auxiliary function called by likfit.proflik
  .temp.list$sigmasq <<- as.vector(sigmasq)
  if(.temp.list$fix.lambda == TRUE) {
    ini.lik <- round(100000000. * apply(.temp.list$ini.grid, 1,
                                        proflik.aux10))
    ini <- as.vector(.temp.list$ini.grid[ini.lik == min(ini.lik),
                                         ])
    if(.temp.list$minimisation.function == "optim") 
      phitausq.rel.res <- optim(ini, proflik.aux10, method="L-BFGS-B",
                               lower = c(.temp.list$
                                 lower.phi, 0),
                               upper=c(.temp.list$upper.phi, 100), ...)$value
    else
      phitausq.rel.res <- nlmP(proflik.aux10, ini,
                               lower = c(.temp.list$
                                 lower.phi, 0),
                               upper=c(.temp.list$upper.phi, 100), ...)$minimum
    
  }
  else {
    ini.lik <- apply(.temp.list$ini.grid, 1, proflik.aux10)
    ini <- as.vector(.temp.list$ini.grid[ini.lik == min(ini.lik),])
    if(ini[2] == 0) ini[2] <- 0.01
    if(.temp.list$minimisation.function == "optim") 
      phitausq.rel.res <- optim(ini, proflik.aux10, method="L-BFGS-B", 
                               lower = c(.temp.list$lower.phi,
                                 0,-2),
                               upper = c(.temp.list$upper.phi,
                                 100, 2), ...)$value
    else
      phitausq.rel.res <- nlmP(proflik.aux10,ini, 
                               lower = c(.temp.list$lower.phi,
                                 0,-2),
                               upper = c(.temp.list$upper.phi,
                                 100, 2), ...)$minimum
  }
  .temp.list$sigmasq <<- NULL
  return( - phitausq.rel.res)
}

"plot.proflik" <-
  function(obj.proflik, pages = c("user", "one", "two"),
           uni.only = FALSE, bi.only = FALSE,
           type.bi = c("contour", "persp"),
           conf.int = c(0.90000000000000002,0.94999999999999996),
           yaxis.lims = c("conf.int", "as.computed"),
           by.col = TRUE, log.scale = FALSE, use.splines = TRUE,
           par.mar.persp = c(0, 0, 0, 0),
           ...)
                                        #
{
                                        #  par.mfrow.ori <- par()$mfrow
                                        #  par.mfcol.ori <- par()$mfcol
  par.ori <- par(no.readonly = TRUE)
  on.exit(par(par.ori))
  pages <- match.arg(pages)
  if(all(is.character(yaxis.lims)))
    yaxis.lims <- match.arg(yaxis.lims)
  type.bi <- match.arg(type.bi)
  n.uni <- obj.proflik$n.uni
  n.bi <- obj.proflik$n.bi
  if(n.bi == 0)
    uni.only <- T
  if((uni.only == FALSE) & (bi.only == FALSE))
    np <- n.uni + n.bi
  if((uni.only == TRUE) & (bi.only == FALSE))
    np <- n.uni
  if((uni.only == FALSE) & (bi.only == TRUE))
    np <- n.bi
  if(n.uni==0 & np > 0) bi.only <- T
  if(n.bi==0 & np > 0) uni.only <- T
  if(pages == "one") {
    if(np >= 1 & np < 4)
      par(mfrow = c(np, 1))
    if(np >= 4) {
      if(by.col == TRUE)
        par(mfcol = c(ceiling(np/2), 2))
      else par(mfrow = c(ceiling(np/2), 2))
    }
  }
  if(pages == "two") {
    if(n.uni > 1 & n.uni < 4)
      par(mfrow = c(n.uni, 1))
    if(n.uni >= 4)
      par(mfrow = c(ceiling(n.uni/2), 2))
  }     
  if(bi.only == FALSE) {
    for(i in 1:n.uni) {
      if(obj.proflik$method == "ML")
        ylabm <- "profile log-likelihood"
      else ylabm <- "profile log-(restricted) likelihood"
      if(all(conf.int) != FALSE) {
        if(!is.numeric(conf.int) | any(conf.int > 1))
          stop("argument conf.int must be numerical (scalar or vector) with values between 0 and 1")
        conf.int.drop <- obj.proflik[[i]][[3]][2] - 0.5 * qchisq(conf.int,1)
      }
      if(all(is.character(yaxis.lims))){
        if(yaxis.lims == "conf.int")
          lik.lims <- c(min(conf.int.drop), 
                        obj.proflik[[i]][[3]][2])
        else lik.lims <- c(min(obj.proflik[[i]][[2]]),
                           obj.proflik[[i]][[3]][2])
      }
      else
        lik.lims <- yaxis.lims
      if(log.scale == TRUE) {
        if(use.splines){
          nxpoints <- 5*length(obj.proflik[[i]][[1]])
          plot(spline(x = log(obj.proflik[[i]][[1]]), 
                      y = obj.proflik[[i]][[2]],
                      n = nxpoints,
                      method="natural"), type = "l",
               xlab = paste("log-",
                 plot.proflik.aux1(names(obj.proflik[[i]])[1])),
               ylab = ylabm, ylim = lik.lims)
        }
        else{
          plot(log(obj.proflik[[i]][[1]]), 
               obj.proflik[[i]][[2]],
               type = "l",
               xlab = paste("log-",
                 plot.proflik.aux1(names(obj.proflik[[i]])[1])),
               ylab = ylabm, ylim = lik.lims)
        }
        lines(log(c(obj.proflik[[i]][[3]][1], obj.proflik[[i]][[3]][1])),
              c(min(lik.lims), obj.proflik[[i]][[3]][2]), lty = 2)
      }
      else {
        if(use.splines){
          nxpoints <- 5*length(obj.proflik[[i]][[1]])
          plot(spline(x = obj.proflik[[i]][[1]],
                      y = obj.proflik[[i]][[2]],
                      n = nxpoints,
                      method="natural"),
               type = "l",
               xlab = plot.proflik.aux1(names(obj.proflik[[i]])[1]),
               ylab = ylabm, ylim = lik.lims)
        }
        else{
          plot(spline(x = obj.proflik[[i]][[1]],
                      y = obj.proflik[[i]][[2]],
                      method="natural"), type = "l", xlab = 
               plot.proflik.aux1(names(obj.proflik[[i]])[1]),
               ylab = ylabm, ylim = lik.lims)
        }
        lines(c(obj.proflik[[i]][[3]][1], 
                obj.proflik[[i]][[3]][1]),
              c(min(lik.lims), obj.proflik[[
                                            i]][[3]][2]), lty = 2)
      }
      abline(h = conf.int.drop, lty = 3)
    }
  }
  if(uni.only == FALSE) {
    if(pages == "two") {
      if(n.bi >= 1 & n.bi < 4)
        par(mfrow = c(n.bi, 1))
      if(n.bi >= 4)
        par(mfrow = c(ceiling(n.bi/2), 2))
    }
    for(i in 1:n.bi) {
      if(type.bi == "contour") {
        if(log.scale == TRUE) {
          contour(log(obj.proflik[[(n.uni + i)]][[1]]),
                  log(obj.proflik[[(n.uni + i)]][[2]]),
                  matrix(obj.proflik[[(n.uni + i)]][[3]],
                         ncol = length(obj.proflik[[(n.uni +i)]][[2]])),
                  xlab = paste("log-", plot.proflik.aux1(names(obj.proflik[[(n.uni + i)]][1]))),
                  ylab = paste("log-", plot.proflik.aux1(names(obj.proflik[[(n.uni + i)]][2]))),
                    ...)
          points(log(t(obj.proflik[[(n.uni + i)]][[4]][1:2])))
        }
        else {
          contour(obj.proflik[[(n.uni + i)]][[1]],
                  obj.proflik[[(n.uni + i)]][[2]],
                  matrix(obj.proflik[[(n.uni + i)]][[3]],
                         ncol = length(obj.proflik[[(n.uni + i)]][[2]])),
                  xlab = plot.proflik.aux1(names(obj.proflik[[(n.uni + i)]][1])),
                  ylab = plot.proflik.aux1(names(obj.proflik[[(n.uni + i)]][2])),
                  ...)
          points(t(obj.proflik[[(n.uni + i)]][[4]][1:2]))
        }
      }
      if(type.bi == "persp") {
        cat("For better visualisation arguments for the funtion `persp` can be passed.\nSome relevant argments are: theta, phi, r, d, among others.\n Type help(persp) for a description of the options\n")
        if(obj.proflik$method == "ML")
          zlabm <- 
            "profile log-likelihood"
        else zlabm <- "profile log-(restricted) likelihood"
        zlimm <- range(obj.proflik[[(n.uni +
                                     i)]][[3]])
        zlimm[1] <- 1.01 * zlimm[1]
        minlik <- min(obj.proflik[[(n.uni + i)]][[3]])
        if(log.scale == TRUE) {
          persp(log(obj.proflik[[(n.uni + i)]][[1]]),
                log(obj.proflik[[(n.uni + i)]][[2]]),
                matrix(obj.proflik[[(n.uni + i)]][[3]],
                       ncol = length(obj.proflik[[(n.uni + i)]][[2]])),
                xlab = plot.proflik.aux1(paste("log-", names(obj.proflik[[(n.uni + i)]][1]))),
                ylab = paste("log-", plot.proflik.aux1(names(obj.proflik[[(n.uni + i)]][2]))),
                zlab = zlabm, box = T, ...)
                                        #          pp1 <- perspp(x = log(c(obj.proflik[[(n.uni + i)]][[4]][1],
                                        #                          min( obj.proflik[[(n.uni + i)]][[1]]))[c(1, 1, 1, 2)]),
                                        #                        y = log(c(obj.proflik[[(n.uni + i)]][[4]][2],
                                        #                          min(obj.proflik[[(n.uni + i)]][[2]]))[c(1, 1, 2, 1)]),
                                        #                        z = c(minlik, obj.proflik[[(n.uni + i)]][[4]][3])[c(1, 2, 1, 1)], pp)
                                        #          segments(log(pp1$x[1]), log(pp1$y[1]), log(pp1$x[2]), log(pp1$y[2]),
                                        #                   lwd = 2)
        }
        else {
          persp(x = obj.proflik[[(n.uni + i)]][[1]],
                y = obj.proflik[[(n.uni + i)]][[2]],
                z = matrix(obj.proflik[[(n.uni + i)]][[3]],
                  ncol = length(obj.proflik[[(n.uni + i)]][[2]])),
                xlab = plot.proflik.aux1(names(obj.proflik[[(n.uni + i)]][1])),
                ylab = plot.proflik.aux1(names(obj.proflik[[(n.uni + i)]][2])),
                zlab = zlabm, box = TRUE, ...)
                                        #          pp1 <- perspp(x = c(obj.proflik[[(n.uni + i)]][[4]][1],
                                        #                          min(obj.proflik[[(n.uni + i)]][[1]]))[c(1, 1, 1, 2)],
                                        #                        y = c(obj.proflik[[(n.uni + i)]][[4]][2],
                                        #                          min(obj.proflik[[(n.uni + i)]][[2]]))[c(1, 1,2, 1)],
                                        #                        z = c(minlik, obj.proflik[[(n.uni + i)]][[4]][3])[c(1, 2, 1, 1)], pp)
                                        #          segments(pp1$x[1], pp1$y[1], pp1$x[2], pp1$y[2], lwd = 2)
        }
      }
    }
  }
  return(invisible())
}

"plot.proflik.aux1" <-
  function(parameter.name)
{
  switch(parameter.name,
         range = expression(phi),
         sill = expression(sigma^2),
         lambda = expression(lambda),
         nugget = expression(tau^2),
         nugget.rel = expression(tau^2[rel]))
}

"proflik.main" <-
  function(tausq, sigmasq, phi, lambda)
{
  z <- .temp.list$z
  n <- .temp.list$n
  if(lambda == 1) {
    log.jacobian <- 0
  }
  else {
    if(any(z <= 0))
      stop("Transformation option not allowed when there are zeros or negative data"
           )
    if(any(z^(lambda - 1) <= 0))
      log.jacobian <- log(prod(z^(lambda - 1)))
    else log.jacobian <- sum(log(z^(lambda - 1)))
    if(lambda == 0)
      z <- log(z)
    else z <- ((z^lambda) - 1)/lambda
  }
  beta.size <- .temp.list$beta.size
  kappa <- .temp.list$kappa
  covinf <- varcov.spatial(dists.lowertri = .temp.list$dists.lowertri,
                           cov.model = .temp.list$cov.model, kappa = kappa,
                           nugget = tausq, cov.pars = c(sigmasq, phi),
                           inv = TRUE, det = TRUE, func.inv = "eigen",
                           only.inv.lower.diag = TRUE)  
  xix <- as.double(rep(0, beta.size*beta.size))
  xix <- .C("bilinearform_XAY",
            as.double(covinf$lower.inverse),
            as.double(covinf$diag.inverse),
            as.double(as.vector(.temp.list$xmat)),
            as.double(as.vector(.temp.list$xmat)),
            as.integer(beta.size),
            as.integer(beta.size),
            as.integer(n),
            res = xix)$res
  attr(xix, "dim") <- c(beta.size, beta.size)
  if(length(as.vector(xix)) == 1) {
    ixix <- 1/xix
    choldet <- 0.5 * log(xix)
  }
  else {
    chol.xix <- chol(xix)
    ixix <- chol2inv(chol.xix)
    choldet <- sum(log(diag(chol.xix)))
  }
  xiy <- as.double(rep(0, beta.size))
  xiy <- .C("bilinearform_XAY",
            as.double(covinf$lower.inverse),
            as.double(covinf$diag.inverse),
            as.double(as.vector(.temp.list$xmat)),
            as.double(as.vector(z)),
            as.integer(beta.size),
            as.integer(1),
            as.integer(n),
            res = xiy)$res
  beta.hat <- as.vector(ixix %*% xiy)
  yiy <- as.double(0.0)
  yiy <- .C("bilinearform_XAY",
            as.double(covinf$lower.inverse),
            as.double(covinf$diag.inverse),
            as.double(as.vector(z)),
            as.double(as.vector(z)),
            as.integer(1),
            as.integer(1),
            as.integer(n),
            res = yiy)$res
  ssresmat <- as.vector(yiy - crossprod(beta.hat,xiy))
  return(list(log.det.to.half = covinf$log.det.to.half,
              ssresmat = ssresmat,
              ixix = ixix, log.jacobian = log.jacobian))
}
