"ksline" <-
  function (geodata, coords=geodata$coords, data=geodata$data, locations,
            cov.model = "matern",
            cov.pars = stop("covariance parameters (sigmasq and phi) needed"), 
            kappa = 0.5, nugget = 0, micro.scale = 0,
            lambda = 1, m0 = "ok", nwin = "full", 
            n.samples.backtransform = 500, 
            trend = 1, d = 2, ktedata = NULL, ktelocations = NULL,
            aniso.pars = NULL,  signal = FALSE,  dist.epsilon = 1e-10,
            messages.screen = TRUE) 
{
  require(mva)
  call.fc <- match.call()
  cov.model <- match.arg(cov.model,
                         choices = c("matern", "exponential", "gaussian",
                           "spherical", "circular", "cubic", "wave", "power",
                           "powered.exponential", "cauchy", "gneiting",
                           "gneiting.matern", "pure.nugget"))
  if(lambda != 1) {
    if(messages.screen)
      cat("ksline: Data transformation (Box-Cox) performed.\n")
    if(lambda == 0)
      data <- log(data)
    else data <- ((data^lambda) - 1)/lambda
  }
  coords <- as.matrix(coords)
  locations <- as.matrix(locations)
  dimnames(coords) <- list(NULL, NULL)
  dimnames(locations) <- list(NULL, NULL)
  if (!is.null(ktedata) & !is.null(ktelocations) & m0 != "kte"){
    cat("ksline: external variable (covariate) provided. Kriging ste to KTE\n")
    m0 <- "kte"
  }
  ##
  ## anisotropy correction
  ##
  if(!is.null(aniso.pars)) {
    if(length(aniso.pars) != 2 | !is.numeric(aniso.pars))
      stop("anisotropy parameters must be provided as a numeric vector with two elements: the rotation angle (in radians) and the anisotropy ratio (a number greater than 1)")
    if(messages.screen)
      cat("ksline: anisotropy correction performed\n")
    coords.c <- coords.aniso(coords = coords, aniso.pars = aniso.pars)
    locations.c <- coords.aniso(coords = locations, aniso.pars = aniso.pars)
  }
  else {
    coords.c <- coords
    locations.c <- locations
  }
  ## 2. Preparing KTE matrices #####
  ##  
  if (m0 == "kte") {
    ktedata <- as.matrix(ktedata)
    ktelocations <- as.matrix(ktelocations)
    dimnames(ktedata) <- list(NULL, NULL)
    dimnames(ktelocations) <- list(NULL, NULL)
  }
  n <- length(data)
  ni <- length(locations[, 1])
  tausq <- nugget
  sigmasq <- cov.pars[1]
  phi <- cov.pars[2]
  if (nwin == "full") {
    est <- rep(0, ni)
    dif <- rep(0, ni)
    kvar <- rep(0, ni)
    sumw <- rep(0, ni)
    wofmean <- rep(0, ni)
    iv <- varcov.spatial(coords = coords.c, cov.model = cov.model, 
                         kappa = kappa, nugget = nugget, cov.pars = cov.pars, 
                         inv = TRUE, det = FALSE, func.inv = "cholesky")$inverse
    av <- mean(data)
    sd <- sqrt(var(data))
    one <- rep(1, n)
    tone <- t(one)
    toneiv <- crossprod(one, iv)
    den <- solve(toneiv %*% one)
    ml <- den %*% toneiv %*% data
    kmsd <- sqrt(den)
    means <- c(average = av, stdev = sd, kmean = ml, kmsd = kmsd)
    if (m0 != "kt") {
      mktlocations <- "Constant trend"
      beta <- ml
    }
    else {
      mktlocations <- rep(0, ni)
      if (m0 == "kt" & trend == 1) {
        if (d == 1) {
          xmat <- cbind(rep(1, n), coords[, 2])
          xmati <- cbind(rep(1, ni), locations[, 2])
        }
        else {
          xmat <- cbind(rep(1, n), coords[, 1], coords[, 2])
          xmati <- cbind(rep(1, ni), locations[, 1], locations[, 
                                                               2])
        }
        iviv <- solve(crossprod(xmat,iv) %*% xmat)
        txiv <- crossprod(xmat,iv)
        beta <- iviv %*% txiv %*% data
        mkt <- xmat %*% beta
      }
      if (m0 == "kt" & trend == 2) {
        if (d == 1) {
          xmat <- cbind(rep(1, n), coords[, 2], (coords[, 2])^2)
          xmati <- cbind(rep(1, ni), locations[, 2], (locations[, 
                                                                2])^2)
        }
        else {
          xmat <- cbind(rep(1, n), coords[, 1], coords[, 2], 
                        (coords[, 1])^2, (coords[, 2])^2, coords[, 1] * coords[, 
                                                                               2])
          xmati <- cbind(rep(1, ni), locations[, 1], locations[, 
                                                               2], (locations[, 1])^2, (locations[, 2])^2, locations[, 
                                                                                                                     1] * locations[, 2])
        }
        iviv <- solve(crossprod(xmat,iv) %*% xmat)
        txiv <- crossprod(xmat,iv)
        beta <- iviv %*% txiv %*% data
        mkt <- xmat %*% beta
      }
    }
    if (m0 != "kte") 
      mktelocations <- "No external trend"
    else {
      if (m0 == "kte") {
        mktelocations <- rep(0, ni)
        xmat <- cbind(rep(1, n), ktedata)
        xmati <- cbind(rep(1, ni), ktelocations)
        iviv <- solve(crossprod(xmat,iv) %*% xmat)
        txiv <- crossprod(xmat,iv)
        beta <- iviv %*% txiv %*% data
        mkte <- xmat %*% beta
      }
    }
    for (i in 1:ni) {
      if (messages.screen) {
        if (ni < 11) 
          cat(paste("ksline: kriging location: ", i, "out of", 
                    ni, "\n"))
        else {
          if (ni < 101 & (i%%10 == 1)) 
            cat(paste("ksline: kriging location: ", i, "out of", 
                      ni, "\n"))
          if (ni > 100 & i%%100 == 1) 
            cat(paste("ksline: kriging location: ", i, "out of", 
                      ni, "\n"))
          if (i == ni) 
            cat(paste("ksline: kriging location: ", i, "out of", 
                      ni, "\n"))
        }
      }
      coords0 <- cbind((coords.c[, 1] - locations.c[i, 1]), (coords.c[, 2] - locations.c[i, 
                                                                                 2]))
      dm0 <- sqrt(coords0[, 1]^2 + coords0[, 2]^2)
      v0 <- cov.spatial(obj = dm0, cov.model = cov.model, 
                        kappa = kappa, cov.pars = cov.pars)
      v0[dm0 < dist.epsilon] <- micro.scale + sigmasq
      tv0 <- t(v0)
      v0iv <- crossprod(v0, iv)
      v0ivv0 <- v0iv %*% v0
      skw <- crossprod(v0,iv)
      wofmean[i] <- 1 - sum(skw)
      ##
      ## 4.2.1 Simple kriging with known mean
      ##
      if (is.numeric(m0) == TRUE) {
        dif[i] <- skw %*% (data - m0)
        est[i] <- m0 + dif[i]
        if (signal == TRUE) 
          kvar[i] <- sigmasq - v0ivv0
        else kvar[i] <- tausq + sigmasq - v0ivv0
        sumw[i] <- sum(skw)
      }
      ##
      ## 4.2.2 Simple kriging with data average mean
      ##
      if (m0 == "av") {
        dif[i] <- skw %*% (data - av)
        est[i] <- av + dif[i]
        if (signal == TRUE) 
          kvar[i] <- sigmasq - v0ivv0
        else kvar[i] <- tausq + sigmasq - v0ivv0
        sumw[i] <- sum(((tone/n) + skw - ((skw %*% one %*% 
                                           tone)/n)))
      }
      ##
      ## 4.2.3 Ordinary kriging (or SK with G.L.S. mean)
      ##
      if(m0 == "ok") {
        dif[i] <- skw %*% (data - ml)
        est[i] <- ml + dif[i]
        redu <- as.vector(1 - toneiv %*% v0)
        if(signal == TRUE)
          kvar[i] <- sigmasq - v0ivv0 + (redu %*%
                                         den %*% redu)
        else kvar[i] <- tausq + sigmasq - v0ivv0 + (
                                                    redu %*% den %*% redu)
        sumw[i] <- sum((den %*% one + tv0 - v0iv %*% 
                        one %*% den %*% tone) %*% iv)
      }
      ##
      ## 4.2.4 Universal Kriging (or Kriging with trend model) 
      ##
      if(m0 == "kt") {
        dif[i] <- skw %*% (data - mkt)
        est[i] <- xmati[i,  ] %*% beta + dif[i]
        redu <- as.vector(xmati[i,  ]) - as.vector(
                                                   txiv %*% v0)
        if(signal == TRUE)
          kvar[i] <- sigmasq - v0ivv0 + (redu %*%
                                         iviv %*% redu)
        else kvar[i] <- tausq + sigmasq - v0ivv0 + (
                                                    redu %*% iviv %*% redu)
        sumw[i] <- sum(skw + xmati[i,  ] %*% iviv %*%
                       txiv - skw %*% xmat %*% iviv %*% txiv)
        mktlocations[i] <- xmati[i,  ] %*% beta
      }
      ##
      ## 4.2.5 Kriging with external trend 
      ##
      if(m0 == "kte") {
        dif[i] <- skw %*% (data - mkte)
        est[i] <- xmati[i,  ] %*% beta + dif[i]
        redu <- as.vector(xmati[i,  ]) - as.vector(
                                                   txiv %*% v0)
        if(signal == TRUE)
          kvar[i] <- sigmasq - v0ivv0 - (redu %*%
                                         iviv %*% redu)
        else kvar[i] <- tausq + sigmasq - v0ivv0 + (
                                                    redu %*% iviv %*% redu)
        sumw[i] <- sum(skw + xmati[i,  ] %*% iviv %*%
                       txiv - skw %*% xmat %*% iviv %*% txiv)
        mktelocations[i] <- xmati[i,  ] %*% beta
      }
      NULL
    }
    message <- "Kriging performed in global neighbourhood"
    if (messages.screen) 
      cat(paste(message,"\n"))
    results <- list(predict = est, krige.var = kvar, dif = dif, summary = means, 
                    ktrend = mktlocations, ktetrend = mktelocations, beta = beta, 
                    wofmean = wofmean)
  }
  else {
    nwin <- min(n, nwin)
    avwin <- rep(0, ni)
    sdwin <- rep(0, ni)
    mlwin <- rep(0, ni)
    kmsdwin <- rep(0, ni)
    estwin <- rep(0, ni)
    difwin <- rep(0, ni)
    kvarwin <- rep(0, ni)
    sumwwin <- rep(0, ni)
    wofmean <- rep(0, ni)
    if (m0 != "kt") 
      mkt <- "Constant position trend"
    else mkt <- rep(0, ni)
    if (m0 != "kte") 
      mkte <- "No external trend"
    else mkte <- rep(0, ni)
    if (m0 != "kt" & m0 != "kte") 
      betawin <- "No polynomial or external trend"
    if (m0 == "kt") {
      if (trend == 1) {
        if (d == 1) 
          xmati <- cbind(rep(1, ni), locations[, 2])
        else xmati <- cbind(rep(1, ni), locations[, 1], locations[, 
                                                                  2])
      }
      if (trend == 2) {
        if (d == 1) 
          xmati <- cbind(rep(1, ni), locations[, 2], locations[, 
                                                               2]^2)
        else xmati <- cbind(rep(1, ni), locations[, 1], locations[, 
                                                                  2], (locations[, 1])^2, (locations[, 2])^2, locations[, 1] * 
                            locations[, 2])
      }
      betawin <- matrix(0, nrow = (ncol(xmati) * ni), ncol = ncol(xmati))
    }
    if (m0 == "kte") {
      xmati <- cbind(rep(1, ni), ktelocations)
      if (is.vector(ktedata) == TRUE) 
        betawin <- matrix(0, nrow = (2 * ni), ncol = 2)
      else betawin <- matrix(0, nrow = ((ncol(ktedata) + 
                                         1) * ni), ncol = (ncol(ktedata) + 1))
    }
    for (i in 1:ni) {
      temp.win <- ksline.aux.1(coords = coords, coords.c = coords.c,
                               data = data, n = n,
                               locations = locations[i,  ],
                               locations.c = locations.c[i,  ],
                               cov.pars = cov.pars, nugget = nugget,
                               cov.model = cov.model, kappa = kappa, m0 = m0,
                               nwin = nwin, trend = trend, d = d, ktedata = 
                               ktedata, ktelocations = ktelocations,
                               micro.scale = micro.scale, 
                               location.number = i, xmati = xmati[i,  ],
                               mkte = NULL, mkt = NULL, betawin = NULL,
                               signal = signal, dist.epsilon = dist.epsilon)
      avwin[i] <- temp.win$avwin
      sdwin[i] <- temp.win$sdwin
      mlwin[i] <- temp.win$mlwin
      kmsdwin[i] <- temp.win$kmsdwin
      estwin[i] <- temp.win$estwin
      difwin[i] <- temp.win$difwin
      kvarwin[i] <- temp.win$kvarwin
      sumwwin[i] <- temp.win$sumwwin
      wofmean[i] <- temp.win$wofmean
      if (m0 == "kt") 
        mkt[i] <- temp.win$mkt
      if (m0 == "kte") 
        mkte[i] <- temp.win$mkte
      if (m0 == "kt" | m0 == "kte") 
        betawin[i, ] <- temp.win$betawin
      if (messages.screen) {
        if (ni < 11) 
          cat(paste("ksline: kriging location: ", i, "out of", 
                    ni, "\n"))
        else {
          if (ni < 101 & (i%%10 == 1)) 
            cat(paste("ksline: kriging location: ", i, "out of", 
                      ni, "\n"))
          if (ni > 100 & i%%100 == 1) 
            cat(paste("ksline: kriging location: ", i, "out of", 
                      ni, "\n"))
          if (i == ni) 
            cat(paste("ksline: kriging location: ", i, "out of", 
                      ni, "\n"))
        }
      }
    }
    message <- "kriging performed in moving neighbourhood"
    if (messages.screen) 
      cat(paste(message,"\n"))
    results <- list(predict = estwin, krige.var = kvarwin, dif = difwin, 
                    avtrend = avwin, sd = sdwin, oktrend = mlwin, oksd = kmsdwin, 
                    ktrend = mkt, ktetrend = mkte, beta = betawin, sumw = sumwwin, 
                    wofmean = wofmean)
  }  
  if(lambda != 1) {
    if(messages.screen)
      cat("Back-transforming the predictions according to the (Box-Cox) parameter lambda\n")
    if(lambda == 0) {
      predict.transf <- results$predict
      results$predict <- exp(predict.transf) - 0.5 * results$krige.var
      results$krige.var <- (exp(2 * predict.transf - results$krige.var)) * (exp(results$krige.var) - 1)
    }
    if(lambda > 0) {
      if(messages.screen)
        cat("Back-transformation done by sampling from the resulting (normal) predictive distribution\n")
      ap.warn <- options()$warn
      options(warn = -1)
      temp.data <- matrix(rnorm(ni * n.samples.backtransform,
                                mean = results$predict, sd = sqrt(results$krige.var)),
                          nrow = ni)
      options(warn = ap.warn)
      temp.data[(results$krige.var == 0),  ] <- results$predict[(results$krige.var == 0)]
      temp.data[temp.data < -1/lambda] <- -1/lambda     
      temp.data <- ((temp.data * lambda) + 1)^(1/lambda)
###      temp.data[is.na(temp.data)] <- Inf
      results$predict <- as.vector(apply(temp.data, 1, mean))
      results$krige.var <- as.vector(apply(temp.data, 1, var))
    }
    if(lambda < 0) {
      cat("Resulting distribution has no mean for lambda < 0 - back transformation not performed\n"
          )
    }
  }
  results$locations <- locations
  results$message <- message
  results$call <- call.fc
  class(results) <- c("kriging")
  return(invisible(results))
}

"ksline.aux.1" <-
  function (coords, coords.c, data, n, locations, locations.c, cov.pars,
            nugget, cov.model, kappa, 
            m0, nwin, trend, d, ktedata, ktelocations, mbased,
            micro.scale, location.number, 
            xmati, mkte, mkt, betawin, signal, dist.epsilon) 
{
  require(mva)
  i <- location.number
  sigmasq <- cov.pars[1]
  phi <- cov.pars[2]
  tausq <- nugget
  coords0 <- cbind((coords.c[, 1] - locations.c[1]), (coords.c[, 2] -
                                                      locations.c[2]))
  dm0 <- sqrt(coords0[, 1]^2 + coords0[, 2]^2)
  coordswin <- coords[order(dm0)[1:nwin],  ]
  coordswin.c <- coords.c[order(dm0)[1:nwin],  ]
  datawin <- data[order(dm0)[1:nwin]]
  ivwin <- varcov.spatial(coords = coordswin.c, cov.model = cov.model,
                          kappa = kappa, nugget = nugget, cov.pars = cov.pars, inv = TRUE,
                          det = FALSE, func.inv = "cholesky", only.decomp = FALSE)$inverse
  avwin <- mean(datawin)
  sdwin <- sqrt(var(datawin))
  onewin <- rep(1, nwin)
  toneivwin <- crossprod(onewin, ivwin)
  denwin <- solve(toneivwin %*% onewin)
  mlwin <- denwin %*% toneivwin %*% datawin
  kmsdwin <- sqrt(denwin)
  coords0win <- cbind((coordswin[, 1] - locations[1]), (coordswin[, 2] -
                                                        locations[2]))
  coords0win.c <- cbind((coordswin.c[, 1] - locations.c[1]), (coordswin.c[
                                                                          , 2] - locations.c[2]))
  dm0win <- sqrt(coords0win.c[, 1]^2 + coords0win.c[, 2]^2)
  v0win <- cov.spatial(obj = dm0win, cov.model = cov.model, kappa = kappa,
                       cov.pars = cov.pars)
  v0win[dm0win < dist.epsilon] <- micro.scale + sigmasq
  skwwin <- crossprod(v0win, ivwin)
  wofmean <- 1 - sum(skwwin)
  if(m0 == "kt" & trend == 1) {
    if(d == 1)
      xmatwin <- cbind(rep(1, nwin), coordswin[, 2])
    else xmatwin <- cbind(rep(1, nwin), coordswin[, 1], coordswin[
                                                                  , 2])
    txivwin <- crossprod(xmatwin, ivwin)
    ivivwin <- solve(txivwin %*% xmatwin)
    betawin <- ivivwin %*% txivwin %*% datawin
    mktwin <- xmatwin %*% betawin
    mkt <- xmati %*% betawin
  }
  if(m0 == "kt" & trend == 2) {
    if(d == 1)
      xmatwin <- cbind(rep(1, nwin), coordswin[, 2], (
                                                      coordswin[, 2])^2)
    else xmatwin <- cbind(rep(1, nwin), coordswin[, 1], (coordswin[
                                                                   , 1])^2, coordswin[, 2], (coordswin[, 2])^
                          2, coordswin[, 1] * coordswin[, 2])
    xmatwin.cent <- xmatwin
    xmatwin.cent[, 2] <- xmatwin.cent[, 2] - mean(xmatwin[, 2])
    xmatwin.cent[, 3] <- xmatwin.cent[, 3] - mean(xmatwin[, 3])
    ivivwin <- solve(crossprod(xmatwin.cent, ivwin) %*% 
                     xmatwin.cent)
    txivwin <- crossprod(xmatwin.cent, ivwin)
    betawin <- ivivwin %*% txivwin %*% datawin
    betawin <- mean(datawin) - crossprod(betawin, c(0, mean(xmatwin[
                                                                    , 2]), mean(xmatwin[, 3])))
    mktwin <- xmatwin %*% betawin
    mkt <- xmati %*% betawin
  }
  if(m0 == "kte") {
    if(is.vector(ktedata))
      ktedatawin <- ktedata[order(dm0)[1:nwin]]
    else ktedatawin <- ktedata[order(dm0)[1:nwin],  ]
    xmatwin <- cbind(rep(1, nwin), ktedatawin)
    ivivwin <- solve(crossprod(xmatwin, ivwin) %*% xmatwin)
    txivwin <- crossprod(xmatwin, ivwin)
    betawin <- ivivwin %*% txivwin %*% datawin
    mktewin <- xmatwin %*% betawin
    mkte <- xmati %*% betawin
  }
  ##
  ##  Simple kriging with known mean
  ##
  if(is.numeric(m0)) {
    difwin <- skwwin %*% (data - m0)
    estwin <- m0win + difwin
    if(signal)
      kvarwin <- sigmasq - crossprod(v0win, ivwin) %*% v0win
    else kvarwin <- tausq + sigmasq - crossprod(v0win, ivwin) %*%
      v0win
    sumwwin <- sum(skwwin)
  }
  ##
  ## 4.2.2 Simple kriging with data average mean
  ##
  if(m0 == "av") {
    difwin <- skwwin %*% (datawin - avwin)
    estwin <- avwin + difwin
    if(signal)
      kvarwin <- sigmasq - crossprod(v0win, ivwin) %*% v0win
    else kvarwin <- tausq + sigmasq - crossprod(v0win, ivwin) %*%
      v0win
    sumwwin <- sum(((t(onewin)/nwin) + skwwin - ((skwwin %*% onewin %*%
                                                  t(onewin))/n)))
  }
  ##
  ## Ordinary kriging (or SK with G.L.S. mean)
  ##
  if(m0 == "ok") {
    difwin <- skwwin %*% (datawin - mlwin)
    estwin <- mlwin + difwin
    redu <- as.vector(1 - toneivwin %*% v0win)
    if(signal)
      kvarwin <- sigmasq - v0win %*% ivwin %*% v0win + (
                                                        redu %*% denwin %*% redu)
    else kvarwin <- tausq + sigmasq - v0win %*% ivwin %*% v0win +
      (redu %*% denwin %*% redu)
    sumwwin <- sum((denwin %*% onewin + t(v0win) - crossprod(v0win,
                                                             ivwin) %*% onewin %*% denwin %*% t(onewin)) %*% ivwin)
  }
  ##
  ## Universal Kriging (or Kriging with trend model) 
  ##
  if(m0 == "kt") {
    difwin <- skwwin %*% (datawin - mktwin)
    estwin <- mkt + difwin
    xmati <- as.vector(xmati)
    redu <- as.vector(xmati) - as.vector(txivwin %*% v0win)
    if(signal)
      kvarwin <- sigmasq - (v0win %*% ivwin %*% v0win) + (redu %*% ivivwin %*% redu)
    else kvarwin <- tausq + sigmasq - (v0win %*% ivwin %*% v0win) + (redu %*% ivivwin %*% redu)
    sumwwin <- sum(skwwin + xmati %*% ivivwin %*% txivwin - skwwin %*%
                   xmatwin %*% ivivwin %*% txivwin)
  }
  ##
  ## Kriging with external trend 
  ##
  if(m0 == "kte") {
    difwin <- skwwin %*% (datawin - mktewin)
    estwin <- mkte + difwin
    xmati <- as.vector(xmati)
    redu <- as.vector(xmati) - as.vector(txivwin %*% v0win)
    if(signal)
      kvarwin <- sigmasq - (v0win %*% ivwin %*% v0win) + (redu %*% ivivwin %*% redu)
    else kvarwin <- tausq + sigmasq - (v0win %*% ivwin %*% v0win) + (redu %*% ivivwin %*% redu)
    sumwwin <- sum(skwwin + xmati %*% ivivwin %*% txivwin - skwwin %*%
                   xmatwin %*% ivivwin %*% txivwin)
  }
  ##
  ##
  ##  
  results <- list(avwin = avwin, sdwin = sdwin, mlwin = mlwin, kmsdwin = 
                  kmsdwin, wofmean = wofmean, betawin = betawin, mkt = mkt, mkte
                  = mkte, difwin = difwin, estwin = estwin, kvarwin = kvarwin,
                  sumwwin = sumwwin)
  return(results)
}

"krige.conv" <-
  function (geodata, coords=geodata$coords, data=geodata$data, locations,
            krige = krige.control(
              type.krige, beta = NULL, trend.d, trend.l,
              cov.model, cov.pars, kappa = 0.5, 
              nugget = 0, micro.scale = 0,
              dist.epsilon = 1e-10,
              aniso.pars = NULL, lambda = 1,
              signal = FALSE,
              n.samples.backtransform = 500, n.sim = 0),
            messages.screen = TRUE)
{
  call.fc <- match.call()
  ##
  ## reading input
  ##
  cov.model <- krige$cov.model
  kappa <- krige$kappa
  lambda <- krige$lambda
  ##
  beta <- krige$beta
  cov.pars <- krige$cov.pars
  nugget <- krige$nugget
  ##
  signal <- krige$signal
  n.sim <- krige$n.sim
  n.samples.backtransform <- krige$n.samples.backtransform
  micro.scale <- krige$micro.scale
  ##
  ## checking input
  ##
  if(micro.scale > nugget)
    stop("krige.conv: krige$micro.scale must be in the interval [0, nugget]")
  if (krige$type.krige != "ok" & krige$type.krige != "OK" & krige$type.krige != "o.k." & krige$type.krige != "O.K." & krige$type.krige != "sk" & krige$type.krige != "SK" & krige$type.krige != "s.k." & krige$type.krige != "S.K.")
    stop("krige.conv: wrong option in the argument type.krige. It should be \"OK\" or \"SK\"(if ordinary or simple kriging is to be performed)")
  if (krige$type.krige == "ok" | krige$type.krige == "OK" | krige$type.krige == "o.k." | krige$type.krige == "O.K.") 
    beta.prior <- "flat"
  if (krige$type.krige == "sk" | krige$type.krige == "SK" | krige$type.krige == "s.k." | krige$type.krige == "S.K."){
    if(is.null(beta) | !is.numeric(beta))
      stop("krige.conv: argument beta must be provided in order to perform simple kriging")
    beta.prior <- "deg"
  }
  ##
  if(is.vector(coords)){
    coords <- cbind(coords, 0)
    warning("krige.conv: vector of coordinates, one spatial dimension assumed")
  }
  coords <- as.matrix(coords)
  if (is.vector(locations)) {
    if (length(locations) == 2) {
      locations <- t(as.matrix(locations))
      if (messages.screen) 
        warning("krige.conv: assuming that there is 1 prediction point")
    }
    else{
      warning("krige.conv: vector of locations: one spatial dimension assumed")
      locations <- as.matrix(cbind(locations, 0))
    }
  }
  else locations <- as.matrix(locations)  
  dimnames(coords) <- list(NULL, NULL)
  dimnames(locations) <- list(NULL, NULL)
  ##
  if(inherits(krige$trend.d, "formula") | inherits(krige$trend.l, "formula")){
    if((inherits(krige$trend.d, "formula") == FALSE) | (inherits(krige$trend.l, "formula") == FALSE))
      stop("krige.conv: krige$trend.d and krige$trend.l must have similar specification")
    if(messages.screen)
      cat("krige.conv: Kriging with external trend to be performed using covariates provided by the user\n")
  }
  else{
    if (krige$trend.d != krige$trend.l){
      stop("krige.conv: krige$trend.l is different from krige$trend.d")
    }
    if(messages.screen){
      if(krige$trend.d == "cte" & krige$type.krige == "sk")
        cat("krige.conv: Simple kriging to be performed with constant mean provided by the user\n"
            )
      if(krige$trend.d == "cte" & krige$type.krige == "ok")
        cat("krige.conv: Ordinary kriging to be performed filtering a constant mean\n")
      if(krige$trend.d == "1st")
        cat("krige.conv: Trend (or universal) kriging to be performed filtering a 1st degree polinomial trend\n")
      if(krige$trend.d == "2nd") 
        cat("krige.conv: Trend (or universal) kriging to be performed filtering a 2nd degree polinomial trend\n")
    }
  }
  trend.d <- trend.spatial(trend=krige$trend.d, coords=coords)
  beta.size <- ncol(trend.d)
  trend.l <- trend.spatial(trend=krige$trend.l, coords=locations)
  ##
  ## Anisotropy correction (should be placed AFTER trend.d/trend.l
  ##
  if(!is.null(krige$aniso.pars)) {
    if(length(krige$aniso.pars) != 2 | !is.numeric(krige$aniso.pars))
      stop("krige.conv: anisotropy parameters must be provided as a numeric vector with two elements: the rotation angle (in radians) and the anisotropy ratio (a number greater than 1)")
    if(messages.screen)
      cat("krige.conv: anisotropy correction performed\n")
    coords <- coords.aniso(coords = coords, aniso.pars = krige$aniso.pars)
    locations <- coords.aniso(coords = locations, aniso.pars = krige$aniso.pars)
  }
  ##
  ## Box-Cox transformation
  ##
  if(lambda != 1) {
    if(messages.screen)
      cat("krige.conv: Box-Cox's transformation of the data was performed.\n")
    if(lambda == 0)
      data <- log(data)
    else data <- ((data^lambda) - 1)/lambda
  }  
  ## 
  ## setting covariance parameters
  ##
  if (is.vector(cov.pars)) {
    sigmasq <- cov.pars[1]
    phi <- cov.pars[2]
  }
  else {
    sigmasq <- cov.pars[, 1]
    phi <- cov.pars[, 2]
  }
  sill.partial <- micro.scale + sum(sigmasq)
  sill.total <- nugget + sum(sigmasq)
  n <- length(data)
  ni <- nrow(trend.l)
  ##
  ## starting kriging calculations
  ##
  kc.result <- list()
### old code
  ##  invcov <- varcov.spatial(coords = coords, cov.model = cov.model, 
  ##                           kappa = kappa, nugget = nugget,
  ##                           cov.pars = cov.pars, inv = TRUE)$inverse
  ##  ttiv <- crossprod(trend.d, invcov)
  ##  ittivtt <- solve(ttiv %*% trend.d)
  invcov <- varcov.spatial(coords = coords, cov.model = cov.model, 
                               kappa = kappa, nugget = nugget,
                               cov.pars = cov.pars, inv = TRUE,
                               only.inv.lower.diag = TRUE)
  temp <- as.double(rep(0, beta.size * beta.size))
  temp <- .C("bilinearform_XAY",
             as.double(invcov$lower.inverse),
             as.double(invcov$diag.inverse),
             as.double(as.vector(trend.d)),
             as.double(as.vector(trend.d)),
             as.integer(beta.size),
             as.integer(beta.size),
             as.integer(n),
             res=temp)$res
  attr(temp, "dim") <- c(beta.size, beta.size)
  ittivtt <- solve(temp)
  remove("temp")
  if (beta.prior == "flat"){
### old code
    ##    beta.flat <- ittivtt %*% ttiv %*% data
    temp <- as.double(rep(0, beta.size))
    temp <- .C("bilinearform_XAY",
               as.double(invcov$lower.inverse),
               as.double(invcov$diag.inverse),
               as.double(as.vector(trend.d)),
               as.double(as.vector(data)),
               as.integer(beta.size),
               as.integer(1),
               as.integer(n),
               res=temp)$res
    beta.flat <- ittivtt %*% temp
    remove("temp")
    
  }
  v0mat <- as.double(rep(0, n*ni))
  .C("loccoords",
     as.double(as.vector(locations[,1])),
     as.double(as.vector(locations[,2])),
     as.double(as.vector(coords[,1])),
     as.double(as.vector(coords[,2])),
     as.integer(ni),
     as.integer(n),
     v0mat, DUP=FALSE)
  attr(v0mat, "dim") <- c(n, ni)
  if(n.sim > 0){
    ## checking data points coincident with prediction locations
    loc.coincide <- apply(v0mat, 2, function(x, min.dist){any(x < min.dist)},min.dist=krige$dist.epsilon)
    if(any(loc.coincide))
      loc.coincide <- (1:ni)[loc.coincide]
    else
      loc.coincide <- NULL
    if(!is.null(loc.coincide)){
      temp.f <- function(x, data, dist.eps){return(data[x < dist.eps])}
      data.coincide <- apply(v0mat[, loc.coincide, drop=FALSE], 2, temp.f, data=data, dist.eps=krige$dist.epsilon)
    }
    else
      data.coincide <- NULL
  }
  else
    remove("locations")
  if(signal){
    v0mat <- ifelse(v0mat < krige$dist.epsilon, sill.partial,
                    cov.spatial(obj = v0mat, cov.model = cov.model, 
                                kappa = kappa, cov.pars = cov.pars))
  }
  else{
    v0mat <- ifelse(v0mat < krige$dist.epsilon, sill.total,
                    cov.spatial(obj = v0mat, cov.model = cov.model, 
                                kappa = kappa, cov.pars = cov.pars))
  }   
### old code
  ##  tv0iv <- t(apply(v0mat, 2, crossprod, y = invcov))
  ##  remove("invcov")	
  ##  if(n.sim > 0) {
  ##    reduce.var <- tv0iv %*% v0mat
  ##    tv0ivv0 <- diag(reduce.var)
  ##  }
  ##  else
  ##    tv0ivv0 <- diag(tv0iv %*% v0mat)
  tv0ivv0 <- as.double(rep(0,ni))
  tv0ivv0 <- .C("diag_quadraticform_XAX",
                as.double(invcov$lower.inverse),
                as.double(invcov$diag.inverse),
                as.double(as.vector(v0mat)),
                as.integer(ni),
                as.integer(n),
                res = tv0ivv0)$res
### old code
  ##  tb <- t(trend.l) - apply(tv0iv, 1, crossprod, y = trend.d)
  ##  tb <- trend.l - tv0iv %*% trend.d
  tb <- as.double(rep(0, ni*beta.size))
  tb <- .C("bilinearform_XAY",
           as.double(invcov$lower.inverse),
           as.double(invcov$diag.inverse),
           as.double(as.vector(v0mat)),
           as.double(as.vector(trend.d)),
           as.integer(ni),
           as.integer(beta.size),
           as.integer(n),
           res=tb)$res
  attr(tb, "dim") <- c(ni, beta.size)
  tb <- trend.l - tb
###  
  if (beta.prior == "deg") {
### old code    
    ##    kc.result$predict <- as.vector((tv0iv %*% data) + (tb %*% beta))
    tv0ivdata <- as.double(rep(0,ni))
    tv0ivdata <- .C("bilinearform_XAY",
                    as.double(invcov$lower.inverse),
                    as.double(invcov$diag.inverse),
                    as.double(as.vector(v0mat)),
                    as.double(data),
                    as.integer(ni),
                    as.integer(1),
                    as.integer(n),
                    res = tv0ivdata)$res
    if(n.sim == 0) remove("v0mat","invcov")
    kc.result$predict <- tv0ivdata + as.vector(tb %*% beta)
    if(n.sim == 0) remove("tb")
    remove("tv0ivdata")
    if (krige$signal) 
      kc.result$krige.var <- as.vector(sill.partial - tv0ivv0)
    else kc.result$krige.var <- as.vector(sill.total - tv0ivv0)
    beta.est <- "Simple kriging performed (beta provided by user)"
  }
  if (beta.prior == "flat"){
### old    
    ##    kc.result$predict <- as.vector((tv0iv %*% data) + (tb %*% beta.flat))
    tv0ivdata <- as.double(rep(0,ni))
    tv0ivdata <- .C("bilinearform_XAY",
                    as.double(invcov$lower.inverse),
                    as.double(invcov$diag.inverse),
                    as.double(as.vector(v0mat)),
                    as.double(data),
                    as.integer(ni),
                    as.integer(1),
                    as.integer(n),
                    res = tv0ivdata)$res
    if(n.sim == 0) remove("v0mat","invcov")
    kc.result$predict <- tv0ivdata + as.vector(tb %*% beta.flat)
    remove("tv0ivdata")
### old code    
    ##    bi <- tb %*% ittivtt
    ##    if(n.sim > 0) {
    ##      ok.add.var <- bi %*% t(tb)
    ##      reduce.var <- reduce.var + ok.add.var
    ##      bitb <- diag(ok.add.var)
    ##    }
    ##    else
    ##      bitb <- diag(bi %*% t(tb))
    if(beta.size == 1)
      bitb <- as.vector(tb^2) * as.vector(ittivtt)
    else{
      bitb <- as.double(rep(0,ni))
      bitb <- .C("diag_quadraticform_XAX",
                 as.double(ittivtt[lower.tri(ittivtt)]),
                 as.double(diag(ittivtt)),
                 as.double(as.vector(t(tb))),
                 as.integer(ni),
                 as.integer(beta.size),
                 res = bitb)$res
    }
    if(n.sim == 0) remove("tb")
    if (krige$signal) 
      kc.result$krige.var <- as.vector(sill.partial - tv0ivv0 + bitb)
    else kc.result$krige.var <- as.vector(sill.total - tv0ivv0 + bitb)
    kc.result$beta.est <- beta.flat
    remove("bitb")
  }
  remove("tv0ivv0")
  if(any(round(kc.result$krige.var, dig=12) < 0))
    cat("krige.conv: negative kriging variance found! Investigate why this is happening.\n")
  message <- "krige.conv: Kriging performed using global neighbourhood"
  if(messages.screen)
    cat(paste(message, "\n"))
############## Sampling from the resulting distribution #####################
  if(n.sim > 0) {
    if(messages.screen)
      cat("krige.conv: sampling from the predictive distribution (conditional simulations)\n")
    if(length(cov.pars) > 2){
      reduce.var <- as.double(rep(0, ni * ni))
      .C("bilinearform_XAY",
         as.double(invcov$lower.inverse),
         as.double(invcov$diag.inverse),
         as.double(as.vector(v0mat)),
         as.double(as.vector(v0mat)),
         as.integer(ni),
         as.integer(ni),
         as.integer(n),
         reduce.var, DUP=FALSE)
      remove("v0mat")
      attr(reduce.var, "dim") <- c(ni, ni)
      if(beta.prior == "flat"){
        if(beta.size == 1)
          ok.add.var <- outer(as.vector(tb),as.vector(tb)) * as.vector(ittivtt)
        else{
          b <- t(tb)
          remove("tb")
          ok.add.var <- as.double(rep(0,ni*ni))
          .C("bilinearform_XAY",
             as.double(ittivtt[lower.tri(ittivtt)]),
             as.double(diag(ittivtt)),
             as.double(as.vector(b)),
             as.double(as.vector(b)),
             as.integer(ni),
             as.integer(ni),
             as.integer(beta.size),
             ok.add.var, DUP=FALSE)
          attr(ok.add.var, "dim") <- c(ni, ni)
          remove("b")
        }
        reduce.var <- reduce.var + ok.add.var
      }
      varcov <- varcov.spatial(coords = locations,
                               cov.model = cov.model,
                               cov.pars = cov.pars,
                               kappa = kappa, nugget = nugget)$varcov - reduce.var
      remove("reduce.var")
      if(is.R()) gc(verbose=FALSE)
      kc.result$simulations <-  kc.result$predict + crossprod(chol(varcov), matrix(rnorm(ni * n.sim), ncol=n.sim))
    }
    else{
      if(((round(1e12 * nugget) == 0) | signal) & (!is.null(loc.coincide))){
        v0mat <- v0mat[,-(loc.coincide)]
        nloc <- ni - length(loc.coincide)
        tmean.coincide <- kc.result$predict[loc.coincide]
        tmean <- kc.result$predict[-(loc.coincide)]
        tb <- tb[-(loc.coincide),]
      }
      else{
        nloc <- ni
        tmean <- kc.result$predict
      }
      normalsc <-  rnorm(nloc*n.sim)
      if (signal) Dval <- 1.0 + micro.scale
      else Dval <-  1.0 + (nugget/cov.pars[1])
      if(beta.size == 1){
        Blower <- 0
        if(beta.prior == "flat")
          Bdiag <- ittivtt
        else
          Bdiag <- 0.0
      }
      else{
        Blower <-  ittivtt[lower.tri(ittivtt)]
        if(beta.prior == "flat")
          Bdiag <-   diag(ittivtt)
        else
          Bdiag <- rep(0, beta.size)
      }
      R0 <- as.double(rep(0.0, (nloc*(nloc+1))/2))
      if(((round(1e12 * nugget) == 0) | signal) & (!is.null(loc.coincide))){
        .C("distdiag",
           as.double(locations[-(loc.coincide),1]),
           as.double(locations[-(loc.coincide),2]),
           as.integer(nloc),
           R0, DUP = FALSE)
      }
      else
        .C("distdiag",
           as.double(locations[,1]),
           as.double(locations[,2]),
           as.integer(ni),
           R0, DUP = FALSE)
      remove("locations")
      R0 <- cov.spatial(R0, cov.pars=cov.pars, cov.model=cov.model, kappa=kappa)
      normalsc <- .C("kb_sim",
                     as.double(tmean),
                     out = as.double(as.vector(normalsc)),
                     as.double(invcov$lower.inverse),
                     as.double(invcov$diag.inverse),
                     as.double(as.vector(v0mat)),
                     as.integer(nloc),
                     as.integer(n),
                     as.double(Dval),
                     as.integer(n.sim),
                     as.double(rep(1, n.sim)),                      
                     as.double(sill.partial),                      
                     as.double(Blower),
                     as.double(Bdiag),
                     as.double(as.vector(t(tb))),
                     as.integer(beta.size),
                     as.double(R0))$out
      attr(normalsc, "dim") <- c(nloc, n.sim)
      remove("v0mat", "R0", "tb", "invcov")
      if(((round(1e12 * nugget) == 0) | signal) & (!is.null(loc.coincide))){
        kc.result$simulations <- matrix(0, nrow=ni, ncol=n.sim)
        kc.result$simulations[-(loc.coincide),] <- normalsc
        kc.result$simulations[loc.coincide,] <- rep(tmean.coincide, n.sim)
      }
      else
        kc.result$simulations <- normalsc
      remove("normalsc")
    }
    if(lambda != 1){
      cat("krige.conv: back-transforming the simulated values\n")
      if(any(kc.result$simulations < -1/lambda))
        warning("Truncation in the back-transformation: there are simulated values less than (- 1/lambda) in the normal scale.")
      if(lambda == 0)
        kc.result$simulations <- ifelse(kc.result$simulations > -1/lambda, exp(kc.result$simulations), -1/lambda)
      if(lambda > 0)
        kc.result$simulations <- ifelse(kc.result$simulations > -1/lambda, ((kc.result$simulations*lambda) + 1)^(1/lambda), -1/lambda)
      if(lambda < 0)
        warning("back transformation not performed (negative value of lambda)")
    }
  }
########### Back - transforming predictions############################
  if(lambda != 1) {
    if(messages.screen)
      cat("krige.conv: back-transforming the predictions according to the (Box-Cox) parameter lambda\n")
    kc.result$transf.predict <- kc.result$predict
    kc.result$transf.krige.var <- kc.result$krige.var
    if(lambda == 0 & beta.prior == "deg") {
## don't change the order of the next two commands!!!
      kc.result$predict <- exp(kc.result$transf.predict + 0.5 * kc.result$krige.var)
      kc.result$krige.var <- (kc.result$predict^2) * (exp(kc.result$krige.var) -1)
    }
    if(lambda > 0 | (lambda == 0 & beta.prior == "flat")) {
      if(messages.screen)
        cat("krige.conv: back-transformation done by sampling from the resulting (normal) predictive distribution (inspect results carefully, run the function more than once and check for stability of the results\n")
      ap.warn <- options()$warn
      options(warn = -1)
      temp.data <- matrix(rnorm(ni * n.samples.backtransform,
                                mean = kc.result$transf.predict,
                                sd = sqrt(kc.result$transf.krige.var)),
                          nrow = ni)
      options(warn = ap.warn)
      ind.zero <- (round(1e12*kc.result$transf.krige.var) == 0)
      temp.data[ind.zero,  ] <- kc.result$transf.predict[ind.zero]
      remove(ind.zero)
      if(lambda == 0)
        temp.data <- exp(temp.data)
      else{
        temp.data[temp.data < -1/lambda] <- -1/lambda
        temp.data <- ((temp.data * lambda) + 1)^(1/lambda)
###      temp.data[is.na(temp.data)] <- Inf
      }
      kc.result$predict <- as.vector(apply(temp.data, 1, mean))
      kc.result$krige.var <- as.vector(apply(temp.data, 1, var))
    }
    if(lambda < 0) {
      cat("krige.conv: resulting distribution has no mean for lambda < 0 - back transformation not performed. Consider quantiles estimators\n"
          )
      kc.result$predict <- "back-transformation not performed"
      kc.result$krige.var <- "back-transformation not performed"
    }
  }
  kc.result <- c(kc.result, list(message = message, call = call.fc))
#####################################
  class(kc.result) <- "kriging"
  return(kc.result)
}

"krige.control" <-
  function (type.krige = "ok", beta = NULL,  
            trend.d = "cte", trend.l = "cte",
            cov.model = "matern",
            cov.pars = stop("covariance parameters (sigmasq and phi) should be provided"), kappa = 0.5,
            nugget = 0, micro.scale = 0, dist.epsilon = 1e-10, 
            aniso.pars = NULL, lambda = 1, 
            signal = FALSE,
            n.samples.backtransform = 500, n.sim = 0)
{
  cov.model <- match.arg(cov.model,
                         choices = c("matern", "exponential", "gaussian",
                           "spherical", "circular", "cubic", "wave", "power",
                           "powered.exponential", "cauchy", "gneiting",
                           "gneiting.matern", "pure.nugget"))
  return(list(type.krige = type.krige, beta = beta,
              trend.d = trend.d, trend.l = trend.l, 
              cov.model = cov.model, 
              cov.pars = cov.pars, kappa = kappa,
              nugget = nugget,
              micro.scale = micro.scale, dist.epsilon = dist.epsilon, 
              aniso.pars = aniso.pars, lambda = lambda,
              signal = signal,
              n.samples.backtransform = n.samples.backtransform,
              n.sim = n.sim))
}


            


