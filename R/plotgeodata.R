plot.geodata <- function (x, coords = x$coords, data = x$data, borders = NULL, 
    trend = "cte", lambda = 1, col.data = 1, weights.divide = NULL, 
    window.new = FALSE, ...) 
{
  if(missing(x)) x <- list(coords=coords, data = data)
  if (is.R()) 
        par.ori <- par(no.readonly = TRUE)
    else par.ori <- par()
    on.exit(par(par.ori))
    coords <- as.matrix(coords)
    data <- as.matrix(data)
    data <- data[, col.data]
    if (window.new) {
        if (is.R()) 
            X11()
        else trellis.device()
    }
    if (!is.null(weights.divide)) {
        if (length(weights.divide) != length(data)) 
            stop("length of weights.divide must be equals to the length of data")
        data <- data/weights.divide
    }
    if (lambda != 1) {
        if (lambda == 0) 
            data <- log(data)
        else data <- ((data^lambda) - 1)/lambda
    }
    xmat <- unclass(trend.spatial(trend = trend, geodata = x))
    if (nrow(xmat) != nrow(coords)) 
      stop("coords and trend have incompatible sizes")
    if (trend != "cte") {
      data <- lm(data ~ xmat + 0)$residuals
      names(data) <- NULL
    }
    par(mfrow = c(2, 2))
    par(mar = c(4, 4, 0, 0.5))
    data.quantile <- quantile(data)
    if (is.null(borders)) 
        coords.lims <- set.coords.lims(coords = coords)
    else {
        if (ncol(borders) != 2) 
            stop("argument \"borders\" must be a 2 columns object with coordinates of the borders of the study area")
        coords.lims <- set.coords.lims(coords = rbind(coords, 
            borders))
    }
    par(pty = "s")
    plot(coords, xlab = "Coord X", ylab = "Coord Y", type = "n", 
        xlim = coords.lims[, 1], ylim = coords.lims[, 2])
    if (!is.null(borders)) 
        polygon(borders)
    if (is.R()) {
        data.breaks <- unique(quantile(data))
        n.breaks <- length(data.breaks)
        data.cut <- cut(data, breaks = data.breaks, include.l = TRUE, 
            labels = FALSE)
        points(coords, pch = (1:4)[data.cut], col = c("blue", 
            "green", "yellow2", "red")[data.cut])
    }
    else {
        points(coords[(data <= data.quantile[2]), ], pch = 1, 
            cex = 0.6, col = 2)
        points(coords[((data > data.quantile[2]) & (data <= data.quantile[3])), 
            ], pch = 2, cex = 1.4, col = 4)
        points(coords[((data > data.quantile[3]) & (data <= data.quantile[4])), 
            ], pch = 3, cex = 1.7, col = 7)
        points(coords[(data > data.quantile[4]), ], pch = 4, 
            cex = 2, col = 8)
    }
    plot(data, coords[, 2], ylab = "Coord Y", cex = 1, ylim = coords.lims[, 2])
    if (!is.R()) 
        par(mar = c(5, 5, 1, 0.5))
    plot(coords[, 1], data, xlab = "Coord X", cex = 1, xlim = coords.lims[, 1], )
    par(pty = "m")
    if (is.R()) 
        par(mar = c(4, 4, 1, 1))
    else par(mar = c(0, 1, 0, 0.5))
    if (is.R()) {
        if (require(scatterplot3d) == FALSE) {
            hist(data)
            cat("plot.geodata: a 3d plot would be drawn instead of the histogram if the package \"scatterplot3d\" was available\n")
        }
        else scatterplot3d(x = coords[, 1], y = coords[, 2], 
            z = data, box = FALSE, type = "h", pch = 16, xlab = "Coord X", 
            ylab = "Coord Y", ...)
    }
    else xyzplot(coords = coords, data = data, ...)
    return(invisible())
}
