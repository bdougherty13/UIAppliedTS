#  Part of the R package, https://www.R-project.org
#
#  Copyright (C) 1999-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  https://www.R-project.org/Licenses/
#
#Much of this package was taken from STATS, TSA, and Forcast


#' ar.yw.mts
#' @description Fit an autoregressive time series model to the data using the Yule-Walker method, by default selecting the complexity by AIC.
#' @param x A univariate or multivariate time series.
#' @param aic Logical flag. If TRUE then the Akaike Information Criterion is used to choose the order of the autoregressive model. If FALSE, the model of order order.max is fitted.
#' @param ar.max Maximum order (or order) of model to fit. Defaults to 10*log10(N) where N is the number of observations.
#' @param na.action function to be called to handle missing values.
#' @param demean should the AR model be for x minus its mean?
#' @param intercept should a seperate intercept term be fitted?
#' @param series names for teh series. Defulats to deparse (subsitute(x)).
#' @param further arguments to be passed to or from methods.
#' @param order.max included for backwards compatibility.
ar.yw.mts <-function (x, aic = TRUE, get.aic=TRUE,order.max = NULL, na.action = na.fail,       #added get.aic=true
            demean = TRUE, series = NULL, var.method = 1L,ar.max=NULL, ...)
  {
    #code for backward compatability
    if(!is.null(order.max)&!is.null(ar.max))return("must supply only order.max or ar.max")
    if(!is.null(ar.max)) order.max <- ar.max


    if (is.null(series)) series <- deparse(substitute(x))
    if (ists <- is.ts(x)) xtsp <- tsp(x)
    x <- na.action(as.ts(x))
    if (anyNA(x)) stop("NAs in 'x'")
    if (ists) xtsp <- tsp(x)
    xfreq <- frequency(x)
    x <- as.matrix(x)
    nser <- ncol(x)
    n.used <- nrow(x)
    if (demean) {
      x.mean <- colMeans(x)
      x <- sweep(x, 2L, x.mean, check.margin=FALSE)
    }
    else x.mean <- rep(0, nser)
    order.max <- if (is.null(order.max)) floor(10 * log10(n.used)) else floor(order.max)
    if (order.max < 1L)
      stop("'order.max' must be >= 1")
    xacf <- acf(x, type = "cov", plot = FALSE, lag.max = order.max)$acf
    z <- .C(stats:::C_multi_yw,
            aperm(xacf, 3:1),
            as.integer(n.used),
            as.integer(order.max),
            as.integer(nser),
            coefs = double((1L + order.max) * nser * nser),
            pacf = double((1L + order.max) * nser * nser),
            var = double((1L + order.max) * nser * nser),
            aic = double(1L + order.max),
            order = integer(1L),
            as.integer(aic))
    partialacf <- aperm(array(z$pacf, dim = c(nser, nser, order.max + 1L)), 3:1)[-1L, , , drop = FALSE]
    var.pred <- aperm(array(z$var, dim = c(nser, nser, order.max + 1L)), 3:1)
    min.aic<-min(z$aic)
    xaic <- setNames(z$aic - min(z$aic), 0:order.max)
    order <- z$order
    resid <- x
    if (order > 0) {
      ar <- -aperm(array(z$coefs, dim = c(nser, nser, order.max + 1L)), 3:1)[2L:(order + 1L), , , drop = FALSE]
      for (i in 1L:order)
        resid[-(1L:order), ] <- resid[-(1L:order),] - x[(order - i + 1L):(n.used - i), ] %*% t(ar[i, , ])
      resid[1L:order, ] <- NA
    }
    else ar <- array(dim = c(0, nser, nser))
    var.pred <- var.pred[order + 1L, , , drop = TRUE] * n.used/(n.used - nser * (demean + order))
    if (ists) {
      attr(resid, "tsp") <- xtsp
      attr(resid, "class") <- c("mts", "ts")
    }
    snames <- colnames(x)
    colnames(resid) <- snames
    dimnames(ar) <- list(seq_len(order), snames, snames)
    dimnames(var.pred) <- list(snames, snames)
    dimnames(partialacf) <- list(1L:order.max, snames, snames)
    res <- list(order = order, ar = ar, var.pred = var.pred,
                x.mean = x.mean, aic = xaic,model.aic=min.aic, n.used = n.used, order.max = order.max,
                partialacf = partialacf, resid = resid, method = "Yule-Walker",
                series = series, frequency = xfreq, call = match.call(),get.aic=get.aic)
    class(res) <- "ar"
    return(res)
  }
