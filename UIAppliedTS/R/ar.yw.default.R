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

#' ar.yw.default
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

ar.yw.default <- function (x, aic = TRUE, order.max = NULL, na.action = na.fail,
            demean = TRUE, series = NULL,get.aic=TRUE,ar.max=NULL, ...)
  {

  #code for backward compatability
  if(!is.null(order.max)&!is.null(ar.max))return("must supply only order.max or ar.max")
  if(!is.null(ar.max)) order.max <- ar.max
    if (is.null(series))
      series <- deparse(substitute(x))
    ists <- is.ts(x)
    x <- na.action(as.ts(x))
    if (ists)
      xtsp <- tsp(x)
    xfreq <- frequency(x)
    x <- as.matrix(x)
    if (!is.numeric(x))
      stop("'x' must be numeric")
    if (any(is.na(x) != is.na(x[, 1])))
      stop("NAs in 'x' must be the same row-wise")
    nser <- ncol(x)
    if (demean) {
      xm <- colMeans(x, na.rm = TRUE)
      x <- sweep(x, 2L, xm, check.margin = FALSE)
    }
    else xm <- rep.int(0, nser)
    n.used <- nrow(x)
    n.obs <- sum(!is.na(x[, 1]))
    order.max <- if (is.null(order.max))
      min(n.obs - 1L, floor(10 * log10(n.obs)))
    else floor(order.max)
    if (order.max < 1L)
      stop("'order.max' must be >= 1")
    else if (order.max >= n.obs)
      stop("'order.max' must be < 'n.obs'")
    xacf <- acf(x, type = "covariance", lag.max = order.max,
                plot = FALSE, demean = demean, na.action = na.pass)$acf
    if (nser > 1L) {
      snames <- colnames(x)
      A <- B <- array(0, dim = c(order.max + 1L, nser, nser))
      A[1L, , ] <- B[1L, , ] <- diag(nser)
      EA <- EB <- xacf[1L, , , drop = TRUE]
      partialacf <- array(dim = c(order.max, nser, nser))
      xaic <- numeric(order.max + 1L)
      solve.yw <- function(m) {
        betaA <- betaB <- 0
        for (i in 0L:m) {
          betaA <- betaA + A[i + 1L, , ] %*% xacf[m + 2L -
                                                    i, , ]
          betaB <- betaB + B[i + 1L, , ] %*% t(xacf[m +
                                                      2L - i, , ])
        }
        KA <- -t(qr.solve(t(EB), t(betaA)))
        KB <- -t(qr.solve(t(EA), t(betaB)))
        EB <<- (diag(nser) - KB %*% KA) %*% EB
        EA <<- (diag(nser) - KA %*% KB) %*% EA
        Aold <- A
        Bold <- B
        for (i in seq_len(m + 1L)) {
          A[i + 1L, , ] <<- Aold[i + 1L, , ] + KA %*% Bold[m +
                                                             2L - i, , ]
          B[i + 1L, , ] <<- Bold[i + 1L, , ] + KB %*% Aold[m +
                                                             2L - i, , ]
        }
      }
      cal.aic <- function(m) {
        logdet <- determinant.matrix(EA)$modulus
        n.obs * logdet + 2 * m * nser * nser
      }
      cal.resid <- function() {
        resid <- array(0, dim = c(n.used - order, nser))
        for (i in 0L:order) resid <- resid + tcrossprod(x[(order -
                                                             i + 1L):(n.used - i), , drop = FALSE], ar[i +
                                                                                                         1L, , ])
        rbind(matrix(NA, order, nser), resid)
      }
      order <- 0L
      for (m in 0L:order.max) {
        xaic[m + 1L] <- cal.aic(m)
        if (!aic || xaic[m + 1L] == min(xaic[seq_len(m +
                                                     1L)])) {
          ar <- A
          order <- m
          var.pred <- EA * n.obs/(n.obs - nser * (m + 1L))
        }
        if (m < order.max) {
          solve.yw(m)
          partialacf[m + 1L, , ] <- -A[m + 2L, , ]
        }
      }
      min.aic<-min(xaic)
      xaic <- setNames(xaic - min(xaic), 0L:order.max)
      resid <- cal.resid()
      if (order) {
        ar <- -ar[2L:(order + 1L), , , drop = FALSE]
        dimnames(ar) <- list(seq_len(order), snames, snames)
      }
      else ar <- array(0, dim = c(0L, nser, nser), dimnames = list(NULL,
                                                                   snames, snames))
      dimnames(var.pred) <- list(snames, snames)
      dimnames(partialacf) <- list(seq_len(order.max), snames,
                                   snames)
      colnames(resid) <- colnames(x)
    }
    else {
      if (xacf[1L] == 0)
        stop("zero-variance series")
      r <- as.double(drop(xacf))
      z <- .Fortran(stats:::C_eureka, as.integer(order.max), r, r,
                    coefs = double(order.max^2), vars = double(order.max),
                    double(order.max))
      coefs <- matrix(z$coefs, order.max, order.max)
      partialacf <- array(diag(coefs), dim = c(order.max, 1L,
                                               1L))
      var.pred <- c(r[1L], z$vars)
      xaic <- n.obs * log(var.pred) + 2 * (0L:order.max) +
        2 * demean
      maic <- min(aic)
      min.aic <- min(xaic)
      xaic <- setNames(if (is.finite(maic))
        xaic - min(xaic)
        else ifelse(xaic == maic, 0, Inf), 0L:order.max)
      order <- if (aic)
        (0L:order.max)[xaic == 0L]
      else order.max
      ar <- if (order)
        coefs[order, seq_len(order)]
      else numeric()
      var.pred <- var.pred[order + 1L]
      var.pred <- var.pred * n.obs/(n.obs - (order + 1L))
      resid <- if (order)
        c(rep.int(NA, order), embed(x, order + 1L) %*% c(1,
                                                         -ar))
      else as.vector(x)
      if (ists) {
        attr(resid, "tsp") <- xtsp
        attr(resid, "class") <- "ts"
      }
    }
    res <- list(order = order, ar = ar, var.pred = var.pred,
                x.mean = drop(xm), aic = xaic, n.used = n.used, n.obs = n.obs,
                order.max = order.max, partialacf = partialacf, resid = resid,
                method = "Yule-Walker", series = series, frequency = xfreq,
                call = match.call(),model.aic=min.aic,get.aic=get.aic)
    if (nser == 1L && order)
      res$asy.var.coef <- var.pred/n.obs * solve(toeplitz(drop(xacf)[seq_len(order)]))

  res
  }
