#' ar.burg
#' @description Fit an autoregressive time series model to the data using the Burg method, by default selecting the complexity by AIC.
#' @param x	a univariate or multivariate time series.
#' @param aic logical. If TRUE then the Akaike Information Criterion is used to choose the order of the autoregressive model. If FALSE, the model of order ar.max is fitted.
#' @param ar.max maximum order (or order) of model to fit. Defaults to the smaller of N-1 and 10*log10(N) where N is the number of non-missing observations except for method = "mle" where it is the minimum of this quantity and 12.
#' @param method character string specifying the method to fit the model. Must be one of the strings in the default argument (the first few characters are sufficient). Defaults to "yule-walker".
#' @param na.action function to be called to handle missing values. Currently, via na.action = na.pass, only Yule-Walker method can handle missing values which must be consistent within a time point: either all variables must be missing or none.
#' @param demean should a mean be estimated during fitting?
#' @param series names for the series. Defaults to deparse (subsititute (x)).
#' @param var.method the method to estimate the innovations variance (see ‘Details’).
#' @param get.aic should the aic of the final model be returned?
#' @param ... additional arguments for specific methods.
#' @param object a fit from ar().
#' @param newdata data to which to apply the prediction.
#' @param n.ahead number of steps ahead at which to predict.
#' @param se.fit logical: return estimated standard errors of the prediction error?
#' @param order.max included for backward compatability.
#' @author ar.burg by B.D. Ripley based on R version by Martyn Plummer. Minor revisions by Brad Dougherty.

#' @return The coefficients of an Autoregressive model, its sigma^2, and possibly the value of its AIC.





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




## ar.burg by B.D. Ripley based on R version by Martyn Plummer



ar.burg <- function(x, ...) UseMethod("ar.burg")
ar.burg.default <-
  function (x, aic = TRUE, na.action = na.fail,     #added get.aic=true
            demean = TRUE, series = NULL, var.method = 1L,ar.max=NULL, order.max = NULL,get.aic=get.aic, ...)
  {
    #code for backward compatability
    if(!is.null(order.max)&!is.null(ar.max))return("must supply only order.max or ar.max")
    if(!is.null(order.max)) ar.max <- order.max


    if(is.null(series)) series <- deparse(substitute(x))
    if (ists <- is.ts(x)) xtsp <- tsp(x)
    x <- na.action(as.ts(x))
    if(anyNA(x)) stop("NAs in 'x'")
    if (ists)  xtsp <- tsp(x)
    xfreq <- frequency(x)
    x <- as.vector(x) # drop attributes including class
    if (demean) {
      x.mean <- mean(x)
      x <- x - x.mean
    } else x.mean <- 0
    n.used <- length(x)
    ar.max <- if (is.null(ar.max))
      min(n.used-1L, floor(10 * log10(n.used)))
    else floor(ar.max)
    if (ar.max < 1L) stop("'ar.max' must be >= 1")
    else if (ar.max >= n.used) stop("'ar.max' must be < 'n.used'")
    xaic <- numeric(ar.max + 1L)
    z <- .Call(stats:::C_Burg, x, ar.max)
    coefs <- matrix(z[[1L]], ar.max, ar.max)
    partialacf <- array(diag(coefs), dim = c(ar.max, 1L, 1L))
    var.pred <- if(var.method == 1L) z[[2L]] else z[[3L]]
    if (any(is.nan(var.pred))) stop("zero-variance series")
    xaic <- n.used * log(var.pred) + 2 * (0L:ar.max) + 2 * demean   # note this leads to very different value than the aic found by creating the equivalant model in arima()
    maic <- min(aic)

    min.aic<-min(xaic)
    xaic <- setNames(if(is.finite(maic)) xaic - min(xaic) else
      ifelse(xaic == maic, 0, Inf), 0L:ar.max)
    order <- if (aic) (0L:ar.max)[xaic == 0] else ar.max
    ar <- if (order) coefs[order, 1L:order] else numeric()
    var.pred <- var.pred[order + 1L]
    resid <- if(order) c(rep(NA, order), embed(x, order+1L) %*% c(1, -ar))
    else x
    if(ists) {
      attr(resid, "tsp") <- xtsp
      attr(resid, "class") <- "ts"
    }
    res <- list(order = order, ar = ar, var.pred = var.pred, x.mean = x.mean,model.aic=min.aic,
                aic.list = xaic, n.used = n.used, ar.max = ar.max,
                partialacf = partialacf, resid = resid,
                method = ifelse(var.method==1L,"Burg","Burg2"),
                series = series, frequency = xfreq, call = match.call(),get.aic=get.aic)
    if(order) {
      xacf <- acf(x, type = "covariance", lag.max = order, plot = FALSE)$acf
      res$asy.var.coef <- solve(toeplitz(drop(xacf)[seq_len(order)]))*var.pred/n.used
    }
    class(res) <- "ar"
    res
  }
