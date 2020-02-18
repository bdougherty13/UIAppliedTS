#'AR
#' @description Fit an autoregressive time series model to the data, by default selecting the complexity by AIC.
#'
#' AR uses code from the [`ar`](/library/stats/help/ar) function found in the [`stats`](/library/stats/help/stats-package) package. The function has been edited to allow for the retrieval of its AIC value, and changes have been made to the names of its arguments.
#' @param x 	a univariate or multivariate time series.
#' @param aic logical. If TRUE then the Akaike Information Criterion is used to choose the order of the autoregressive model. If FALSE, the model of order order.max is fitted.
#' @param ar.max maximum order (or order) of model to fit. Defaults to the smaller of N-1 and 10*log10(N) where N is the number of non-missing observations except for method = "mle" where it is the minimum of this quantity and 12.
#' @param method character string specifying the method to fit the model. Must be one of the strings in the default argument (the first few characters are sufficient). Defaults to "yule-walker".
#' @param na.action function to be called to handle missing values. Currently, via na.action = na.pass, only Yule-Walker method can handle missing values which must be consistent within a time point: either all variables must be missing or none.
#' @param demean should a mean be estimated during fitting?
#' @param series names for the series. Defaults to deparse (subsititute (x)).
#' @param var.method the method to estimate the innovations variance (see ‘Details’).
#' @param get.aic logical. If TRUE the aic of the final model will be returned.
#' @param ... additional arguments for specific methods.
#' @param object a fit from ar().
#' @param newdata data to which to apply the prediction.
#' @param n.ahead number of steps ahead at which to predict.
#' @param se.fit logical: return estimated standard errors of the prediction error?
#' @param order.max included for backward compatibility
#' @examples AR(lh)

#' AR(lh, method = "burg")
#' AR(lh, method = "ols")
#' AR(lh, FALSE, 4) # fit ar(4)
#'
#'(sunspot.ar <- AR(sunspot.year))
#' predict(sunspot.ar, n.ahead = 25)
#' ## try the other methods too
#'
#' AR(ts.union(BJsales, BJsales.lead))
#' ## Burg is quite different here, as is OLS (see ar.ols)
#' AR(ts.union(BJsales, BJsales.lead), method = "burg")

#' @inherit stats::ar
#' @author Martyn Plummer. Univariate case of ar.yw, ar.mle and C code for univariate case of ar.burg by B. D. Ripley. Minor revisions to R code by Brad Dougherty.
#' @return The coefficients of an Autoregressive model, its sigma^2, and possibly the value of its AIC.
#' @importFrom knitr knit_print

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




AR<-function (x, aic = TRUE,get.aic=NULL, ar.max = NULL,     #added  get.aic=TRUE, changed order.max to ar.max
            method = c("yule-walker","burg", "ols", "mle", "yw"),
            na.action = na.fail, series = deparse(substitute(x)), order.max=NULL,...)
  {
    res <- switch(match.arg(method),
                  yw =,
                  "yule-walker" = ar.yw(x, aic = aic,get.aic=get.aic, ar.max = ar.max,  #added get.aic=get.aic
                                        na.action = na.action, series = series, order.max = order.max, ...),
                  "burg" = ar.burg(x, aic = aic,get.aic=get.aic, ar.max = ar.max,      #added get.aic=get.aic
                                   na.action = na.action, series = series, order.max = order.max,...),
                  "ols" = ar.ols(x, aic = aic,get.aic=get.aic, ar.max = ar.max,      #added get.aic=get.aic
                                 na.action = na.action, series = series, order.max = order.max,...),
                  "mle" = ar.mle(x, aic = aic,get.aic=get.aic, ar.max = ar.max,       #added get.aic=get.aic
                                 na.action = na.action, series = series,order.max = order.max, ...)
    )
    res$call <- match.call()
    class(res)<-"AR"
  res
  }


