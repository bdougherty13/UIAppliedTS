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

#' predict.ar
#' @export
#' @description Creates a times series prediction based off of an AR object.

predict.ar <- function(object, newdata, n.ahead = 1L, se.fit = TRUE, ...)
{
  if (n.ahead < 1L) stop("'n.ahead' must be at least 1")
  if(missing(newdata)) {
    newdata <- eval.parent(parse(text=object$series))
    if (!is.null(nas <- object$call$na.action))
      newdata <- eval.parent(call(nas, newdata))
  }
  nser <- NCOL(newdata)
  ar <- object$ar
  p <- object$order
  st <- tsp(as.ts(newdata))[2L]
  dt <- deltat(newdata)
  xfreq <- frequency(newdata)
  tsp(newdata) <- NULL
  class(newdata) <- NULL
  if(NCOL(ar) != nser)
    stop("number of series in 'object' and 'newdata' do not match")
  n <- NROW(newdata)
  if(nser > 1L) {
    if(is.null(object$x.intercept)) xint <- rep.int(0, nser)
    else xint <- object$x.intercept
    x <- rbind(sweep(newdata, 2L, object$x.mean, check.margin = FALSE),
               matrix(rep.int(0, nser), n.ahead, nser, byrow = TRUE))
    pred <- if(p) {
      for(i in seq_len(n.ahead)) {
        x[n+i,] <- ar[1L,,] %*% x[n+i-1L,] + xint
        if(p > 1L) for(j in 2L:p)
          x[n+i,] <- x[n+i,] + ar[j,,] %*% x[n+i-j,]
      }
      x[n + seq_len(n.ahead), ]
    } else matrix(xint, n.ahead, nser, byrow = TRUE)
    pred <- pred + matrix(object$x.mean, n.ahead, nser, byrow = TRUE)
    colnames(pred) <- colnames(object$var.pred)
    if(se.fit) {
      warning("'se.fit' not yet implemented for multivariate models")
      se <- matrix(NA, n.ahead, nser)
    }
  } else {
    if(is.null(object$x.intercept)) xint <- 0
    else xint <- object$x.intercept
    x <- c(newdata - object$x.mean, rep.int(0, n.ahead))
    if(p) {
      for(i in seq_len(n.ahead))
        x[n+i] <- sum(ar * x[n+i - seq_len(p)]) + xint
      pred <- x[n + seq_len(n.ahead)]
      if(se.fit) {
        psi <- .Call(C_ar2ma, ar, n.ahead - 1L)
        vars <- cumsum(c(1, psi^2))
        se <- sqrt(object$var.pred*vars)[seq_len(n.ahead)]
      }
    } else {
      pred <- rep.int(xint, n.ahead)
      if (se.fit) se <- rep.int(sqrt(object$var.pred), n.ahead)
    }
    pred <- pred + rep.int(object$x.mean, n.ahead)
  }
  pred <- ts(pred, start = st + dt, frequency = xfreq)
  if(se.fit)
    list(pred = pred, se = ts(se, start = st + dt, frequency = xfreq))
  else pred
}
