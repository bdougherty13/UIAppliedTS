#' BOXCOX
#' @export
#' @description Determine the appropriate power transformation for time-series data. The objective is to estimate the power transformation so that the transformed time series is approximately a Gaussian AR process. This function is very similar to the armasubsets function in TSA, with changes having been made to some of the argument names.

#' @param x univariate time series (must be positive)
#' @param order AR order for the data; if missing, the order is determined by AIC for the log-transformed data
#' @param lambda a vector of candidate power transformation values; if missing, it is set to be from -2 to 2, with increment .01
#' @param plotit logical value, if true, plot the profile log-likelihood for the power estimator
#' @param method method of AR estimation; default is "mle"
#' @param ... other parameters to be passed to the ar function
#' @param y Included for backward compatability.
#' @author Kung-Sik Chan, minor revisions were made by Brad Dougherty.

BOXCOX.AR<-function (x, order, lambda = seq(-2, 2, 0.01), plotit = TRUE,
                     method = c("mle", "yule-walker", "burg",
                                "ols", "yw"), y,...)
  #This code was taken from The University of Iowa's, Dr. Kung-Sik Chan's BoxCox.ar() function in his TSA package.

{
  #included for standardization
  if(!is.null(x)&!is.null(y))return("must supply only x or y")
if(!is.null(x)) y <- x
  x <- NULL


  if (missing(method))
    method = "mle"
  y = as.vector(y/(max(abs(y)) + 1))
  if (any(y <= 0))
    stop("Data values must be positive")
  order = ar(log(y), method = "mle")$order
  nlngmy <- sum(log(y))
  if (!missing(lambda))
    xl <- lambda
  else xl <- seq(-2, 2, 0.1)
  loglik <- as.vector(xl)
  for (i in 1:length(xl)) if (abs(xl[i]) > 0) {
    if (missing(order))
      ar.result = ar((y^xl[i] - 1)/xl[i], method = method)
    else ar.result = ar((y^xl[i] - 1)/xl[i], method = method,
                        order.max = order)
    n = length(y) - ar.result$order
    ar.res = ar.result$resid
    n = length(y)
    loglik[i] <- -n/2 * log(ar.result$var.pred) + (xl[i] -
                                                     1) * nlngmy
  }
  else {
    if (missing(order))
      ar.result = ar(log(y), method = method)
    else ar.result = ar(log(y), method = method, order.max = order)
    n = length(y) - ar.result$order
    ar.res = ar.result$resid
    n = length(y)
    loglik[i] <- -n/2 * log(ar.result$var.pred) - nlngmy
  }
  if (plotit) {
    plot(xl, loglik, xlab = expression(lambda), ylab = "Log Likelihood",
         type = "l", ylim = c(min(loglik), max(loglik)))
    lambdahat <- loglik[loglik == max(loglik)]
    limit <- lambdahat - 0.5 * qchisq(0.95, 1)
    in.interval = xl[loglik >= limit]
    lower = in.interval[1]
    upper = rev(in.interval)[1]
    mle = (xl[loglik == max(loglik)])[1]
    lines(x = c(lower, lower), y = c(min(loglik), limit),
          lty = 2)
    lines(x = c(upper, upper), y = c(min(loglik), limit),
          lty = 2)
    lines(x = c(mle, mle), y = c(min(loglik), max(loglik)),
          lty = 2)
    abline(limit, 0, lty = 2)
    scal <- (par("usr")[4] - par("usr")[3])/par("pin")[2]
    text(c(xl[1]) + 0.1, limit + 0.08 * scal, " 95%")
  }
  invisible(list(lambda = xl, loglike = loglik, mle = mle,
                 ci = c(lower, upper)))
}
