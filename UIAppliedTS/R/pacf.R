#' pacf

#' @description Used for the partial autocorrelations.
#'
#' pacf uses code from the [`pacf`](/library/stats/help/pacf) function found in the [`stats`](/library/stats/help/stats-package) package.
#' @param lag.max maximum number of lags at which to calculate the acf. Default is 10*log10(N/m) where N is the number of observations and m the number of series.
#' @param type character string giving the type of acf to be computed. Allowed values are "correlation" (the default), "covariance" or "partial".
#' @param plot logical. If TRUE (the default) the acf is plotted.
#' @param na.action function to be called to handle missing values. na.pass can be used.
#' @param demean logical. Should the covariances be about the sample means?
#' @param drop.lag.0 logical. Should lag 0 be dropped.
#' @inherit stats::pacf
#' @export

pacf<-function (x, lag.max = NULL, plot = TRUE, na.action = na.fail,
                 ...)
{
  series <- deparse(substitute(x))
  x <- drop(na.action(as.ts(x)))
  if (!is.numeric(x))
    stop("'x' must be numeric")
  x.freq <- frequency(x)
  sampleT <- NROW(x)
  if (is.null(lag.max))
    lag.max <- if (is.matrix(x))
      floor(10 * (log10(sampleT) - log10(ncol(x))))
  else floor(10 * (log10(sampleT)))
  lag.max <- min(lag.max, sampleT - 1)
  if (lag.max < 1)
    stop("'lag.max' must be at least 1")
  if (is.matrix(x)) {
    if (anyNA(x))
      stop("NAs in 'x'")
    nser <- ncol(x)
    x <- sweep(x, 2, colMeans(x), check.margin = FALSE)
    lag <- matrix(1, nser, nser)
    lag[lower.tri(lag)] <- -1
    pacf <- ar.yw(x, order.max = lag.max)$partialacf
    lag <- outer(1L:lag.max, lag/x.freq)
    snames <- colnames(x)
  }
  else {
    x <- scale(x, TRUE, FALSE)
    acf <- drop(acf(x, lag.max = lag.max, plot = FALSE, na.action = na.action)$acf)
    pacf <- .Call(C_pacf1, acf, lag.max)
    lag <- array((1L:lag.max)/x.freq, dim = c(lag.max, 1L,
                                              1L))
    snames <- NULL
  }
  acf.out <- structure(.Data = list(acf = pacf, type = "partial",
                                    n.used = sampleT, lag = lag, series = series, snames = snames),
                       class = "acf")
  if (plot) {
    plot.acf(acf.out, ...)
    invisible(acf.out)
  }
  else acf.out
}
