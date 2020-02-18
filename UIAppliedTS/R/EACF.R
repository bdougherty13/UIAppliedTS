#' EACF
#' @export
#' @description Computes the sample extended acf (ESACF) for the time series stored in z. The matrix of ESACF with the AR order up to ar.max and the MA order up to ma.max is stored in the matrix EACFM. This function is very similar to the eacf function in TSA, with changes having been made to some of the argument names.
#' @param x 	a time series.
#' @param ar.max maximum AR order; default=7
#' @param ma.max maximum MA order; default=13
#' @param z Included for backward compatability



EACF<-function (x, lag.max = NULL, type = c("correlation", "covariance","partial"),
                plot = TRUE, na.action = na.fail, demean = TRUE,z,...)

#This code was taken from The University of Iowa's, Dr. Kung-Sik Chan's eacf() function in his TSA package.

{
  if(!is.null(x)&!is.null(z))return("must only supply x or z")
  if(!is.null(z)) x <- z  #included for backward compatability.

  type <- match.arg(type)
  if (type == "partial") {
    m <- match.call()
    m[[1L]] <- quote(stats::pacf)
    m$type <- NULL
    return(eval(m, parent.frame()))
  }
  series <- deparse(substitute(x))
  x <- na.action(as.ts(x))
  x.freq <- frequency(x)
  x <- as.matrix(x)
  if (!is.numeric(x))
    stop("'x' must be numeric")
  sampleT <- as.integer(nrow(x))
  nser <- as.integer(ncol(x))
  if (is.na(sampleT) || is.na(nser))
    stop("'sampleT' and 'nser' must be integer")
  if (is.null(lag.max))
    lag.max <- floor(10 * (log10(sampleT) - log10(nser)))
  lag.max <- as.integer(min(lag.max, sampleT - 1L))
  if (is.na(lag.max) || lag.max < 0)
    stop("'lag.max' must be at least 0")
  if (demean)
    x <- sweep(x, 2, colMeans(x, na.rm = TRUE), check.margin = FALSE)
  lag <- matrix(1, nser, nser)
  lag[lower.tri(lag)] <- -1
  acf <- .Call(C_acf, x, lag.max, type == "correlation")
  lag <- outer(0:lag.max, lag/x.freq)
  acf.out <- structure(list(acf = acf, type = type, n.used = sampleT,
                            lag = lag, series = series, snames = colnames(x)), class = "acf")
  if (plot) {
    plot.acf(acf.out, ...)
    invisible(acf.out)
  }
  else acf.out
}
