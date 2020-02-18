#' ARMASUBSETS
#' @export
#' @description This function finds a number of subset ARMA models. A "long" AR model is fitted to the data y to compute the residuals which are taken as a proxy of the error process. Then, an ARMA model is approximated by a regression model with the the covariates being the lags of the time series and the lags of the error process. Subset ARMA models may then be selected using the subset regression technique by leaps and bounds, via the regsubsets function of the leaps package in R. This function is very similar to the armasubsets function in TSA, with changes having been made to some of the argument names.
#' @param x time-series data.
#' @param ar.max maximum AR order.
#' @param ma.max maximum MA order.
#' @param y.name label of the time series.
#' @param method method used for fitting the long AR model; default is ols with the AR order determined by AIC.
#' @param ... arguments passed to the plot.armasubsets function
#' @param y included for backward compatability.
#' @param nar included for backward compatability.
#' @param nma included for backward compatability.
#' @author This code was orignially implemented by Kung-Sik Chan. Minor changes have been made by Brad Dougherty.

ARMASUBSETS<-function (x, ar.max, ma.max, y.name = "Y", ar.method = "ols",y,
                       ...)


{
  #code for backward compatibility
  if(!is.null(x)&!is.null(y))return("must supply only x or y")
  if(!is.null(x)) y <- x

  if(!is.null(nar)&!is.null(ar.max))return("must supply only ar.max or nar")
  if(!is.null(nar))  ar.max <- nar

  if(!is.null(nma)&!is.null(ma.max))return("must supply only ma.max or nma")
  if(!is.null(nma))  ma.max <- nma


  lab = NULL
  if (ar.max > 1)
    lab = c(lab, paste(y.name, 1:ar.max, sep = "-lag"))
  if (ma.max > 1)
    lab = c(lab, paste("error", 1:ma.max, sep = "-lag"))
  res.ar = ar(y, method = ar.method)
  resid = res.ar$resid
  x = NULL
  if (ar.max > 1) {
    for (i in 1:ar.max) {
      x = cbind(x, zlag(y, d = i))
    }
  }
  if (ma.max > 1) {
    for (j in 1:ma.max) {
      x = cbind(x, zlag(resid, d = j))
    }
  }
  x = na.omit(cbind(y, x))
  y = x[, 1]
  x = x[, -1]
  x = data.frame(x)
  colnames(x) = lab
  regobj = leaps::regsubsets(y ~ ., data = x, ...)
  class(regobj) = c("armasubsets", "regsubsets")
  invisible(regobj)
}
