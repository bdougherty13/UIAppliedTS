#' ar.burg.mts
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
#' @export

ar.burg.mts<-function (x, aic = TRUE, ar.max = NULL, na.action = na.fail,
          demean = TRUE, series = NULL, var.method = 1L, get.aic=TRUE, order.max=NULL, ...)
{
  #code for backward compatability
  if(!is.null(order.max)&!is.null(ar.max))return("must supply only order.max or ar.max")
  if(!is.null(order.max)) ar.max <- order.max

  if (is.null(series))
    series <- deparse(substitute(x))
  if (ists <- is.ts(x))
    xtsp <- tsp(x)
  x <- na.action(as.ts(x))
  if (anyNA(x))
    stop("NAs in 'x'")
  if (ists)
    xtsp <- tsp(x)
  xfreq <- frequency(x)
  x <- as.matrix(x)
  nser <- ncol(x)
  n.used <- nrow(x)
  if (demean) {
    x.mean <- colMeans(x)
    x <- sweep(x, 2L, x.mean, check.margin = FALSE)
  }
  else x.mean <- rep(0, nser)
  ar.max <- floor(if (is.null(ar.max)) 10 * log10(n.used) else ar.max)
  z <- .C(stats:::C_multi_burg, as.integer(n.used), resid = as.double(x),
          as.integer(ar.max), as.integer(nser), coefs = double((1L +
                                                                     ar.max) * nser * nser), pacf = double((1L + ar.max) *
                                                                                                                nser * nser), var = double((1L + ar.max) * nser *
                                                                                                                                             nser), aic = double(1L + ar.max), order = integer(1L),
          as.integer(aic), as.integer(var.method))
  partialacf <- aperm(array(z$pacf, dim = c(nser, nser, ar.max +
                                              1L)), 3:1)[-1L, , , drop = FALSE]
  var.pred <- aperm(array(z$var, dim = c(nser, nser, ar.max +
                                           1L)), 3:1)
  min.aic <- min(z$aic)
  xaic <- setNames(z$aic - min(z$aic), 0:ar.max)
  order <- z$order
  ar <- if (order)
    -aperm(array(z$coefs, dim = c(nser, nser, ar.max +
                                    1L)), 3:1)[2L:(order + 1L), , , drop = FALSE]
  else array(dim = c(0, nser, nser))
  var.pred <- var.pred[order + 1L, , , drop = TRUE]
  resid <- matrix(z$resid, nrow = n.used, ncol = nser)
  if (order)
    resid[seq_len(order), ] <- NA
  if (ists) {
    attr(resid, "tsp") <- xtsp
    attr(resid, "class") <- "mts"
  }
  snames <- colnames(x)
  colnames(resid) <- snames
  dimnames(ar) <- list(seq_len(order), snames, snames)
  dimnames(var.pred) <- list(snames, snames)
  dimnames(partialacf) <- list(seq_len(ar.max), snames,
                               snames)
  res <- list(order = order, ar = ar, var.pred = var.pred,
              x.mean = x.mean, aic = xaic, n.used = n.used, n.obs = n.used,
              ar.max = ar.max, partialacf = partialacf, resid = resid,
              method = ifelse(var.method == 1L, "Burg", "Burg2"),
              series = series, frequency = xfreq, model.aic = min.aic,get.aic=get.aic, call = match.call())
  class(res) <- "ar"
  res
}
