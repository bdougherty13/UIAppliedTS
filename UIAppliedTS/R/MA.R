#' ma()
#' @export
#' @description computes a simple moving average smoother of a given time series.
#' @param x Univariate time series
#' @param order Order of moving average smoother
#' @param centre If TRUE, then the moving average is centered for even orders.
#' @author Rob J Hydman.



ma<-function (x, order, centre = TRUE)
{
  if (abs(order - round(order)) > 1e-08) {
    stop("order must be an integer")
  }
  if (order%%2 == 0 && centre) {
    w <- c(0.5, rep(1, order - 1), 0.5)/order
  }
  else {
    w <- rep(1, order)/order
  }
  return(filter(x, w))
}
