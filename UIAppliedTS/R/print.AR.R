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

#' print.AR
#' @description Prints output from the AR function.
#' @export

print.AR <- function(x, digits = max(3L, getOption("digits") - 3L), ...)
{
  cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
  nser <- NCOL(x$var.pred)
  if(nser > 1L) {
    res <- x[c("ar", if(!is.null(x$x.intercept)) "x.intercept", "var.pred")]
    res$ar <- aperm(res$ar, c(2L,3L,1L))
    print(res, digits = digits)
  } else { ## univariate case
    if(x$order) {
      cat("Coefficients:\n")
      coef <- setNames(round(drop(x$ar), digits = digits),
                       seq_len(x$order))
      print.default(coef, print.gap = 2L)
    }
    if(!is.null(xint <- x$x.intercept) && !is.na(xint))
      cat("\nIntercept: ", format(xint, digits = digits),
          ## FIXME? asy.se.coef  *only* exists for  ar.ols (??)
          " (", format(x$asy.se.coef$x.mean, digits = digits),
          ") ", "\n", sep = "")
    cat("\nOrder selected", x$order, " sigma^2 estimated as ",
        format(x$var.pred, digits = digits))
    if(!is.null(x$get.aic)) cat("\nAIC =",format(x$model.aic, digits = digits))
    cat("\n")

  }
  invisible(x)
}
