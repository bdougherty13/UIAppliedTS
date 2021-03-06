% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AR.R
\name{AR}
\alias{AR}
\title{AR}
\usage{
AR(
  x,
  aic = TRUE,
  get.aic = NULL,
  ar.max = NULL,
  method = c("yule-walker", "burg", "ols", "mle", "yw"),
  na.action = na.fail,
  series = deparse(substitute(x)),
  order.max = NULL,
  ...
)
}
\arguments{
\item{x}{a univariate or multivariate time series.}

\item{aic}{logical. If TRUE then the Akaike Information Criterion is used to choose the order of the autoregressive model. If FALSE, the model of order order.max is fitted.}

\item{get.aic}{logical. If TRUE the aic of the final model will be returned.}

\item{ar.max}{maximum order (or order) of model to fit. Defaults to the smaller of N-1 and 10*log10(N) where N is the number of non-missing observations except for method = "mle" where it is the minimum of this quantity and 12.}

\item{method}{character string specifying the method to fit the model. Must be one of the strings in the default argument (the first few characters are sufficient). Defaults to "yule-walker".}

\item{na.action}{function to be called to handle missing values. Currently, via na.action = na.pass, only Yule-Walker method can handle missing values which must be consistent within a time point: either all variables must be missing or none.}

\item{series}{names for the series. Defaults to deparse (subsititute (x)).}

\item{order.max}{included for backward compatibility}

\item{...}{additional arguments for specific methods.}

\item{demean}{should a mean be estimated during fitting?}

\item{var.method}{the method to estimate the innovations variance (see ‘Details’).}

\item{object}{a fit from ar().}

\item{newdata}{data to which to apply the prediction.}

\item{n.ahead}{number of steps ahead at which to predict.}

\item{se.fit}{logical: return estimated standard errors of the prediction error?}
}
\value{
The coefficients of an Autoregressive model, its sigma^2, and possibly the value of its AIC.
}
\description{
Fit an autoregressive time series model to the data, by default selecting the complexity by AIC.

AR uses code from the \href{/library/stats/help/ar}{\code{ar}} function found in the \href{/library/stats/help/stats-package}{\code{stats}} package. The function has been edited to allow for the retrieval of its AIC value, and changes have been made to the names of its arguments.
}
\details{
For definiteness, note that the AR coefficients have the sign in

  \deqn{x_t - \mu = a_1(x_{t-1} - \mu) + \cdots +  a_p(x_{t-p} - \mu) + e_t}{x[t] - m = a[1]*(x[t-1] - m) + \\dots +  a[p]*(x[t-p] - m) + e[t]}

  \code{ar} is just a wrapper for the functions \code{ar.yw},
  \code{ar.burg}, \code{\link[stats]{ar.ols}} and \code{ar.mle}.

  Order selection is done by AIC if \code{aic} is true. This is
  problematic, as of the methods here only \code{ar.mle} performs
  true maximum likelihood estimation. The AIC is computed as if the variance
  estimate were the MLE, omitting the determinant term from the
  likelihood. Note that this is not the same as the Gaussian likelihood
  evaluated at the estimated parameter values.  In \code{ar.yw} the
  variance matrix of the innovations is computed from the fitted
  coefficients and the autocovariance of \code{x}.

  \code{ar.burg} allows two methods to estimate the innovations
  variance and hence AIC. Method 1 is to use the update given by
  the Levinson-Durbin recursion (Brockwell and Davis, 1991, (8.2.6)
  on page 242), and follows S-PLUS. Method 2 is the mean of the sum
  of squares of the forward and backward prediction errors
  (as in Brockwell and Davis, 1996, page 145). Percival and Walden
  (1998) discuss both. In the multivariate case the estimated
  coefficients will depend (slightly) on the variance estimation method.

  Remember that \code{ar} includes by default a constant in the model, by
  removing the overall mean of \code{x} before fitting the AR model,
  or (\code{ar.mle}) estimating a constant to subtract.
}
\examples{
AR(lh)
AR(lh, method = "burg")
AR(lh, method = "ols")
AR(lh, FALSE, 4) # fit ar(4)

(sunspot.ar <- AR(sunspot.year))
predict(sunspot.ar, n.ahead = 25)
## try the other methods too

AR(ts.union(BJsales, BJsales.lead))
## Burg is quite different here, as is OLS (see ar.ols)
AR(ts.union(BJsales, BJsales.lead), method = "burg")
}
\references{
Brockwell, P. J. and Davis, R. A. (1991).
  \emph{Time Series and Forecasting Methods}, second edition.
  Springer, New York.
  Section 11.4.

  Brockwell, P. J. and Davis, R. A. (1996).
  \emph{Introduction to Time Series and Forecasting}.
  Springer, New York.
  Sections 5.1 and 7.6.

  Percival, D. P. and Walden, A. T. (1998).
  \emph{Spectral Analysis for Physical Applications}.
  Cambridge University Press.

  Whittle, P. (1963).
  On the fitting of multivariate autoregressions and the approximate
  canonical factorization of a spectral density matrix.
  \emph{Biometrika}, \bold{40}, 129--134.
  \ifelse{{text}{doi: 10.2307/2333753 (URL: https://doi.org/10.2307/2333753)}{\ifelse{{latex}{doi:\out{\\nobreakspace{}}\href{https://doi.org/10.2307/2333753}{10.2307\out{\\slash{}}2333753}}{doi: \href{https://doi.org/10.2307/2333753}{10.2307/2333753}}}}}.
}
\seealso{
\code{\link[stats]{ar.ols}}, \code{\link[stats]{arima}} for ARMA models;
  \code{\link[stats]{acf2AR}}, for AR construction from the ACF.

  \code{\link[stats]{arima.sim}} for simulation of AR processes.
}
\author{
Martyn Plummer. Univariate case of ar.yw, ar.mle and C code for univariate case of ar.burg by B. D. Ripley. Minor revisions to R code by Brad Dougherty.
}
