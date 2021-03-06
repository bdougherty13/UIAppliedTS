% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ar.burg.R
\name{ar.burg}
\alias{ar.burg}
\title{ar.burg}
\usage{
ar.burg(x, ...)
}
\arguments{
\item{x}{a univariate or multivariate time series.}

\item{...}{additional arguments for specific methods.}

\item{aic}{logical. If TRUE then the Akaike Information Criterion is used to choose the order of the autoregressive model. If FALSE, the model of order ar.max is fitted.}

\item{ar.max}{maximum order (or order) of model to fit. Defaults to the smaller of N-1 and 10*log10(N) where N is the number of non-missing observations except for method = "mle" where it is the minimum of this quantity and 12.}

\item{method}{character string specifying the method to fit the model. Must be one of the strings in the default argument (the first few characters are sufficient). Defaults to "yule-walker".}

\item{na.action}{function to be called to handle missing values. Currently, via na.action = na.pass, only Yule-Walker method can handle missing values which must be consistent within a time point: either all variables must be missing or none.}

\item{demean}{should a mean be estimated during fitting?}

\item{series}{names for the series. Defaults to deparse (subsititute (x)).}

\item{var.method}{the method to estimate the innovations variance (see ‘Details’).}

\item{get.aic}{should the aic of the final model be returned?}

\item{object}{a fit from ar().}

\item{newdata}{data to which to apply the prediction.}

\item{n.ahead}{number of steps ahead at which to predict.}

\item{se.fit}{logical: return estimated standard errors of the prediction error?}

\item{order.max}{included for backward compatability.}
}
\value{
The coefficients of an Autoregressive model, its sigma^2, and possibly the value of its AIC.
}
\description{
Fit an autoregressive time series model to the data using the Burg method, by default selecting the complexity by AIC.
}
\author{
ar.burg by B.D. Ripley based on R version by Martyn Plummer. Minor revisions by Brad Dougherty.
}
