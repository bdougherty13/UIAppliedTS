% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EACF.R
\name{EACF}
\alias{EACF}
\title{EACF}
\usage{
EACF(
  x,
  lag.max = NULL,
  type = c("correlation", "covariance", "partial"),
  plot = TRUE,
  na.action = na.fail,
  demean = TRUE,
  z,
  ...
)
}
\arguments{
\item{x}{a time series.}

\item{z}{Included for backward compatability}

\item{ar.max}{maximum AR order; default=7}

\item{ma.max}{maximum MA order; default=13}
}
\description{
Computes the sample extended acf (ESACF) for the time series stored in z. The matrix of ESACF with the AR order up to ar.max and the MA order up to ma.max is stored in the matrix EACFM. This function is very similar to the eacf function in TSA, with changes having been made to some of the argument names.
}
