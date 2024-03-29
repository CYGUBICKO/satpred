% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gbm_satpred.R
\name{gbm3Hazard}
\alias{gbm3Hazard}
\title{Compute predicted hazard for gbm3}
\usage{
gbm3Hazard(
  y,
  x = NULL,
  wt = rep(1, NROW(y)),
  risk = NULL,
  survtype = NULL,
  vartype = NULL
)
}
\arguments{
\item{object}{fitted \code{\link[gbm3]{gbm}} object.}
}
\value{
A list of S3 objects.
\item{n}{number of observations used in the fit.}
\item{events}{total number of events of interest in the fit.}
\item{time}{time points defined by the risk set.}
\item{n.risk}{the number of individuals at risk at time \code{t}.}
\item{n.event}{the number of events that occur at time \code{t}.}
\item{n.censor}{the number of subjects who exit the risk set, without an event, at time \code{t}.}
\item{surv}{a vector or a matrix of estimated survival function.}
\item{chaz, hazard}{a vector or a matrix of estimated cumulative hazard.}
}
\description{
This code is borrowed from internal function agsurv from survival package.
}
\keyword{internal}
