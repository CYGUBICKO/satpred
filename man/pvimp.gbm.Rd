% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gbm_satpred.R
\name{pvimp.gbm}
\alias{pvimp.gbm}
\title{Permutation variable importance method for gbm}
\usage{
\method{pvimp}{gbm}(
  model,
  newdata,
  nrep = 20,
  parallelize = TRUE,
  nclusters = parallel::detectCores(),
  estimate = c("mean", "quantile"),
  probs = c(0.025, 0.5, 0.975),
  seed = NULL,
  ...
)
}
\description{
Permutation variable importance method for gbm
}
\keyword{internal}
