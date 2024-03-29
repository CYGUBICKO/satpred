% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{varimp.default}
\alias{varimp.default}
\title{Compute variable importance for default models}
\usage{
\method{varimp}{default}(
  object,
  type = c("coef", "perm", "model"),
  relative = TRUE,
  newdata,
  nrep = 20,
  parallelize = TRUE,
  nclusters = parallel::detectCores(),
  ...
)
}
\description{
Compute variable importance for default models
}
\keyword{internal}
