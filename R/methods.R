#' @keywords internal
predtidy <- function(x)UseMethod("predtidy") 

#' @keywords internal
cverror <- function(x, y = NULL, ...)UseMethod("cverror")

#' @export
modfit <- function(object, return_data = FALSE, ...)UseMethod("modfit")

#' @export
modtidy <- function(x)UseMethod("modtidy") 

#' @export
get_avesurv <- function(object, ...)UseMethod("get_avesurv")

#' @export
get_indivsurv <- function(object, newdata)UseMethod("get_indivsurv")

#' @keywords internal
survconcord <- function(object, newdata = NULL, stats = FALSE, ...)UseMethod("survconcord")

#' @keywords internal
pvimp <- function(model, newdata, nrep = 20, nclusters = parallel::detectCores(), ...)UseMethod("pvimp")

#' @keywords internal
varimp <- function(object, type = c("coef", "perm", "model"), relative = TRUE, newdata, nrep = 20, nclusters = parallel::detectCores(), ...)UseMethod("varimp")
