#' Predictions
#'
#' @export
predict.satpred <- function(object, ...) {
	new_args <- list(...)
	new_args$object <- object
	pred <- do.call("predict", new_args)
	return(pred)
}

#' Average survival
#'
#' @export
get_avesurv.rfsrc <- function(object, ...) {
	object <- modtidy(object)
#	surv <- as.vector(sapply(data.frame(object$surv), mean))
	surv <- as.vector(colMeans(data.frame(object$surv)))
	# chaz <- as.vector(sapply(data.frame(object$chaz), mean))
	chaz <- -log(surv)
	time <- object$time
	out <- list(time = time, surv = surv, chaz=chaz)
	out$call <- match.call()
	class(out) <- "satsurv"
	out
}

#' Average survival for glmnetsurv
#'
#' @export

get_avesurv.glmnetsurv <- function(object, ...) {
	pred <- glmnetsurvfit(object, ...)
	surv <- rowMeans(pred$surv)
	chaz <- -log(surv)
	time <- pred$time
	out <- list(time = time, surv = surv, chaz=chaz)
	out$call <- match.call()
	class(out) <- "satsurv"
	out
}

#' Individual survival
#'
#' @export
get_indivsurv.rfsrc <- function(object, newdata) {
	pred <- predict(object, newdata = newdata)
	out <- predtidy(pred)
	out <- list(time = out$time, surv = out$surv, chaz = out$chaz)
	out$call <- match.call()
	class(out) <- "satsurv"
	out
}

#' Predict survival probabilities at various time points
#'
#' The function extracts the survival probability predictions from a \code{satpred} model.
#'
#' @aliases predictSurvProb
#'
#'
#' @return a matrix of probabilities with as many rows as the rows of the \code{newdata} and as many columns as number of time points (\code{times}).
#'
#'
#' @importFrom prodlim sindex
#' @importFrom pec predictSurvProb
#' @export predictSurvProb
#' @export

predictSurvProb.satpred <- function(object, newdata, times, ...){
	N <- NROW(newdata)
	sfit <- get_indivsurv(object, newdata = newdata)
	S <- sfit$surv
	Time <- sfit$time
	if(N == 1) S <- matrix(S, nrow = 1)
	p <-  cbind(1, S)[, 1 + prodlim::sindex(Time, times),drop = FALSE]
	if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
		stop("Prediction failed")
	p
}

#' Predict survival probabilities at various time points
#'
#' The function extracts the survival probability predictions from a \code{satpred} model.
#'
#' @aliases predictSurvProb
#'
#'
#' @return a matrix of probabilities with as many rows as the rows of the \code{newdata} and as many columns as number of time points (\code{times}).
#'
#'
#' @importFrom prodlim sindex
#' @importFrom pec predictSurvProb
#' @export predictSurvProb
#' @export

predictSurvProb.gbm.satpred <- function(object, newdata, times, ...){
	N <- NROW(newdata)
	sfit <- get_indivsurv(object, newdata = newdata)
	S <- sfit$surv
	Time <- sfit$time
	if(N == 1) S <- matrix(S, nrow = 1)
	p <-  cbind(1, S)[, 1 + prodlim::sindex(Time, times),drop = FALSE]
	if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
		stop("Prediction failed")
	p
}

#' Extract predictions from satpred model
#'
#' Extract event probabilities from the fitted model.
#'
#' @aliases predictRisk
#'
#' @details
#' For survival outcome, the function predicts the risk, \eqn{1 - S(t|x)}, where \eqn{S(t|x)} is the survival chance of an individual characterized by \eqn{x}.
#'
#' @importFrom riskRegression predictRisk
#' @export predictRisk
#' @export

predictRisk.satpred <- function(object, newdata, times, ...){
	p <- 1 - predictSurvProb.satpred(object, newdata, times)
	p
}

#' Compute the concordance statistic for the pcoxtime model
#'
#' The function computes the agreement between the observed response and the predictor.
#'
#' @export

get_survconcord <- function(object, newdata = NULL, stats = FALSE, ...) {
	concord <- survconcord(object, newdata, stats, ...)
	return(concord)
}

#' Permutation variable importance
#'
#' Computes the relative importance based on random permutation of focal variable for various survival models.
#'
#' @export

get_pvimp <- function(model, newdata, nrep = 20, parallelize = TRUE, nclusters = parallel::detectCores(), ...) {
	vi <- pvimp(model, newdata, nrep, parallelize = parallelize, nclusters = nclusters, ...)
	return(vi)
}

#' Compute variable importance of various survival models object
#'
#' @aliases get_varimp
#'
#' @details
#' Absolute value of the coefficients (coef) corresponding the tuned model are used \code{type = perm}. Otherwise, variable level importance is computed from the model. 
#'
#' @export

get_varimp <- function(object, type = c("coef", "perm", "model"), relative = TRUE, newdata, nrep = 20
	, modelname = "model1", parallelize = TRUE, nclusters = parallel::detectCores(), ...) {
	imp <- varimp(object, type, relative, newdata, nrep, parallelize = parallelize, nclusters = nclusters, ...)
	imp$terms <- rownames(imp)
	rownames(imp) <- NULL
	out <- imp[, c("terms", "Overall", "sign")]
	out$model <- modelname
	class(out) <- c("varimp", class(out))
	attr(out, "estimate") <- "mean"
	return(out)
}
