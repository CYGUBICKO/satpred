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
get_avesurv.satpred <- function(object, ...) {
	object <- modtidy(object)
	surv <- sapply(data.frame(object$surv), mean)
	chaz <- sapply(data.frame(object$chaz), mean)
	time <- object$time
	df <- data.frame(time = time, surv = surv, chaz=chaz)
	rownames(df) <- NULL
	return(df)
}

#' Individual survival
get_indivsurv.satpred <- function(object, newdata) {
	pred <- predict(object, newdata = newdata)
	pred <- predtidy(pred)
	return(pred)
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

get_survconcord <- function(object, newdata = NULL, stats = FALSE) {
	concord <- survconcord(object, newdata, stats)
	return(concord)
}
