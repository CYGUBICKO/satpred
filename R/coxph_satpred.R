
#' Implemented packages
#'
#' @export
coxph.satpred <- function(formula = NULL, train_df = NULL, test_df = NULL, param_grid = NULL, method = c("efron","breslow","exact"), finalmod = FALSE, ...) {

	factor_to_characters <- function(x) {
   	factor_cols = sapply(x, is.factor)
    	x[factor_cols] = lapply(x[factor_cols], as.character)
    	x
	}
	coxph_args <- list(formula=formula, data=train_df)
	if (is.null(param_grid)) {
		param <- expand.grid(method = method, stringsAsFactors = FALSE)
	} else {
		param <- factor_to_characters(param_grid)
	}
	param_args <- as.list(param)
	coxph_args[names(param_args)] <- param_args
	new_args <- list(...)
	if (length(new_args)) coxph_args[names(new_args)] <- new_args
	
	if (!finalmod) {
		args_match <- match(colnames(param), names(coxph_args), nomatch = FALSE)
		param_match <- match(names(coxph_args), colnames(param), nomatch = FALSE)
		error <- lapply(1:NROW(param), function(x){
			coxph_args[args_match] <- param[x, param_match]
			fit <- do.call(survival::coxph, coxph_args)
			if (is.null(test_df)) test_df <- train_df
			pred <- predict(fit, test_df)
			class(pred) <- c(class(pred), "coxph")
			all_params <- names(param_args)
			all_params <- union("method", all_params)
			param_temp <- fit[all_params]
			names(param_temp) <- all_params
			y  <- model.extract(model.frame(formula, data = test_df), "response")
			error_list <- list(param_temp, error = 1-cverror(pred, y))
			error_df <- as.data.frame(error_list)
			return(error_df)
		})
		error <- do.call("rbind", error)
		return(error)
	} else {
		coxph_args[["x"]] <- TRUE
		fit <- do.call(survival::coxph, coxph_args)
		return(fit)
	}
}

#' Cross-validation error
#'
#' Method for coxph
#'
#' @keywords internal

cverror.coxph <- function(x, y = NULL, ...){
	score <- survival::concordancefit(y, -x)$concordance
	return(score)
}


#' Cross-validation plots
#'
#' @import ggplot2
#' @export

plot.coxph.satpred <- function(x, ..., show_best = TRUE, lsize = 0.3, pshape = "O") {
	tune_df <- x$result
	best_df <- x$besTune
p1 <- (ggplot(tune_df, aes(x = reorder(method, -error), y = error, group=1))
	+ geom_point()
	+ geom_linerange(aes(ymin = min(error), ymax = error), size = 0.3)
	+ labs(x = "Method", y = "Error (1 - C)")
)
	if (show_best) {
		p1 <- (p1
			+ geom_point(data=best_df, aes(x = method, y = error), colour="red", size=2)
			+ geom_hline(data=best_df, aes(yintercept=error), lty=2)
		)
	}
	return(p1)
}

#' Average survival for coxph
#'
#' @export
get_avesurv.coxph <- function(object, ...) {
	pred <- survfit(object, se=FALSE, ...)
	surv <- rowMeans(pred$surv)
	chaz <- -log(surv)
	time <- pred$time
	out <- list(time = time, surv = surv, chaz=chaz)
	out$call <- match.call()
	class(out) <- "satsurv"
	out
}

#' Average survival
#'
#' @export
get_avesurv.coxph <- function(object, ...) {
	object <- get_indivsurv(object)
	surv <- as.vector(colMeans(object$surv))
	chaz <- -log(surv)
	time <- object$time
	out <- list(time = time, surv = surv, chaz=chaz)
	out$call <- match.call()
	class(out) <- "satsurv"
	out
}

#' Individual survival
#'
#' @export
get_indivsurv.coxph <- function(object, newdata) {
	out <- survfit(object, newdata = newdata)
	out <- list(time = out$time, surv = t(out$surv), chaz = t(out$cumhaz))
	out$call <- match.call()
	class(out) <- "satsurv"
	out
}

