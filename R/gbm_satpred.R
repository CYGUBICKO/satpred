#' Implemented packages
#'
#' @export
gbm.satpred <- function(formula = NULL, train_df = NULL, test_df = NULL, distribution = "coxph", param_grid = NULL, n.trees = 1000, interaction.depth = 1, n.minobsinnode = 10, shrinkage = 0.1, finalmod = FALSE, error.method=c("auto", "OOB", "cv", "test"), ...) {
	
	gbm_args <- list(formula=formula, data=train_df, distribution=distribution)
	if (is.null(param_grid)) {
		if (is.null(shrinkage)) {
			param <- expand.grid(n.trees=n.trees, n.minobsinnode=n.minobsinnode, stringsAsFactors=FALSE)
		} else if (is.null(interaction.depth)) {
			param <- expand.grid(n.trees=n.trees, shrinkage=shrinkage, n.minobsinnode=n.minobsinnode, stringsAsFactors=FALSE)
		} else {
			param <- expand.grid(n.trees=n.trees, shrinkage=shrinkage, interaction.depth=interaction.depth, n.minobsinnode=n.minobsinnode, stringsAsFactors=FALSE)
		}
	} else {
		param <- param_grid
	}
	param_args <- as.list(param)
	gbm_args[names(param_args)] <- param_args
	new_args <- list(...)
	if (length(new_args)) gbm_args[names(new_args)] <- new_args

	error.method <- match.arg(error.method)
	if (!finalmod) {
		args_match <- match(colnames(param), names(gbm_args), nomatch = FALSE)
		param_match <- match(names(gbm_args), colnames(param), nomatch = FALSE)
		error <- lapply(1:NROW(param), function(x){
			gbm_args[args_match] <- param[x, param_match]
			fit <- do.call(gbm::gbm, gbm_args)
			if (is.null(test_df)) test_df <- train_df
			suppressMessages(	
				if (error.method=="auto") {
					if (fit$train.fraction<1) {
						..n.trees <- gbm::gbm.perf(fit, method="test", plot.it=FALSE)
					} else if (fit$cv.folds>1 & fit$train.fraction==1) {
						..n.trees <- gbm::gbm.perf(fit, method="cv", plot.it=FALSE)
					} else {
						..n.trees <- gbm::gbm.perf(fit, method="OOB", plot.it=FALSE)
					}
				} else {
						..n.trees <- gbm::gbm.perf(fit, method=error.method, plot.it=FALSE)
				}
			)
			fit$n.trees <- ..n.trees
			pred <- predict(fit, test_df, fit$n.trees)
			class(pred) <- c(class(pred), "gbm")
			all_params <- names(param_args)
			all_params <- union(c("shrinkage", "n.trees", "interaction.depth", "n.minobsinnode"), all_params)
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
		fit <- do.call(gbm::gbm, gbm_args)
		return(fit)
	}
}

#' Implemented packages
#'
#' @export
gbm3.satpred <- function(formula = NULL, train_df = NULL, test_df = NULL, distribution = "coxph", param_grid = NULL, n.trees = 1000, interaction.depth = 1, n.minobsinnode = 10, shrinkage = 0.1, finalmod = FALSE, error.method=c("auto", "OOB", "cv", "test"), ...) {
	
	gbm_args <- list(formula=formula, data=train_df, distribution=distribution)
	if (is.null(param_grid)) {
		if (is.null(shrinkage)) {
			param <- expand.grid(n.trees=n.trees, n.minobsinnode=n.minobsinnode, stringsAsFactors=FALSE)
		} else if (is.null(interaction.depth)) {
			param <- expand.grid(n.trees=n.trees, shrinkage=shrinkage, n.minobsinnode=n.minobsinnode, stringsAsFactors=FALSE)
		} else {
			param <- expand.grid(n.trees=n.trees, shrinkage=shrinkage, interaction.depth=interaction.depth, n.minobsinnode=n.minobsinnode, stringsAsFactors=FALSE)
		}
	} else {
		param <- param_grid
	}
	param_args <- as.list(param)
	gbm_args[names(param_args)] <- param_args
	new_args <- list(...)
	if (length(new_args)) gbm_args[names(new_args)] <- new_args

	if (!finalmod) {
		args_match <- match(colnames(param), names(gbm_args), nomatch = FALSE)
		param_match <- match(names(gbm_args), colnames(param), nomatch = FALSE)
		error <- lapply(1:NROW(param), function(x){
			gbm_args[args_match] <- param[x, param_match]
			fit <- do.call(gbm3::gbm, gbm_args)
			suppressMessages(	
				fit <- gbm3::to_old_gbm(fit)
			)
			if (is.null(test_df)) test_df <- train_df
			suppressMessages(	
				if (error.method=="auto") {
					if (fit$train.fraction<1) {
						..n.trees <- gbm::gbm.perf(fit, method="test", plot.it=FALSE)
					} else if (fit$cv.folds>1 & fit$train.fraction==1) {
						..n.trees <- gbm::gbm.perf(fit, method="cv", plot.it=FALSE)
					} else {
						..n.trees <- gbm::gbm.perf(fit, method="OOB", plot.it=FALSE)
					}
				} else {
						..n.trees <- gbm::gbm.perf(fit, method=error.method, plot.it=FALSE)
				}
			)
			fit$n.trees <- ..n.trees
			pred <- predict(fit, test_df, fit$n.trees)
			class(pred) <- c(class(pred), "gbm")
			all_params <- names(param_args)
			all_params <- union(c("shrinkage", "n.trees", "interaction.depth", "n.minobsinnode"), all_params)
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
		fit <- do.call(gbm3::gbm, gbm_args)
		suppressMessages(	
			fit <- gbm3::to_old_gbm(fit)
		)
		return(fit)
	}
}


#' Cross-validation error
#'
#' Method for gbm
#'
#' @keywords internal

cverror.gbm <- function(x, y = NULL, ...){
#	score <- Hmisc::rcorr.cens(-x, y)[[1]]
	score <- survival::concordancefit(y, -x)$concordance
	return(score)
}


#' Compute the concordance statistic for the pcoxtime model
#'
#' The function computes the agreement between the observed response and the predictor.
#'
#' @keywords internal

survconcord.gbm <- function(object, newdata = NULL, stats = FALSE, ...) {
	if (is.null(newdata)) newdata <- object$modelData
	pred <- predict(object, newdata=newdata, n.trees=object$n.trees)
	class(pred) <- c(class(pred), "gbm")
	y  <- model.extract(model.frame(terms(object), data = newdata), "response")
	concord <- cverror(pred, y)
	return(concord)
}



#' Cross-validation plots
#'
#' @import ggplot2
#' @export

plot.gbm.satpred <- function(x, ..., show_best = TRUE, lsize = 0.3, pshape = "O") {
	tune_df <- x$result
	tune_df$n.minobsinnode <- factor(tune_df$n.minobsinnode, labels=paste0("nodesize: ", unique(tune_df$n.minobsinnode)))
	best_df <- x$besTune
	best_df$n.minobsinnode <- factor(best_df$n.minobsinnode, labels=paste0("nodesize: ", unique(best_df$n.minobsinnode)))
	p1 <- (ggplot(tune_df, aes(x = as.factor(n.trees), y = error, group=as.factor(interaction.depth), colour = as.factor(interaction.depth)))
		+ geom_point(shape = pshape)
		+ geom_line(size = lsize)
		+ facet_grid(shrinkage~n.minobsinnode)
		+ labs(x = "# boosting iterations", y = "Error (1 - C)", colour = "Max tree depth")
	)
	if (show_best) {
		p1 <- (p1
			+ geom_point(data=best_df, aes(x = as.factor(n.trees), y = error), colour="red", size=2)
			+ geom_hline(data=best_df, aes(yintercept=error), lty=2)
		)
	}
	return(p1)
}


#' Average survival
#'
#' @export
get_avesurv.gbm <- function(object, ...) {
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
#' Adopted from riskRegression predictRisk.gbm
#' @export
get_indivsurv.gbm <- function(object, newdata) {
	n.trees <- object$n.trees
	traindata <- object$modelData
	y <- model.extract(model.frame(terms(object), traindata), "response")
	xb.train <- predict(object, newdata = traindata, n.trees = n.trees)
	
	if (NCOL(y)==2) {
		times <- sort(unique(drop(y[,1])))
		H2 <- gbm::basehaz.gbm(t = y[,1], delta = y[,2]
			, f.x = xb.train
			, t.eval = times
			, cumulative = TRUE
		)
	} else {
		H2 <- gbm3basehaz.gbm(object)
		times <- H2$time
		H2 <- H2$hazard
	}
	if (!missing(newdata)) {
		xb.test <- predict(object, newdata = newdata , n.trees = n.trees)
	} else {
		xb.test <- xb.train
	}
	s <- matrix(0, nrow = length(xb.test), ncol = length(H2))
	
	for (i in 1:length(xb.test)) s[i,] <- exp(-exp(xb.test[i])*H2)
	out <- list(time = times, surv = s, chaz = -log(s))
	out$call <- match.call()
	class(out) <- "satsurv"
	out
}

#' Permutation variable importance method for gbm 
#'
#' @keywords internal

pvimp.gbm <- function(model, newdata, nrep = 20, parallelize = TRUE, nclusters = parallel::detectCores(), estimate=c("mean", "quantile"), probs=c(0.025, 0.5, 0.975), seed=NULL, ...){
	# Overall score
	if (is.null(seed) || !is.numeric(seed)) {
		seed <- 911
		set.seed(seed)
	}

	estimate <- match.arg(estimate)
	overall_c <- survconcord.gbm(model, newdata = newdata, stats = FALSE, ...)
	xvars <- all.vars(formula(delete.response(terms(model))))
	N <- NROW(newdata)
	if (parallelize) {
		## Setup parallel because serial takes a lot of time. Otherwise you can turn it off
		nn <- min(parallel::detectCores(), nclusters)
		if (nn < 2){
			foreach::registerDoSEQ()
		} else{
			cl <-  parallel::makeCluster(nn)
			doParallel::registerDoParallel(cl)
			on.exit(parallel::stopCluster(cl))
		}

		x <- NULL
		vi <- foreach(x = xvars, .export="survconcord.gbm") %dopar% {
			permute_df <- newdata[rep(seq(N), nrep), ]
			if (is.factor(permute_df[,x])) {
				permute_var <- as.vector(replicate(nrep, sample(newdata[,x], N, replace = FALSE)))
				permute_var <- factor(permute_var, levels = levels(permute_df[,x]))
			} else {
				permute_var <- as.vector(replicate(nrep, sample(newdata[,x], N, replace = FALSE)))
			}
			index <- rep(1:nrep, each = N)
			permute_df[, x] <- permute_var
			perm_c <- unlist(lapply(split(permute_df, index), function(d){
				survconcord.gbm(model, newdata = droplevels(d), stats = FALSE, ...)
			}))
			est <- (overall_c - perm_c)/overall_c
			out <- switch(estimate
				, mean={
					out2 <- mean(est)
					names(out2) <- x
					out2
				}
				, quantile={
					out2 <- cbind.data.frame(...xx=x, t(quantile(est, probs = probs, na.rm = TRUE)))
					colnames(out2) <- c("terms", "lower", "estimate", "upper")
					out2
				}
			)
			out
		}
	} else {
		vi <- sapply(xvars, function(x){
			permute_df <- newdata[rep(seq(N), nrep), ]
			if (is.factor(permute_df[,x])) {
				permute_var <- as.vector(replicate(nrep, sample(newdata[,x], N, replace = FALSE)))
				permute_var <- factor(permute_var, levels = levels(permute_df[,x]))
			} else {
				permute_var <- as.vector(replicate(nrep, sample(newdata[,x], N, replace = FALSE)))
			}
			index <- rep(1:nrep, each = N)
			permute_df[, x] <- permute_var
			perm_c <- unlist(lapply(split(permute_df, index), function(d){
				survconcord.gbm(model, newdata = droplevels(d), stats = FALSE, ...)
			}))
			est <- (overall_c - perm_c)/overall_c
			out <- switch(estimate
				, mean={
					out2 <- mean(est)
					out2
				}
				, quantile={
					out2 <- cbind.data.frame(...xx=x, t(quantile(est, probs = probs, na.rm = TRUE)))
					colnames(out2) <- c("terms", "lower", "estimate", "upper")
					out2
				}
			)
			out
		}, simplify=FALSE)
	}
	out <- switch(estimate
		, mean={
			unlist(vi)
		}
		, quantile={
			est <- do.call("rbind", vi)
			rownames(est) <- NULL
			est
		}
	)
	return(out)
}


#' Compute variable importance gbm
#'
#' @keywords internal

varimp.gbm <- function(object, type = c("coef", "perm", "model"), relative = TRUE, newdata, nrep = 20, parallelize = TRUE, nclusters = parallel::detectCores(), estimate=c("mean", "quantile"), probs=c(0.025, 0.5, 0.975), seed=NULL, ...){
	type <- match.arg(type)
	if (type=="perm"){
		out <- data.frame(Overall = get_pvimp(object, newdata, nrep, parallelize = parallelize, nclusters = nclusters, estimate=estimate, probs=probs, seed=seed, ...))
	} else {
		new_args <- list(...)
		new_args$object <- object
		new_args$plotit <- FALSE
		new_args$n.trees <- object$n.trees
		out1 <- do.call("summary", new_args)
		out <- data.frame(Overall = out1$rel.inf)
		rownames(out) <- out1$var
	}
	out$sign <- sign(out$Overall)
	out$Overall <- abs(out$Overall)
	if (relative){
		out$Overall <- out$Overall/sum(out$Overall, na.rm = TRUE)
	}
	return(out)
}

#' Extract predictions from gbm model
#'
#' Extract event probabilities from the fitted model.
#'
#' @aliases predictRisk
#'
#' @details
#' For survival outcome, the function predicts the risk, \eqn{1 - S(t|x)}, where \eqn{S(t|x)} is the survival chance of an individual characterized by \eqn{x}. riskRegression::predictRisk.gbm seems to have issues reconstructing the data.
#'
#' @importFrom riskRegression predictRisk
#' @export predictRisk
#' @export

predictRisk.gbm.satpred <- function(object, newdata, times, ...){
	p <- 1 - predictSurvProb.gbm.satpred(object, newdata, times)
	p
}

#' Compute predicted hazard for gbm3
#'
#' This code is borrowed from internal function agsurv from survival package.
#'
#' @param object fitted \code{\link[gbm3]{gbm}} object.
#' @return A list of S3 objects.
#' \item{n}{number of observations used in the fit.}
#' \item{events}{total number of events of interest in the fit.}
#' \item{time}{time points defined by the risk set.}
#' \item{n.risk}{the number of individuals at risk at time \code{t}.}
#' \item{n.event}{the number of events that occur at time \code{t}.}
#' \item{n.censor}{the number of subjects who exit the risk set, without an event, at time \code{t}.}
#' \item{surv}{a vector or a matrix of estimated survival function.}
#' \item{chaz, hazard}{a vector or a matrix of estimated cumulative hazard.}
#' @keywords internal

gbm3Hazard <- function(y, x = NULL, wt = rep(1, NROW(y)), risk=NULL, survtype=NULL, vartype=NULL){
	status <- y[, ncol(y)]
	dtime <- y[, ncol(y) - 1]
	death <- (status == 1)
	time <- sort(unique(dtime))
	nevent <- as.vector(rowsum(wt * death, dtime))
	ncens <- as.vector(rowsum(wt * (!death), dtime))
	wrisk <- wt * risk
	rcumsum <- function(x) rev(cumsum(rev(x)))
	nrisk <- rcumsum(rowsum(wrisk, dtime))
	irisk <- rcumsum(rowsum(wt, dtime))
	if (NCOL(y) != 2){
		delta <- min(diff(time))/2
		etime <- c(sort(unique(y[, 1])), max(y[, 1]) + delta)
		indx <- approx(etime, 1:length(etime), time, method = "constant", rule = 2, f = 1)$y
		esum <- rcumsum(rowsum(wrisk, y[, 1]))
		nrisk <- nrisk - c(esum, 0)[indx]
		irisk <- irisk - c(rcumsum(rowsum(wt, y[, 1])), 0)[indx]
	}
	haz <- nevent/nrisk
	result <- list(n = NROW(y), time = time, n.event = nevent
																, n.risk = irisk, n.censor = ncens, hazard = haz
																, chaz = cumsum(haz)
	)
	return(result)
}

#' Compute survival curve and cumulative hazard from a gbm3 through gbm model
#'
#' Compute the predicted survivor and cumulative hazard function for a penalized Cox proportional hazard model.
#'
#' @aliases gbm3survfit
#'
#' @export

gbm3survfit.gbm <- function(fit, newdata, ...) {
	y <- model.extract(model.frame(terms(fit), data=fit$modelData), "response")
	risk <- exp(predict(fit, newdata=fit$modelData, n.trees=fit$n.trees))
	afit <- gbm3Hazard(y = y, risk = risk)
	chaz <- afit$chaz
	surv.est <- exp(-chaz)
	if (!missing(newdata)) {
		lp <- predict(fit, newdata = newdata, n.trees=fit$n.trees)
		surv.est <- sapply(surv.est, function(x) x^exp(lp))
		chaz <- -log(surv.est)
	} 
	out <- list(n = afit$n
		, events = sum(afit$n.event)
		, time = afit$time
		, n.risk = afit$n.risk
		, n.event = afit$n.event
		, n.censor = afit$n.censor
		, surv = surv.est
		, chaz = chaz
	)
	out$call <- match.call()
	class(out) <- "satsurv"
	out
}

#' Compute survival curve and cumulative hazard from a gbm3 through gbm model
#'
#' Compute the baseline hazard
#'
#' @aliases gbm3basehaz
#'
#' @export

gbm3basehaz.gbm <- function(fit, centered = TRUE){
	sfit <- gbm3survfit.gbm(fit = fit)
	chaz <- sfit$chaz
	surv.est <- exp(-chaz)
	out <- list(time = sfit$time, hazard = chaz, surv = surv.est)
	class(out) <- "satsurv"
	out
}
