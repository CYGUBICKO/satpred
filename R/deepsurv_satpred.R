#' Implemented packages
#'
#' @export

deepsurv.satpred <- function(formula = NULL, train_df = NULL
	, test_df = NULL, num_nodes = list(h1 = c(32, 32))
	, param_grid = NULL, activation = "relu", dropout = 0.1
	, learning_rate = 1e-5, epochs = 10, batch_norm = TRUE
	, frac = 0, early_stopping = TRUE, finalmod = FALSE
	, seed = 9999, ...) {
	survivalmodels:::set_seed(seed)
	deepsurv_args <- list(formula=formula, data=train_df, activation=activation, batch_norm = batch_norm, frac = frac, early_stopping = early_stopping)
	if (is.null(param_grid)) {
		if (is.null(epochs)) {
			param <- expand.grid(learning_rate=learning_rate, stringsAsFactors=FALSE)
		} else if (is.null(dropout)) {
			param <- expand.grid(epochs=epochs, learning_rate=learning_rate, stringsAsFactors=FALSE)
		} else {
			param <- expand.grid(epochs=epochs, dropout=dropout, learning_rate=learning_rate, stringsAsFactors=FALSE)
		}
	} else {
		param <- param_grid
	}
	param_args <- as.list(param)
	deepsurv_args[names(param_args)] <- param_args
	new_args <- list(...)
	if (length(new_args)) deepsurv_args[names(new_args)] <- new_args

	if (!finalmod) {
		args_match <- match(colnames(param), names(deepsurv_args), nomatch = FALSE)
		param_match <- match(names(deepsurv_args), colnames(param), nomatch = FALSE)
		errors <- lapply(num_nodes, function(num_node){
			error <- lapply(1:NROW(param), function(x){
				deepsurv_args[args_match] <- param[x, param_match]
				deepsurv_args$num_nodes <- num_node
				fit <- do.call("deepsurv", deepsurv_args)
				if (is.null(test_df)) test_df <- train_df
				pred <- predict(fit, test_df, type = "risk")
				class(pred) <- c(class(pred), "deepsurv")
				deepsurv_args$num_nodes <- paste0("[", paste0(num_node, collapse=", "), "]")
				all_params <- c("num_nodes", "learning_rate", "dropout", "epochs")
				param_temp <- deepsurv_args[all_params]
				names(param_temp) <- all_params
				y  <- model.extract(model.frame(formula, data = test_df), "response")
				error_list <- list(param_temp, error = 1-cverror(pred, y))
				error_df <- as.data.frame(error_list)
				rownames(error_df) <- NULL
				return(error_df)
			})
			error <- do.call("rbind", error)
			return(error)
		})
		errors <- do.call("rbind", errors)
		return(errors)
	} else {
		fit <- do.call("deepsurv", deepsurv_args)
		return(fit)
	}
}


#' Cross-validation error
#'
#' Method for deepsurv
#'
#' @keywords internal

cverror.deepsurv <- function(x, y = NULL, ...){
	score <- survival::concordance(y~x, reverse=TRUE, ...)$concordance
	return(score)
}


#' Compute the concordance statistic for the deepsurv model
#'
#' The function computes the agreement between the observed response and the predictor.
#'
#' @keywords internal

survconcord.deepsurv <- function(object, newdata = NULL, stats = FALSE, ...) {
	if (is.null(newdata)) newdata <- object$modelData
	pred <- predict(object, newdata=newdata, type="risk")
	class(pred) <- c(class(pred), "deepsurv")
	y  <- model.extract(model.frame(terms(object), data = newdata), "response")
	concord <- cverror(pred, y)
	return(concord)
}


#' Cross-validation plots
#'
#' @import ggplot2
#' @export

plot.deepsurv.satpred <- function(x, ..., show_best = TRUE, lsize = 0.3, pshape = "O") {
	tune_df <- x$result
	tune_df$epochs <- factor(tune_df$epochs, labels=paste0("epochs: ", unique(tune_df$epochs)))
	tune_df$dropout <- factor(tune_df$dropout, labels=paste0("Dropout: ", unique(tune_df$dropout)))
	best_df <- x$besTune
	best_df$epochs <- factor(best_df$epochs, labels=paste0("epochs: ", unique(best_df$epochs)))
	best_df$dropout <- factor(best_df$dropout, labels=paste0("Dropout: ", unique(best_df$dropout)))
	p1 <- (ggplot(tune_df, aes(x = as.factor(learning_rate), y = error, group=as.factor(num_nodes), colour = as.factor(num_nodes)))
		+ geom_point(shape = pshape)
		+ geom_line(size = lsize)
		+ facet_grid(dropout~epochs)
		+ labs(x = "Learning rate (LR)", y = "Error (1 - C)", colour = "Hidden layer size")
	)
	if (show_best) {
		p1 <- (p1
			+ geom_point(data=best_df, aes(x = as.factor(learning_rate), y = error), colour="red", size=2)
			+ geom_hline(data=best_df, aes(yintercept=error), lty=2)
		)
	}
	return(p1)
}


#' Average survival
#'
#' @export
get_avesurv.deepsurv <- function(object, ...) {
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
#' Adopted from riskRegression predictRisk.deepsurv
#' @export
get_indivsurv.deepsurv <- function(object, newdata) {
	if (missing(newdata)) {
		newdata <- object$modelData
	}
	surv <- predict(object, newdata = newdata, type = "survival")
	times <- as.numeric(colnames(surv))
	times <- times[1:(length(times)-1)]
	ss <- surv[,1:length(times), drop=FALSE]
	out <- list(time = times, surv = ss, chaz = -log(ss))
	out$call <- match.call()
	class(out) <- "satsurv"
	out
}

#' Permutation variable importance method for deepsurv 
#'
#' @keywords internal

pvimp.deepsurv <- function(model, newdata, nrep = 20, parallelize = TRUE, nclusters = parallel::detectCores(), model_matrix = FALSE, rhs_formula = formula(model)[c(1,3)], scale = TRUE, center = TRUE, ...){
	# Overall score
	Terms <- terms(model)
	xvars <- all.vars(formula(delete.response(Terms)))
	ynames <- deparse(formula(Terms)[[2]])
	eventlab <- trimws(gsub(".*\\,|\\)", "", ynames))
	timelab <- gsub(".*\\(|\\,.*", "", ynames)
	form <- rhs_formula 
	N <- NROW(newdata)
	if (model_matrix) {
		xvars <- all.vars(delete.response(form))
		temp_df <- model.matrix(form, newdata)[,-1, drop=FALSE]
		if (center!=TRUE || center!=FALSE) {
			temp_df <- temp_df[, names(center), drop=FALSE]
		} 
		temp_df <- as.data.frame(scale(temp_df, center=center, scale=scale))
		temp_df[[eventlab]] <- newdata[[eventlab]]
		temp_df[[timelab]] <- newdata[[timelab]]
	} else {
		temp_df <- newdata
	}
	overall_c <- survconcord.deepsurv(model, newdata = temp_df, stats = FALSE, ...)
	if (parallelize) {
		## Setup parallel because serial takes a lot of time. Otherwise you can turn it off
		## FIXME: Parallelize breaks, setting nclusters to 1
		nclusters <- 1
		message("Ignores parallelize, running in serial...\n")
		nn <- min(parallel::detectCores(), nclusters)
		if (nn < 2){
			foreach::registerDoSEQ()
		} else{
			cl <-  parallel::makeCluster(nn)
			doParallel::registerDoParallel(cl)
			on.exit(parallel::stopCluster(cl))
		}

		x <- NULL
		vi <- foreach(x = xvars, .packages="survivalmodels", .export=c("deepsurv", "survconcord.deepsurv")) %dopar% {
			permute_df <- newdata[rep(seq(N), nrep), ]
			events <- permute_df[[eventlab]] 
			times <-  permute_df[[timelab]]
			if (is.factor(permute_df[,x])) {
				permute_var <- as.vector(replicate(nrep, sample(newdata[,x], N, replace = FALSE)))
				permute_var <- factor(permute_var, levels = levels(permute_df[,x]))
			} else {
				permute_var <- as.vector(replicate(nrep, sample(newdata[,x], N, replace = FALSE)))
			}
			index <- rep(1:nrep, each = N)
			permute_df[, x] <- permute_var
			if (model_matrix) {
				permute_df <- model.matrix(form, permute_df)[,-1, drop=FALSE]
				if (center!=TRUE || center!=FALSE) {
					permute_df <- permute_df[, names(center), drop=FALSE]
				} 
				permute_df <- as.data.frame(scale(permute_df, center=center, scale=scale))
				permute_df[[eventlab]] <- events 
				permute_df[[timelab]] <- times 
			}
			perm_c <- unlist(lapply(split(permute_df, index), function(d){
				out <- survconcord.deepsurv(model, newdata = droplevels(d), stats = FALSE, ...)
				return(out)
			}))
			est <- mean((overall_c - perm_c)/overall_c)
			names(est) <- x
			est
		}
	} else {
		vi <- sapply(xvars, function(x){
			permute_df <- newdata[rep(seq(N), nrep), ]
			events <- permute_df[[eventlab]] 
			times <-  permute_df[[timelab]]
			if (is.factor(permute_df[,x])) {
				permute_var <- as.vector(replicate(nrep, sample(newdata[,x], N, replace = FALSE)))
				permute_var <- factor(permute_var, levels = levels(permute_df[,x]))
			} else {
				permute_var <- as.vector(replicate(nrep, sample(newdata[,x], N, replace = FALSE)))
			}
			index <- rep(1:nrep, each = N)
			permute_df[, x] <- permute_var
			if (model_matrix) {
				permute_df <- model.matrix(form, permute_df)[,-1, drop=FALSE]
				if (center!=TRUE || center!=FALSE) {
					permute_df <- permute_df[, names(center), drop=FALSE]
				} 
				permute_df <- as.data.frame(scale(permute_df, center=center, scale=scale))
				permute_df[[eventlab]] <- events 
				permute_df[[timelab]] <- times 
			}
			perm_c <- unlist(lapply(split(permute_df, index), function(d){
				survconcord.deepsurv(model, newdata = droplevels(d), stats = FALSE, ...)
			}))
			mean((overall_c - perm_c)/overall_c)
		})
	}
	return(unlist(vi))
}

#' Compute variable importance deepsurv
#'
#' @keywords internal

varimp.deepsurv <- function(object, type = c("coef", "perm", "model"), relative = TRUE, newdata, nrep = 20, parallelize = TRUE, nclusters = parallel::detectCores(), ...){
	type <- match.arg(type)
	if (type!="perm")warning(paste0("type = ", type, " not implemented yet, using 'perm'"))
	if (type=="perm"){
		out <- data.frame(Overall = get_pvimp(object, newdata, nrep, parallelize = parallelize, nclusters = nclusters, ...))
	} else {
		out <- data.frame(Overall = get_pvimp(object, newdata, nrep, parallelize = parallelize, nclusters = nclusters, ...))
	}
	out$sign <- sign(out$Overall)
	out$Overall <- abs(out$Overall)
	if (relative){
		out$Overall <- out$Overall/sum(out$Overall, na.rm = TRUE)
	}
	return(out)
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

predictSurvProb.deepsurv <- function(object, newdata, times, ...){
	N <- NROW(newdata)
	sfit <- get_indivsurv(object, newdata = newdata)
	S <- sfit$surv
	Time <- sfit$time
	S <- matrix(S, ncol=length(Time))
	if(N == 1) S <- matrix(S, nrow = 1)
	p <-  cbind(1, S)[, 1 + prodlim::sindex(Time, times),drop = FALSE]
	if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
		stop("Prediction failed")
	p
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

predictRisk.deepsurv <- function(object, newdata, times, ...){
	p <- 1 - predictSurvProb.deepsurv(object, newdata, times)
	p
}


