#' Implemented packages
#'
#' @export
rfsrc.satpred <- function(formula = NULL, train_df = NULL, test_df = NULL, param_grid = NULL
	, ntree = 1000, mtry = NULL, nodesize = NULL, splitrule = "logrank", forest = FALSE
	, finalmod = FALSE, ...) {
	
	rfsrc_args <- list(formula=formula, data=train_df, forest=forest)
	if (is.null(param_grid)) {
		if (is.null(mtry)) {
			param <- expand.grid(ntree=ntree, splitrule=splitrule, stringsAsFactors=FALSE)
		} else if (is.null(nodesize)) {
			param <- expand.grid(ntree=ntree, mtry=mtry, splitrule=splitrule, stringsAsFactors=FALSE)
		} else {
			param <- expand.grid(ntree=ntree, mtry=mtry, nodesize=nodesize, splitrule=splitrule, stringsAsFactors=FALSE)
		}
	} else {
		param <- param_grid
	}
	param_args <- as.list(param)
#	names(param_args) <- param_args
	rfsrc_args[names(param_args)] <- param_args
	new_args <- list(...)
	if (length(new_args)) rfsrc_args[names(new_args)] <- new_args

	if (!finalmod) {
		args_match <- match(colnames(param), names(rfsrc_args), nomatch = FALSE)
		param_match <- match(names(rfsrc_args), colnames(param), nomatch = FALSE)
		error <- lapply(1:NROW(param), function(x){
			rfsrc_args[args_match] <- param[x, param_match]
			check_ndepth <- rfsrc_args$nodedepth
			if (!is.null(check_ndepth)) {
				if (check_ndepth==0) {
					rfsrc_args['nodedepth'] <- list(NULL)
				}
			}
			fit <- do.call("rfsrc.fast", rfsrc_args)
			if (is.null(test_df)) test_df <- train_df
			pred <- predict(fit, test_df)
			all_params <- names(param_args)
			all_params <- union(c("mtry", "ntree", "nodesize", "splitrule"), all_params)
			param_temp <- fit[all_params]
			names(param_temp) <- all_params
			if (check_ndepth==0) {
					param_temp$nodedepth <- 0
			}
			error_list <- list(param_temp, error = cverror(pred))
			error_df <- as.data.frame(error_list)
			return(error_df)
		})
		error <- do.call("rbind", error)
		return(error)
	} else {
		check_ndepth <- rfsrc_args$nodedepth
		if (!is.null(check_ndepth)) {
			if (check_ndepth==0) {
				rfsrc_args['nodedepth'] <- list(NULL)
			}
		}
		fit <- do.call("rfsrc", rfsrc_args)
		return(fit)
	}
}

#' Cross-validation error
#'
#' @keywords internal

cverror.rfsrc <- function(x, y = NULL, ...){
	return(x$err.rate[x$ntree])
}

#' Tidy model objects
#'
#' @keywords internal

modtidy.rfsrc <- function(x) {
	x$time <- x$time.interest
	x$surv <- x$survival
	x$chaz <- x$chf
	return(x)
}

#' Tidy prediction objects
#'
#' @keywords internal

predtidy.rfsrc <- function(x) {
	x$time <- x$time.interest
	x$surv <- x$survival
	x$chaz <- x$chf
	return(x)
}


#' Compute the concordance statistic for the randomRoresrSRC model
#'
#' The function computes the agreement between the observed response and the predictor.
#'
#' @keywords internal

survconcord.rfsrc <- function(object, newdata = NULL, stats = FALSE, ...) {
	if (is.null(newdata)) {
		pred <- predict(object, ...)
	} else {
		pred <- predict(object, newdata, ...)
	}
	concord <- 1 - pred$err.rate[pred$ntree]
	return(concord)
}

#' Compute the concordance statistic for the default model (coxph, glmnetsurv)
#'
#' The function computes the agreement between the observed response and the predictor.
#'
#' @keywords internal

survconcord.default <- function(object, newdata = NULL, stats = FALSE, ...){
	if (is.null(newdata)) {
		risk <- predict(object, type = "risk", ...)
		y <- object$y
	} else {
		risk <- predict(object, newdata = newdata, type = "risk", ...)
		y <- model.extract(model.frame(object$terms, data = newdata), "response")
	}
	conindex <- survival::survConcordance(y ~ risk)
	if (!stats){
		conindex <- conindex$concordance
	}
	return(conindex)
}


#' Permutation variable importance method for rfsrc 
#'
#' @keywords internal

pvimp.rfsrc <- function(model, newdata, nrep = 20, parallelize = TRUE, nclusters = parallel::detectCores(), ...){
	# Overall score
	overall_c <- survconcord.rfsrc(model, newdata = newdata, stats = FALSE, ...)
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
		vi <- foreach(x = xvars, .export="survconcord.rfsrc") %dopar% {
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
			#	survconcord.rfsrc(model, newdata = droplevels(d), stats = FALSE, ...)
				survconcord.rfsrc(model, newdata = d, stats = FALSE, ...)
			}))
			est <- mean((overall_c - perm_c)/overall_c)
			names(est) <- x
			est
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
			#	survconcord.rfsrc(model, newdata = droplevels(d), stats = FALSE, ...)
				survconcord.rfsrc(model, newdata = d, stats = FALSE, ...)
			}))
			mean((overall_c - perm_c)/overall_c)
		})
	}
	return(unlist(vi))
}

#' Permutation variable importance method for default models 
#'
#' @keywords internal

pvimp.default <- function(model, newdata, nrep = 20, parallelize = TRUE, nclusters = parallel::detectCores(), ...){
	# Overall score
	overall_c <- survconcord.default(model, newdata = newdata, stats = FALSE, ...)
	Terms <- terms(model)
	yvars <- formula(Terms)[[2]]
	y <- with(newdata, eval(yvars))
	xvars <- all.vars(formula(delete.response(Terms)))
	N <- NROW(newdata)
	newdata <- newdata[, xvars, drop = FALSE]
	
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
		vi <- foreach(x = xvars, .export = "survConcordance") %dopar% {
			permute_df <- newdata[rep(seq(N), nrep), ]
			if (is.factor(permute_df[,x])) {
				permute_var <- as.vector(replicate(nrep, sample(newdata[,x], N, replace = FALSE)))
				permute_var <- factor(permute_var, levels = levels(permute_df[,x]))
			} else {
				permute_var <- as.vector(replicate(nrep, sample(newdata[,x], N, replace = FALSE)))
			}
			index <- rep(1:nrep, each = N)
			permute_df[, x] <- permute_var
			risk <- predict(model, newdata = permute_df, type = "risk")
			perm_c <- tapply(risk, index, function(r){
				survConcordance(y~r)$concordance
			})
			est <- mean((overall_c - perm_c)/overall_c)
			names(est) <- x
			est
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
			risk <- predict(model, newdata = permute_df, type = "risk")
			perm_c <- tapply(risk, index, function(r){
				survConcordance(y~r)$concordance
			})
			mean((overall_c - perm_c)/overall_c)
		})
	}
	return(unlist(vi))
}

#' Compute variable importance random forest
#'
#' @keywords internal

varimp.rfsrc <- function(object, type = c("coef", "perm", "model"), relative = TRUE, newdata, nrep = 20, parallelize = TRUE, nclusters = parallel::detectCores(), ...){
	type <- match.arg(type)
	if (type=="perm"){
		out <- data.frame(Overall = get_pvimp(object, newdata, nrep, parallelize = parallelize, nclusters = nclusters, ...))
	} else {
		new_args <- list(...)
		new_args$object <- object
		out <- do.call("vimp", new_args)$importance
		out <- data.frame(Overall = out)
	}
	out$sign <- sign(out$Overall)
	out$Overall <- abs(out$Overall)
	if (relative){
		out$Overall <- out$Overall/sum(out$Overall, na.rm = TRUE)
	}
	return(out)
}

#' Compute variable importance for default models
#'
#' @keywords internal

varimp.default <- function(object, type = c("coef", "perm", "model"), relative = TRUE, newdata, nrep = 20, parallelize = TRUE, nclusters = parallel::detectCores(), ...){
	type <- match.arg(type)
	if (type=="perm"){
		out <- data.frame(Overall = get_pvimp(object, newdata, nrep, parallelize = parallelize, nclusters = nclusters, ...))
	} else {
		if (inherits(object, "glmnetsurv")){
			s <- object$s
			beta <- predict(object$fit, s = s, type = "coef")
			if(is.list(beta)) {
				out <- do.call("cbind", lapply(beta, function(x) x[,1]))
				out <- as.data.frame(out)
			} else out <- data.frame(Overall = beta[,1])	
		} else {
			beta <- coef(object)
			out <- data.frame(Overall = beta)
		}
		out <- out[rownames(out) != "(Intercept)",,drop = FALSE]
	}
	out$sign <- sign(out$Overall)
	out$Overall <- abs(out$Overall)
	if (relative){
		out$Overall <- out$Overall/sum(out$Overall, na.rm = TRUE)
	}
	return(out)
}

#' Get best tune
#'
#' @keywords internal
getbesTune <- function(x) {
	x <- x[which.min(x$error),]
	return(x)
}

