#' Implemented packages
#'
#' @export
rfsrc.satpred <- function(formula = NULL, train_df = NULL, test_df = NULL, param_grid = NULL
	, ntree = 1000, mtry = NULL, nodesize = NULL, splitrule = "logrank", finalmod = FALSE, ...) {
	
	rfsrc_args <- list(formula=formula, data=train_df, forest=TRUE)
	if (is.null(param_grid)) {
		if (is.null(mtry)) {
			param <- expand.grid(ntree=ntree, splitrule=splitrule)
		} else if (is.null(nodesize)) {
			param <- expand.grid(ntree=ntree, mtry=mtry, splitrule=splitrule)
		} else {
			param <- expand.grid(ntree=ntree, mtry=mtry, nodesize=nodesize, splitrule=splitrule)
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
			fit <- do.call("rfsrc.fast", rfsrc_args)
			if (is.null(test_df)) test_df <- train_df
			pred <- predict(fit, test_df)
			all_params <- names(param_args)
			all_params <- union(c("mtry", "ntree", "nodesize", "splitrule"), all_params)
			param_temp <- fit[all_params]
			names(param_temp) <- all_params
			error_list <- list(param_temp, error = cverror(pred))
			error_df <- as.data.frame(error_list)
			return(error_df)
		})
		error <- do.call("rbind", error)
		return(error)
	} else {
		fit <- do.call("rfsrc", rfsrc_args)
		return(fit)
	}
}

#' Cross-validation error
#' @export
cverror.rfsrc <- function(x){
	return(x$err.rate[x$ntree])
}

#' @export
cverror <- function(x)UseMethod("cverror", x)

#' Get best tune
#'
#' @keywords internal
getbesTune <- function(x) {
	x <- x[which.min(x$error),]
	return(x)
}

