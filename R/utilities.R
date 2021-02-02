#' Implemented packages
#'
#' @export
rfsrc.satpred <- function(formula = NULL, train_df = NULL, test_df = NULL, ntree = 100, mtry = NULL, ...) {
	if (is.null(mtry)) {
		param <- expand.grid(ntree=ntree)
	} else {
		param <- expand.grid(ntree=ntree, mtry = mtry)
	}
	new_args <- list(...)
	rfsrc_args <- list(formula = formula, data = train_df, ntree = ntree, mtry = mtry, forest=TRUE)
	if (length(new_args)) rfsrc_args[names(new_args)] <- new_args
	error <- lapply(1:NROW(param), function(x){
		rfsrc_args$ntree <- param$ntree[[x]]
		rfsrc_args$mtry <- param$mtry[[x]]
		fit <- do.call("rfsrc.fast", rfsrc_args)
		if (is.null(test_df)) test_df <- train_df
		pred <- predict(fit, test_df)
		error_df <- data.frame(error = cverror(pred), ntree = fit$ntree, mtry=fit$mtry
			, nodesize = fit$nodesize, nodedepth = fit$nodedepth, splitrule = fit$splitrule
		)
		return(error_df)
	})
	error <- do.call("rbind", error)
#	class(error) <- "rfsrc"
	return(error)
}

#' Cross-validation error
#' @export
cverror.rfsrc <- function(x){
	return(x$err.rate[x$ntree])
}

#' @export
cverror <- function(x)UseMethod("cverror", x)

