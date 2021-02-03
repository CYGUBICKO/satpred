#' Model training
#'
#' @export

modtrain <- function(formula = formula(data), data = sys.parent()
	, modfun, param_grid = NULL, nfolds = 10, foldids = NULL, finalmod = TRUE, ...) {
	new_args <- list(...)
	mod_args <- list(formula = formula, data = data, modfun = modfun, param_grid = param_grid, nfolds = nfolds, foldids = foldids)
	if (length(new_args)) mod_args[names(new_args)] <- new_args
	cv <- do.call("modcv", mod_args)
	hyper <- colnames(cv)[!colnames(cv) %in% c("fold", "error")]
	form <- as.formula(paste0("error~", paste0(hyper, collapse = "+")))
	result <- aggregate(form, cv, FUN = function(x)mean(x, na.rm=TRUE))
	besTune <- getbesTune(result)
	finalModel <- NULL
	if (finalmod) {
		mod_args <- list(formula=formula, data=data)
		mod_args$finalmod <- TRUE
		mod_args[names(new_args)] <- new_args
		best_args <- as.list(besTune[, hyper])
		mod_args[names(best_args)] <- best_args
		finalModel <- do.call(modfun, mod_args) 
	}
	out <- list(result=result, besTune=besTune, finalModel=finalModel)
	out$call <- call()
	class(out) <- "satpred"
	return(out)
}
