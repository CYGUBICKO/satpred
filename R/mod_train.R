#' Model tuning
#'
#' @export

modtune <- function(formula = formula(data), data = sys.parent()
	, modfun, param_grid = NULL, nfolds = 10, foldids = NULL, ...) {
	
	Terms <- terms(formula, data = data)
	new_args <- list(...)
	mod_args <- list(formula = formula, data = data, modfun = modfun, param_grid = param_grid, nfolds = nfolds, foldids = foldids)
	if (length(new_args)) mod_args[names(new_args)] <- new_args
	cv <- do.call("modcv", mod_args)
	hyper <- colnames(cv)[!colnames(cv) %in% c("fold", "error")]
	form <- as.formula(paste0("error~", paste0(hyper, collapse = "+")))
	result <- aggregate(form, cv, FUN = function(x)mean(x, na.rm=TRUE))
	besTune <- getbesTune(result)
	mod_args <- list(formula=formula, data=data)
	mod_args$finalmod <- TRUE
	mod_args[names(new_args)] <- new_args
	best_args <- as.list(besTune[, hyper])
	mod_args[names(best_args)] <- best_args
	out <- list(result=result, besTune=besTune, modelfun = modfun, modelargs = mod_args)
	out$terms <- Terms
	out$call <- match.call()
	class(out) <- c("satpred", out$call[["modfun"]])
	return(out)
}

#' Model fitting
#'
#' @export

modfit.satpred <- function(object, return_data = TRUE, ...) {
	new_args <- list(...)
	modfun <- object$modelfun
	mod_args <- object$modelargs
	if (length(new_args)) {
		new_args <- new_args[!names(new_args) %in% names(mod_args)]
		mod_args[names(new_args)] <- new_args
	}
	finalModel <- do.call(modfun, mod_args)
	finalModel$terms <- object$terms
	if (return_data) {
		finalModel$modelData <- mod_args$data
	}
	finalModel$call <- match.call()
	if (inherits(finalModel, "gbm")){
		class(finalModel) <- c("gbm.satpred", class(finalModel), "satpred")
	} else {
		class(finalModel) <- c(class(finalModel), "satpred")
	}
	return(finalModel)
}
