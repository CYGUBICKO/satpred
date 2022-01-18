#' Model tuning
#'
#' @export

modtune <- function(formula = formula(data), data = sys.parent()
	, modfun, param_grid = NULL, nfolds = 10, foldids = NULL, ...) {
	
	Terms <- terms(formula, data = data)
	new_args <- list(...)
	mod_args <- list(formula = formula, data = data, modfun = modfun, param_grid = param_grid, nfolds = nfolds, foldids = foldids)
	if (length(new_args)) mod_args[names(new_args)] <- new_args
	default_args <- as.list(args(modfun))
	default_args <- default_args[!names(default_args) %in% c("train_df", "test_df", "...", "", names(mod_args))]
	if (length(default_args)) mod_args[names(default_args)] <- default_args
	cv <- do.call("modcv", mod_args)
	hyper <- colnames(cv)[!colnames(cv) %in% c("fold", "error")]
	form <- as.formula(paste0("error~", paste0(hyper, collapse = "+")))
	result <- aggregate(form, cv, FUN = function(x)mean(x, na.rm=TRUE))
	besTune <- getbesTune(result)
	mod_args <- list(formula=formula, data=data)
	mod_args[names(new_args)] <- new_args
	best_args <- as.list(besTune[, hyper])
	mod_args[names(best_args)] <- best_args
	default_args <- default_args[!names(default_args) %in% names(best_args)]
	if (length(default_args)) mod_args[names(default_args)] <- default_args 
	mod_args$finalmod <- TRUE
	mod_args <- mod_args[!names(mod_args) %in% c("parallelize", "nclusters", "nfolds", "foldids")]
	out <- list(result=result, besTune=besTune, modelfun = modfun, modelargs = mod_args)
	out$terms <- Terms
	out$call <- match.call()
	ff_call <- out$call[["modfun"]]
	if (ff_call=="gbm3.satpred") {
		ff_call <- "gbm.satpred"
	}
	class(out) <- c("satpred", ff_call)
	return(out)
}

#' Model fitting
#'
#' @export

modfit.satpred <- function(object, return_data = FALSE, ...) {
	new_args <- list(...)
	modfun <- object$modelfun
	mod_args <- object$modelargs
	if (length(new_args)) {
		if (any(names(new_args) %in% "early_stopping")) mod_args$early_stopping <- new_args$early_stopping
		new_args <- new_args[!names(new_args) %in% names(mod_args)]
		mod_args[names(new_args)] <- new_args
	}
	if (any(names(mod_args) %in% "num_nodes")) {
		temp_nodes <- as.numeric(strsplit(as.character(mod_args$num_nodes), "\\[|\\]$|\\,")[[1]])
		mod_args$num_nodes <- temp_nodes[!is.na(temp_nodes)]
	}
	finalModel <- do.call(modfun, mod_args)
	finalModel$terms <- object$terms
	if (return_data) {
		finalModel$modelData <- mod_args$data
	}
	finalModel$call <- match.call()
	if (inherits(finalModel, "gbm")) {
		class(finalModel) <- c("gbm.satpred", class(finalModel), "satpred")
	} else {
		class(finalModel) <- c(class(finalModel), "satpred")
	}
	return(finalModel)
}
