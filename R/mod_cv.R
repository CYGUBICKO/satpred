#' Cross-validation
#'
#' @export

modcv <- function(formula = formula(data), data = sys.parent(), modfun, nfolds = 10, foldids = NULL, ...) {
	
	if (is.null(foldids)){
	N <- NROW(data)
		foldids <- sample(rep(seq(nfolds), length.out = N))
	} else {
		nfolds = max(foldids)
	}
	if (nfolds < 3)stop("Number of folds should be at least 3: nfolds = 10 recommended")
	new_args <- list(...)
	mod_args <- list(formula = formula, data = data)
	if (length(new_args)) mod_args[names(new_args)] <- new_args
	out <- list()
	for (f in 1:nfolds) {
		index <- which(foldids==f)
		train_df <- data[-index, ]
		heldout_df <- data[index, ]
		mod_args$train_df <- train_df
		mod_args$test_df <- heldout_df
		est <- do.call(modfun, mod_args)
		est[, "fold"] <- paste0("fold", f)
		out[[paste0("fold", f)]] <- est
	}
	out <- do.call("rbind", out)
	rownames(out) <- NULL
	out <- out[, union("fold", colnames(out))]
	return(out)
}
