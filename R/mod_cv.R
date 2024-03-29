#' Cross-validation
#'
#' @export
#' @import parallel
#' @import foreach
#' @import doParallel

modcv <- function(formula = formula(data), data = sys.parent(), modfun, nfolds = 10, foldids = NULL, parallelize = FALSE, nclusters = parallel::detectCores(), cvseed=911, ...) {
	
	set.seed(cvseed)
	if (is.null(foldids)){
	N <- NROW(data)
		foldids <- sample(rep(seq(nfolds), length.out = N))
	} else {
		nfolds = max(foldids)
	}
	## 2022 Jan 25 (Tue): For compatibility with gbm
#	if (nfolds < 3)stop("Number of folds should be at least 3: nfolds = 10 recommended")
	new_args <- list(...)
	mod_args <- list(formula = formula)
	if (length(new_args)) mod_args[names(new_args)] <- new_args

	# Perform CV
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

		f <- NULL
		out <- foreach(f = 1:nfolds, .combine = rbind, .multicombine = TRUE) %dopar% {
			set.seed(cvseed)
			index <- which(foldids==f)
			if (nfolds==1) {
				train_df <- data
				heldout_df <- NULL
			} else {
				train_df <- data[-index, ]
				heldout_df <- data[index, ]
			}
			mod_args$train_df <- train_df
			mod_args$test_df <- heldout_df
			est <- do.call(modfun, mod_args)
			est[, "fold"] <- paste0("fold", f)
			est
		}
	} else if (!parallelize || nfolds==1) {
		out <- list()
		for (f in 1:nfolds) {
			index <- which(foldids==f)
			if (nfolds==1) {
				train_df <- data
				heldout_df <- NULL
			} else {
				train_df <- data[-index, ]
				heldout_df <- data[index, ]
			}
			mod_args$train_df <- train_df
			mod_args$test_df <- heldout_df
			est <- do.call(modfun, mod_args)
			est[, "fold"] <- paste0("fold", f)
			out[[paste0("fold", f)]] <- est
		}
		out <- do.call("rbind", out)
		rownames(out) <- NULL
	}
	out <- out[, union("fold", colnames(out))]
	return(out)
}
