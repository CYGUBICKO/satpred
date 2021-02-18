#' Implemented packages
#'
#' @export
gbm.satpred <- function(formula = NULL, train_df = NULL, test_df = NULL, distribution = "coxph", param_grid = NULL, n.trees = 1000, interaction.depth = 1, n.minobsinnode = 10, shrinkage = 0.1, finalmod = FALSE, ...) {
	
	gbm_args <- list(formula=formula, data=train_df, distribution=distribution)
	if (is.null(param_grid)) {
		if (is.null(shrinkage)) {
			param <- expand.grid(n.trees=n.trees, n.minobsinnode=n.minobsinnode)
		} else if (is.null(interaction.depth)) {
			param <- expand.grid(n.trees=n.trees, shrinkage=shrinkage, n.minobsinnode=n.minobsinnode)
		} else {
			param <- expand.grid(n.trees=n.trees, shrinkage=shrinkage, interaction.depth=interaction.depth, n.minobsinnode=n.minobsinnode)
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
			fit <- do.call("gbm", gbm_args)
			if (is.null(test_df)) test_df <- train_df
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
		fit <- do.call("gbm", gbm_args)
		return(fit)
	}
}


#' Cross-validation error
#'
#' Method for gbm
#'
#' @keywords internal

cverror.gbm <- function(x, y = NULL, ...){
	score <- Hmisc::rcorr.cens(-x, y)[[1]]
	return(score)
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

