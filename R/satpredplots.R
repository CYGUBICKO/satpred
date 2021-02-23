#' Cross-validation plots
#'
#' @import ggplot2
#' @export

plot.rfsrc.satpred <- function(x, ..., show_best = TRUE, lsize = 0.3, pshape = "O") {
	tune_df <- x$result
	best_df <- x$besTune
	tune_df$nodesize <- factor(tune_df$nodesize, labels=paste0("nodesize: ", unique(tune_df$nodesize)))
	best_df <- x$besTune
	best_df$nodesize <- factor(best_df$nodesize, labels=paste0("nodesize: ", unique(best_df$nodesize)))
	p1 <- (ggplot(tune_df, aes(x = as.factor(mtry), y = error, group=as.factor(ntree), colour = as.factor(ntree)))
		+ geom_point(shape = pshape)
		+ geom_line(size = lsize)
		+ facet_grid(splitrule~nodesize)
		+ labs(x = "# randomly selected predictors", y = "Error (1 - C)", colour = "# trees")
	)
	if (show_best) {
		p1 <- (p1
			+ geom_point(data=best_df, aes(x = as.factor(mtry), y = error), colour="red", size=2)
			+ geom_hline(data=best_df, aes(yintercept=error), lty=2)
		)
	}
	return(p1)
}

#' Plot survival and cumulative hazard curves
#'
#' Plot estimated survival and cumulative  hazard curves for \code{satpred} model.
#'
#' @import ggplot2
#' @export

plot.satsurv <- function(x, ..., type = c("surv", "cumhaz"), lsize = 0.3, lcol = "grey", ltype = 2, compare = FALSE, compare_lab = "satpred") {
	type <- match.arg(type)

	cumhaz <- x$chaz
	surv <- x$surv
	avesurv <- surv
	avechaz <- cumhaz
	if (NCOL(surv)>1) {
		avesurv <- as.vector(sapply(data.frame(surv), mean))
		avechaz <- as.vector(sapply(data.frame(cumhaz), mean))
		cumhaz <- t(cumhaz)
		surv <- t(surv)
	}
	time <- x$time
	if (min(time)>0) {
		if (NCOL(surv) > 1) {
			cumhaz <- rbind(0, cumhaz)
			surv <- rbind(1, surv)
		} else {
			cumhaz <- c(0, cumhaz)
			surv <- c(1, surv)
		}
		avechaz <- c(0, avechaz)
		avesurv <- c(1, avesurv)
		time <- c(0, time)
	}
	if (NCOL(surv) > 1){
		nindivs <- NCOL(surv)
		individ <- as.factor(rep(1:nindivs, each = length(time)))
		surv <- as.vector(surv)
		cumhaz <- as.vector(cumhaz)
		time <- rep(time, nindivs)
		plot_df <- data.frame(id = individ, time = time, surv = surv, cumhaz = cumhaz, avechaz = avechaz, avesurv = avesurv)
	} else {
		plot_df <- data.frame(id = 1, time = time, surv = surv, cumhaz = cumhaz, avechaz = avechaz, avesurv = avesurv)
	}
	id <- NULL
	p0 <- ggplot(plot_df, aes(x = time, group = id)) + labs(x = "Time")
	if (type == "surv"){
		p1 <- (p0 
			+ geom_step(aes(y = avesurv), colour = "red")
			+ geom_step(aes(y = surv), colour = lcol, size = lsize, lty = ltype) + labs(y = "Survival prob.")
		)
		if (compare){
			p1 <- p0 + geom_step(aes(y = surv, colour = compare_lab), size = lsize) + labs(y = "Survival prob.")
		}
	} else {
		p1 <- (p0 
			+ geom_step(aes(y = avechaz), colour = "red")
			+ geom_step(aes(y = cumhaz), colour = lcol, size = lsize, lty = ltype) + labs(y = "Cumulative hazard")
		)
		if (compare){
			p1 <- p0 + geom_step(aes(y = cumhaz, colour = compare_lab), size = lsize) + labs(y = "Cumulative hazard")
		}
	}
	return(p1)
}



#' Prediction performance
#'
#' Plots predictive performance comparison of survival models. It uses risk scoring from \code{\link[riskRegression]{Score}}. 
#'
#' @details
#' Implements plot method for \code{\link[riskRegression]{Score}} for time-dependent Brier score, AUC and ROC. However, currently, no support for time-dependent covariate models.
#'
#' @param x \code{\link[riskRegression]{Score}} object. See examples.
#' @param ... for future implementations.
#' @param type metric to return. Choices are \code{"roc", "auc", "brier"}.
#' @param pos spacing between the lines.
#'
#' @return a \code{\link[ggplot2]{ggplot}} object.
#'
#' @import ggplot2
#' @export

plot.Score <- function(x, ..., type = c("roc", "auc", "brier"), pos = 0.3){
	if (!inherits(x, "Score"))
		stop("Object should be score. See ?riskRegression::Score")
	type <- match.arg(type)
	if (type == "roc"){
		df <- x$ROC$plotframe
		df$times <- as.factor(df$times)
		FPR <- TPR <- model <- AUC <- lower <- upper <- Brier <- NULL
		model_cols <- unique(df$model)
		p1 <- (ggplot(df, aes(x = FPR, y = TPR, color = model))
			+ geom_line(size = 1)
			+ geom_abline(size = 1, colour = "grey")
			+ facet_wrap(~times)
			+ labs(x = "1-Specificity", y = "Sensitivity", colour = "Time")
#			+ scale_colour_viridis_d(option = "inferno")
			+ scale_color_manual(breaks = model_cols
				, values = rainbow(n = length(model_cols))
			)
			+ theme(legend.position = "right")
		)
	} else if (type == "auc"){
		df <- x$AUC$score
		df$times <- as.factor(df$times)
		p1 <- (ggplot(df, aes(x = times, y = AUC, group = model, colour = model))
			+ geom_point(position = position_dodge(pos))
			+ geom_pointrange(aes(ymin = lower, ymax = upper, colour = model), position = position_dodge(pos))
#			+ scale_colour_viridis_d(option = "inferno")
			+ scale_color_manual(breaks = model_cols
				, values = rainbow(n = length(model_cols))
			)
			+ labs(x = "Time", y = "AUC", colour = "Model")
			+ theme(legend.position = "right")
		)
	} else {
		df <- x$Brier$score
		df$times <- as.factor(df$times)
		p1 <- (ggplot(df, aes(x = times, y = Brier, group = model, colour = model))
			+ geom_point(position = position_dodge(pos))
			+ geom_pointrange(aes(ymin = lower, ymax = upper, colour = model), position = position_dodge(pos))
#			+ scale_colour_viridis_d(option = "inferno")
			+ scale_color_manual(breaks = model_cols
				, values = rainbow(n = length(model_cols))
			)
			+ labs(x = "Time", y = "Brier", colour = "Model")
			+ theme(legend.position = "right")
		)
	}
	return(p1)
}


#' Plotting prediction error curves
#'
#' @import ggplot2
#' @export

plotpec <- function(x, ..., lsize = 0.3, ltype = 2, xlab = "Time", ylab = "Prediction error") {
	if (!inherits(x, "pec")) stop("Needs a pec object. See ?pec::pec")
	df <- do.call("data.frame", list(x$AppErr, times=x$time))
	vnames <- colnames(df)[!colnames(df) %in% "times"]
	df <- reshape(df, timevar = "model", v.names = "score"
		, varying = vnames, times = vnames, direction = "long"
	)
	rownames(df) <- NULL
	model_cols <- unique(df$model)
	p1 <- (ggplot(df, aes(x = times, y = score, colour = model))
		+ geom_line()
		+ scale_color_manual(breaks = model_cols
			, values = rainbow(n = length(model_cols))
		)
		+ labs(x = xlab, y = ylab, colour = "Model")
		+ theme(legend.position = "right")
	)
	return(p1)
}

#' Generic method for plotting variable importance of various models
#'
#' @import ggplot2
#' @export

plot.varimp <- function(x, ..., pos = 0.3, drop_zero = TRUE){
	x$sign <- ifelse(x$sign==1, "+", ifelse(x$sign==-1, "-", "0"))
	x <- x[order(x$Overall), ]
	if (drop_zero){
		x <- x[x$Overall!=0, ]
		x <- droplevels(x)
	}
	Overall <- NULL
	nmods <- unique(x$model)
	pos <- position_dodge(width = pos)
	if (length(nmods)==1) {
		p0 <- ggplot(x, aes(x = reorder(terms, Overall), y = Overall)) 
	} else {
		p0 <- (ggplot(x, aes(x = reorder(terms, Overall), y = Overall, colour = model))
			+ scale_color_manual(breaks = nmods
				, values = rainbow(n = length(nmods))
			)
			+ labs(colour = "Model")
		)
	}
	p1 <- (p0
		+ geom_point(aes(shape=sign), position = pos)
		+ geom_linerange(aes(ymin = 0, ymax = Overall, lty = sign), position = pos)
		+ scale_shape_manual(name = "Sign", values=c(1,16, 15))
		+ labs(x = "", y = "Importance", linetype = "Sign")
		+ coord_flip(clip = "off", expand = TRUE)
		+ theme_minimal()	
	)
	return(p1)
}

#' Customized theme for satpred plots
#'
#' Sets a theme for satpred and other ggplot objects
#'
#' @import ggplot2
#' @export

satpredtheme <- function(){
   theme_set(theme_bw() +
      theme(panel.spacing = grid::unit(0,"lines")
      	, plot.title = element_text(hjust = 0.5)
			, legend.position = "bottom"
			, axis.ticks.y = element_blank()
			, axis.text.x = element_text(size = 12)
			, axis.text.y = element_text(size = 12)
			, axis.title.x = element_text(size = 12)
			, axis.title.y = element_text(size = 12)
			, legend.title = element_text(size = 13, hjust = 0.5)
			, legend.text = element_text(size = 13)
			, panel.grid.major = element_blank()
			, legend.key.size = unit(0.8, "cm")
			, legend.key = element_rect(fill = "white")
			, panel.spacing.y = unit(0.3, "lines")
			, panel.spacing.x = unit(1, "lines")
			, strip.background = element_blank()
			, panel.border = element_rect(colour = "grey"
				, fill = NA
				, size = 0.8
			)
			, strip.text.x = element_text(size = 11
				, colour = "black"
				, face = "bold"
			)
      )
   )
}


