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
	p0 <- ggplot(plot_df, aes(x = time, group = id), colour = lcol) + labs(x = "Time")
	if (type == "surv"){
		p1 <- (p0 
			+ geom_step(aes(y = avesurv), colour = "red")
			+ geom_step(aes(y = surv), size = lsize, lty = ltype) + labs(y = "Survival prob.")
		)
		if (compare){
			p1 <- p0 + geom_step(aes(y = surv, col = compare_lab), size = lsize) + labs(y = "Survival prob.")
		}
	} else {
		p1 <- (p0 
			+ geom_step(aes(y = avechaz), colour = "red")
			+ geom_step(aes(y = cumhaz), size = lsize, lty = ltype) + labs(y = "Cumulative hazard")
		)
		if (compare){
			p1 <- p0 + geom_step(aes(y = cumhaz, col = compare_lab), size = lsize) + labs(y = "Cumulative hazard")
		}
	}
	return(p1)
}
