#' Generate (TDC) ROC and AUC
#'
#' @keywords internal
one_tdcroc <- function(mod
	, newdata
	, pred=NULL
	, pred.type="lp"
	, predict.time=NULL
	, cut.values=NULL
	, method="NNE"
	, lambda=NULL
	, span=NULL
	, window="symmetric", ...) {
	y <- model.extract(model.frame(terms(mod), data=newdata), "response")
	if (NCOL(y)==2) {
		entry <- NULL
		Stime <- y[,1]
		status <- y[,2]
	} else {
		entry <- y[,1]
		Stime <- y[,2]
		status <- y[,3]
	}
	if (is.null(predict.time)) {
		predict.time <- max(Stime, na.rm=TRUE)
	}
	if (is.null(pred)) {
		pred <- predict(mod, newdata=newdata, type=pred.type, ...)
	} else if (is.function(pred)) {
		pred <- do.call(pred, list(mod=mod, newdata=newdata, ...))
	}
	if (is.null(span)) {
		span <- 0.25 * NROW(newdata)^(-0.20)
	}
	out <- tryCatch({
		survivalROC(Stime=Stime
			, status=status
			, entry=entry
			, marker=pred
			, predict.time=predict.time
			, cut.values=cut.values
			, method=method
			, lambda=lambda
			, span=span
			, window=window
		)
	}, error=function(e)return(list(cut.values=NA, TP=NA, FP=NA, predict.time=predict.time, Survival=NA, AUC=NA)))
}

#' Generate AUC and ROC
#'
#' @export

get_tdcroc <- function(mod
	, newdata
	, idvar=NULL
	, time.eval=NULL
	, pred=NULL
	, cut.values=NULL
	, method="NNE"
	, lambda=NULL
	, span=NULL
	, window="symmetric"
	, nreps = 50
	, prop = c(0.025, 0.5, 0.975)
	, parallelize = FALSE
	, nclusters = 1
	, ...
	, modelname=class(mod)[1]) {
	
	if (is.null(idvar)) {
		newdata$..newidvar <- 1:NROW(newdata)
		idvar <- "..newidvar"
	}
	resamples_index <- resamplesTDC(newdata, nreps, idvar, return_index=TRUE, ...)
	if (nreps==1) {
		resamples_index <- list(list(index=1:NROW(newdata), newID=newdata[,idvar]))
		plusone <- 0
	} else {
		resamples_index[nreps+1] <- list(list(index=1:NROW(newdata), newID=newdata[,idvar]))
		plusone <- 1
	}
	if (parallelize) {
		nn <- min(parallel::detectCores(), nclusters)
		if (nn < 2){
			foreach::registerDoSEQ()
		} else{
			cl <-  parallel::makeCluster(nn)
			doParallel::registerDoParallel(cl)
			on.exit(parallel::stopCluster(cl))
		}

		out <- list()
		for (i in 1:length(time.eval)) {
			est <- foreach(x = 1:length(resamples_index), .export=c("one_tdcroc"), .packages="survivalROC") %dopar% {
				index <- resamples_index[[x]]$index
				df <- newdata[index,,drop=FALSE]
				est <- one_tdcroc(mod=mod
					, newdata=df
					, pred=pred
					, predict.time=time.eval[[i]]
					, cut.values=cut.values
					, method=method
					, lambda=lambda
					, span=span
					, window=window
				)
				est
			}
			est <- do.call("rbind", est)
			AUC0 <- as.data.frame(est[, c("predict.time", "Survival")])
			colnames(AUC0) <- c("times", "AUC")
			AUC0[, "model"] <- "Null model"
			AUC1 <- as.data.frame(est[, c("predict.time", "AUC")])
			colnames(AUC1) <- c("times", "AUC")
			AUC1[, "model"] <- modelname
			AUC <- rbind(AUC0, AUC1)
			AUC <- sapply(AUC, unlist, simplify=FALSE)
			AUC <- cbind.data.frame(AUC)
			AUC <- aggregate(AUC~model+times, data=AUC, FUN=function(x){quantile(x, prop, na.rm=TRUE)})
			AUC <- cbind.data.frame(AUC[,c("model", "times")], AUC$AUC)
			colnames(AUC) <- c("model", "times", "lower", "estimate", "upper")
			ROC <- est[nreps+plusone, c("TP", "FP")]
			ROC$times <- AUC$times[1]
			out[[i]] <- list(ROC=as.data.frame(ROC), AUC=AUC)
		}
	} else {
		out <- list()
		for (i in 1:length(time.eval)) {
			est <- lapply(resamples_index, function(x){
				index <- x$index
				df <- newdata[index,,drop=FALSE]
				est <- one_tdcroc(mod=mod
					, newdata=df
					, pred=pred
					, predict.time=time.eval[[i]]
					, cut.values=cut.values
					, method=method
					, lambda=lambda
					, span=span
					, window=window
				)
				return(est)
			})
			est <- do.call("rbind", est)
			AUC0 <- as.data.frame(est[, c("predict.time", "Survival")])
			colnames(AUC0) <- c("times", "AUC")
			AUC0[, "model"] <- "Null model"
			AUC1 <- as.data.frame(est[, c("predict.time", "AUC")])
			colnames(AUC1) <- c("times", "AUC")
			AUC1[, "model"] <- modelname
			AUC <- rbind(AUC0, AUC1)
			AUC <- sapply(AUC, unlist, simplify=FALSE)
			AUC <- cbind.data.frame(AUC)
			AUC <- aggregate(AUC~model+times, data=AUC, FUN=function(x){quantile(x, prop, na.rm=TRUE)})
			AUC <- cbind.data.frame(AUC[,c("model", "times")], AUC$AUC)
			colnames(AUC) <- c("model", "times", "lower", "estimate", "upper")
			ROC <- est[nreps+plusone, c("TP", "FP")]
			ROC$times <- AUC$times[1]
			out[[i]] <- list(ROC=as.data.frame(ROC), AUC=AUC)
		}
	}
	out <- unlist(out, recursive=FALSE)
	AUC <- do.call("rbind", out[names(out) %in% "AUC"])
	rownames(AUC) <- NULL
	ROC <- do.call("rbind", out[names(out) %in% "ROC"])
	rownames(ROC) <- NULL
	out <- list(AUC=AUC, ROC=ROC)
	class(out) <- c("roctdc", class(out))
	return(out)
}

#' Plot AUC and ROC supporting TDC
#'
#' @import ggplot2
#' @export

plot.roctdc <- function(x, ...
	, which.plot=c("both", "auc", "roc")
	, auc.type=c("line", "point")
	, include.null=FALSE
	, facet.roc=TRUE
	, ci=TRUE
	, pos=0.3
	, size=0.2) {
	which.plot <- match.arg(which.plot)
	if (which.plot=="auc") {
		AUC <- x$AUC
		if (!include.null) {
			AUC <- AUC[AUC$model!="Null model", ]
		}
		auc.type <- match.arg(auc.type)
		if (length(unique(AUC$times))==1) auc.type <- "point"
		if (auc.type=="point") {
			p1 <- (ggplot(AUC, aes(x=as.factor(times), y=estimate, colour=model))
				+ geom_pointrange(aes(ymin=lower, ymax=upper, colour=model), position=position_dodge(pos))
			) 
		} else {
			p1 <- (ggplot(AUC, aes(x=times, colour=model, group=model))
				+ geom_line(aes(y=estimate), size=size)
			)
			if (ci) {
				p1 <- (p1
					+ geom_line(aes(y=lower), lty=2, size=size)
					+ geom_line(aes(y=upper), lty=2, size=size)
				)
			}
		}
		p1 <- (p1
			+ labs(x="Time", y="AUC")
		)
	} else if (which.plot=="roc"||which.plot=="both") {
		ROC <- x$ROC
		ROC$times <- round(ROC$times, 2)
		if (which.plot=="both") {
			AUC <- x$AUC
			AUC <- AUC[AUC$model!="Null model", ]
			ROC <- merge(ROC, AUC)
			ROC$times <- paste0(round(ROC$times,2), "; ", round(ROC$estimate,2), "[", round(ROC$lower,2), ",", round(ROC$upper,2), "]")
		}
		ROC <- ROC[order(ROC$TP, ROC$FP),]
		p1 <- (ggplot(ROC, aes(x=FP, y=TP, colour=as.factor(times)))
			+ geom_step()
		)
		if (facet.roc) {
			p1 <- (p1
				+ facet_wrap(~times)
				+ theme(legend.position="none")
			)
		}
		p1 <- (p1
			+  scale_colour_viridis_d(option = "inferno")
			+ labs(x="False positve rate", y="True positve rate")
		)
	} 
	return(p1)
}
