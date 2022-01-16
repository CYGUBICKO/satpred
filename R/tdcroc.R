one_tdcroc <- function(mod
	, newdata
	, pred=NULL
	, predict.time=NULL
	, cut.values=NULL
	, method="NNE"
	, lambda=NULL
	, span=NULL
	, window="symmetric") {
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
		pred <- predict(mod, newdata=newdata, type="lp")
	} else if (is.function(pred)) {
		pred <- do.call(pred, list(newdata=newdata, type="lp"))
	}
	if (is.null(span)) {
		span <- 0.25 * NROW(newdata)^(-0.20)
	}
	out <- try(
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
	)
}

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
	} else {
		resamples_index[nreps+1] <- list(list(index=1:NROW(newdata), newID=newdata[,idvar]))
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
			AUC <- aggregate(AUC~model+times, data=AUC, FUN=function(x){quantile(x, prob)})
			AUC <- cbind.data.frame(AUC[,c("model", "times")], AUC$AUC)
			colnames(AUC) <- c("model", "times", "lower", "estimate", "upper")
			ROC <- est[nreps+1, c("cut.values", "TP", "FP")]
			ROC$times <- AUC$times[1]
			out[[i]] <- list(ROC=ROC, AUC=AUC)
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
			AUC <- aggregate(AUC~model+times, data=AUC, FUN=function(x){quantile(x, prop)})
			AUC <- cbind.data.frame(AUC[,c("model", "times")], AUC$AUC)
			colnames(AUC) <- c("model", "times", "lower", "estimate", "upper")
			ROC <- est[nreps+1, c("cut.values", "TP", "FP")]
			ROC$times <- AUC$times[1]
			out[[i]] <- list(ROC=ROC, AUC=AUC)
		}
	}
	out <- unlist(out, recursive=FALSE)
#	out <- do.call("rbind", out[names(out) %in% "AUC"])
	out <- out[names(out) %in% "ROC"]
	return(out)
}
