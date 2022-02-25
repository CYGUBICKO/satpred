#' Generate (TDC) AUC based on riskstAUC from risksetROC package
#'
#' @keywords internal
one_tdcauc <- function(mod
	, newdata
	, pred=NULL
	, pred.type="lp"
	, predict.time=NULL
	, method="Cox"
	, span=NULL
	, order=1
	, window="asymmetric"
	, weight="rescale", ...) {
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
		risksetROC::risksetAUC(Stime=Stime
			, status=status
			, entry=entry
			, marker=pred
			, tmax=predict.time
			, method=method
			, span=span
			, window=window
			, weight=weight
			, plot=FALSE
		)
	}, error=function(e)return(list(utimes=predict.time, St=NA, AUC=NA, Cindex=predict.time)))
	return(out)
}

#' Generate (TDC) AUC based on riskstAUC from risksetROC package
#'
#' @export

get_tdcauc <- function(mod
	, newdata
	, idvar=NULL
	, time.eval=NULL
	, pred=NULL
	, pred.type="lp"
	, method="Cox"
	, span=NULL
	, order=1
	, window="asymmetric"
	, weight="rescale"
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
		est <- foreach(x = 1:length(resamples_index), .export=c("one_tdcauc"), .packages="risksetROC") %dopar% {
			index <- resamples_index[[x]]$index
			df <- newdata[index,,drop=FALSE]
			est <- one_tdcauc(mod=mod
				, newdata=df
				, pred=pred
				, pred.type=pred.type
				, predict.time=time.eval
				, method=method
				, span=span
				, order=order
				, window=window
				, ...
			)
			est <- as.data.frame(est)
			est
		}
	} else {
		est <- lapply(resamples_index, function(x){
			index <- x$index
			df <- newdata[index,,drop=FALSE]
			est <- one_tdcauc(mod=mod
				, newdata=df
				, pred=pred
				, pred.type=pred.type
				, predict.time=time.eval
				, method=method
				, span=span
				, order=order
				, window=window
				, ...
			)
			est <- as.data.frame(est)
			return(est)
		})
	}
	est <- do.call("rbind", est)
	est[, "model"] <- modelname
	AUC <- aggregate(AUC~model+utimes, data=est, FUN=function(x){quantile(x, prop, na.rm=TRUE)})
	AUC <- do.call("cbind.data.frame", AUC)
	colnames(AUC) <- c("model", "times", "lower", "estimate", "upper")
	St <- aggregate(St~model+utimes, data=est, FUN=function(x){quantile(x, prop, na.rm=TRUE)})
	St <- do.call("cbind.data.frame", St)
	colnames(St) <- c("model", "times", "lower", "estimate", "upper")
	if (!any(St$times==0)) {
		St0 <- data.frame(model=modelname, times=0, lower=1, estimate=1, upper=1)
		St <- rbind.data.frame(St0, St)
	}
	Cindex <- aggregate(Cindex~model, data=est, FUN=function(x){quantile(x, prop, na.rm=TRUE)})
	Cindex <- do.call("cbind.data.frame", Cindex)
	colnames(Cindex) <- c("model", "lower", "estimate", "upper")
	est <- list(AUC=AUC, Survival=St, Cindex=Cindex)
	class(est) <- "tdcauc"
	return(est)
}

#' Plot AUC and survival curves supporting TDC using risksetROC package
#'
#' @import ggplot2
#' @export

plot.tdcauc <- function(x, ...
	, which.plot=c("auc", "cindex", "surv")
	, auc.type=c("line", "point")
	, include.null=FALSE
	, ci=TRUE
	, pos=0.3
	, size=0.2) {
	which.plot <- match.arg(which.plot)
	if (which.plot=="auc") {
		AUC <- x$AUC
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
	} else if (which.plot=="cindex") {
		Cindex <- x$Cindex
		p1 <- (ggplot(Cindex, aes(x=model, y=estimate))
			+ geom_pointrange(aes(ymin=lower, ymax=upper))
			+ labs(x="", y="C-index")
		)
	} else {
		St <- x$Survival
		p1 <- (ggplot(St, aes(x=times, colour=model, group=model))
			+ geom_line(aes(y=estimate), size=size)
		)
		if (ci) {
			p1 <- (p1
				+ geom_line(aes(y=lower), lty=2, size=size)
				+ geom_line(aes(y=upper), lty=2, size=size)
			)
		}
	}
	return(p1)
}
