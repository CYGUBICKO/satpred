#' Generate probabilities for TDC brier score
#'
#' @keywords internal

brierProb <- function(mod, newdata, time.eval, idvar
		, ...
		, coefficients=mod$coefficientsd
		, method=mod$method
	) {
	if(!inherits(mod, "coxph")) {
		class(mod) <- c(class(mod), "coxph")
		mod$x <- NULL
	}
	if (is.null(coefficients)) {
		mod$coefficients <- coef(mod)
		mod$coef <- NULL
	}
	if (is.null(method)) mod$method <- "breslow"
	if (is.null(mod$y)) {
		mod$y <- mod$Y
		mod$Y <- NULL
	}
	time.eval <- union(0, time.eval)
	newdata <- as.data.frame(newdata)
	pred <- lapply(1:NROW(newdata), function(i){survival::survfit(mod, newdata=newdata[i,], se.fit=FALSE)})
	ynames <- deparse(formula(mod)[[2]])
	timevars <- trimws(strsplit(gsub(".*\\(|\\).*", "", ynames), "\\,")[[1]])
	newdata <- newdata[, c(timevars, idvar)]
	N <- length(unique(newdata[, idvar]))
	time.tau <- rep(max(time.eval), N)
	Shat <- sapply(1:N, function(Ni) .shatfunc2(Ni, data = newdata, pred = pred, tpnt = time.eval, tau = time.tau, idvar=idvar))
	obj <- Surv(newdata[, timevars[1]]
		, newdata[, timevars[2]]
		, newdata[, timevars[3]]
	)
	RES <- list(survival.probs = Shat
		, survival.times = time.eval
		, survival.tau = time.tau
		, survival.obj = obj
		, survival.id = newdata[,idvar]
	)
	return(RES)
}

#' Helper function for the TDC brier score
#'
#' @keywords internal

## Directly copied (with some few edits) from LTRCforests package
.shatfunc2 <- function(Ni, data, pred, tpnt, tau, idvar){
	## This function is to compute the estimated survival probability of the Ni-th subject
	id.seu <- data[, idvar] # id
	id.sub <- unique(id.seu)
	
	## the i-th data
	TestData <- data[id.seu == id.sub[Ni], ]
	
	TestT <- c(TestData[1, 1], TestData[, 2])
	TestTIntN <- nrow(TestData)
	
	tpnt <- tpnt[tpnt <= tau[Ni]]
	
	################ Changes at July 29th
	tpntL <- c(TestT, tpnt)
	torder <- order(tpntL)
	tpntLod <- tpntL[torder]
	tlen <- length(tpntLod)
	
	## Compute the estimated survival probability of the Ni-th subject
	Shat_temp <- matrix(0, nrow = 1, ncol = tlen)
	
	r.ID <- findInterval(tpntLod, TestT)
	r.ID[r.ID > TestTIntN] <- TestTIntN
	
	jall <- unique(r.ID[r.ID > 0])
	nj <- length(jall)
	
	## Deal with left-truncation
	Shat_temp[1, r.ID == 0] <- 1
	if(nj == 1){
		## Get the index of the Pred to compute Shat
		II <- which(id.seu == id.sub[Ni])[jall[nj]]
		Shat_i = ipred::getsurv(pred[[II]], tpntLod[r.ID == jall[nj]])
		Shat_temp[1, r.ID == jall[nj]] <- Shat_i / Shat_i[1]
	} else if (nj > 1) {
		# c(1, S_{1}(R_{1}), ..., S_{nj}(R_{nj}))
		ShatR_temp <- matrix(0, nrow = 1, ncol = nj + 1)
		ShatR_temp[1, 1] <- 1
		
		# S_1(L_1), S_2(L_2), S_3(L_3), ..., S_{nj}(L_{nj})
		qL = rep(0, nj)
		for (j in 1:nj){
			## Get the index of the Pred to compute Shat
			II <- which(id.seu == id.sub[Ni])[1] + jall[j] - 1
			Shat_j = ipred::getsurv(pred[[II]], tpntLod[r.ID == jall[j]])
			
			qL[j] <- Shat_j[1]
			# S_{j}(R_{j}), j=1,...nj-1
			jR = ipred::getsurv(pred[[II]], TestT[j + 1])
			ShatR_temp[1, j + 1] = jR / qL[j]
			Shat_temp[1, r.ID == jall[j]] <- Shat_j / qL[j]
		}
		
		ql0 <- which(qL == 0)
		if (length(ql0) > 0){
			if (any(qL > 0)){
				maxqlnot0 <- max(which(qL > 0))
				
				ql0lmax <- ql0[ql0 < maxqlnot0]
				ql0mmax <- ql0[ql0 >= maxqlnot0]
				ShatR_temp[1, ql0lmax + 1] <- 1
				Shat_temp[1, r.ID %in% jall[ql0lmax]] <- 1
				ShatR_temp[1, ql0mmax + 1] <- 0
				Shat_temp[1, r.ID %in% jall[ql0mmax]] <- 0
				# for(j in ql0){
				#   if (j < maxqlnot0) {
				#     ShatR_temp[1, j + 1] <- 1
				#     Shat_temp[1, r.ID == jall[j]] <- 1
				#   } else{
				#     ShatR_temp[1, j + 1] <- 0
				#     Shat_temp[1, r.ID == jall[j]] <- 0
				#   }
				# }
			} else {
				ShatR_b[1, 2:(nj + 1)] <- 0
				Shat_temp[1, r.ID %in% jall] <- 0
			}
		}
		m <- cumprod(ShatR_temp[1, 1:nj])
		for (j in 1:nj){
			Shat_temp[1, r.ID == jall[j]] <- Shat_temp[1, r.ID == jall[j]] * m[j]
		}
	}
	
	# since: tpntLod[torder == 1] == TestData[1, 1]
	return(Shat_temp[1, -match(TestT, tpntLod)])
}

#' Helper function for the TDC brier score
#'
#' @keywords internal

one_ibsTDC <- function(mod
	, newdata
	, time.eval=NULL
	, idvar
	, type=c("BS", "IBS", "both")
	, ...
	, coefficients=mod$coefficientsd
	, method=mod$method) {
	type <- match.arg(type)
	y <- model.extract(model.frame(terms(mod), data=newdata), "response")
	if (is.null(time.eval)) {
		time.eval <- sort(unique(y[,2]))
	}
	time.eval <- union(0, time.eval)
	if (NCOL(y)==2)stop("Only handles Surv(time1, time2, event) format models. For other formats use see ?pec::pec or ?riskRegression::Score!", call.=FALSE)
	pred <- brierProb(mod=mod, newdata=newdata, time.eval=time.eval, idvar=idvar, coefficients=coefficients, method=method, ...)
	
	if (type=="both") {
		BS <- LTRCforests::sbrier_ltrc(obj=y, id=pred$survival.id, pred=pred, type="BS")
		IBS <- LTRCforests::sbrier_ltrc(obj=y, id=pred$survival.id, pred=pred, type="IBS")
		score <- list(BS=BS, IBS=IBS)
		attr(score, "type") <- type
		class(score) <- "ibsTDC"
	} else {
		score <- LTRCforests::sbrier_ltrc(obj=y, id=pred$survival.id, pred=pred, type=type)
		attr(score, "type") <- type
		class(score) <- c("ibsTDC", class(score))
	}
	if (type=="IBS"||type=="both") {
		attr(score, "ibslabel") <- paste0("IBS[", min(time.eval), ";", max(time.eval), "]")
	}
	return(score)
}

#' Create resample indices or data frame for a TDC survival data
#'
#' @export

resamplesTDC <- function(df, nreps, idvar, return_index=TRUE, ...) {
	mergefun <- function(df, indices, idvar, return_index) {
		df$..oldindex <- 1:NROW(df)
		newdf <- data.frame(..xx=indices)
		colnames(newdf) <- idvar
		newdf$..newid <- ave(rep(1, length(newdf[,idvar])), newdf[,idvar], FUN = seq_along)
		out <- merge(newdf, df, by=idvar, sort=FALSE)
		out <- out[order(out$..oldindex), ]
		out <- out[order(out[,idvar], out$..newid), ]
		out[,paste0("_..new",idvar, "_")] <- as.numeric(factor(paste0(out[,idvar], out$..newid)))
		if (return_index) {
			out <- list(index=out$..oldindex, newID=out[, paste0("_..new",idvar, "_")])
		} else {
			out <- out[, colnames(out)[!colnames(out) %in% c("..oldindex", "..newid")]]
		}
		return(out)
	}
	ind <- unique(df[, idvar])
	N <- length(ind)
	out <- replicate(nreps
		, mergefun(df=df, sort(sample(ind,size=N,replace=TRUE)), idvar=idvar, return_index=return_index)
		, simplify=FALSE
	)
	return(out)
}

#' Compute brier or integrated brier score for a TDC survival data
#'
#' @export

get_ibsTDC <- function(mod
	, newdata
	, idvar
	, time.eval=NULL
	, type=c("BS", "IBS", "both")
	, nreps = 50
	, prop = c(0.025, 0.5, 0.975)
	, parallelize = FALSE
	, nclusters = 1
	, ...
	, coefficients=mod$coefficientsd
	, method=mod$method
	, modelname=class(mod)) {
	
	if (missing(idvar)) {
		stop("Specify id column (idvar) or if not in the newdata, add id column corresponding to row number.")
	}
	type <- match.arg(type)
	if (nreps==1) {
		resamples_index <- list(list(index=1:NROW(newdata), newID=newdata[,idvar]))
	} else {
		resamples_index <- resamplesTDC(newdata, nreps, idvar, return_index=TRUE, ...)
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
		est <- foreach(x = 1:length(resamples_index), .export=c("one_ibsTDC", ".shatfunc2", "brierProb")) %dopar% {
			index <- resamples_index[[x]]$index
			newID <- resamples_index[[x]]$newID
			df <- newdata[index,,drop=FALSE]
			df[, idvar] <- newID
			est <- one_ibsTDC(mod=mod
				, newdata=df
				, time.eval=time.eval
				, idvar=idvar
				, type=type
				, coefficients=coefficients
				, method=method
				, ...
			)
			est
		}

	} else {
		est <- lapply(resamples_index, function(x){
			index <- x$index
			newID <- x$newID
			df <- newdata[index,,drop=FALSE]
			df[, idvar] <- newID
			est <- one_ibsTDC(mod=mod
				, newdata=df
				, time.eval=time.eval
				, idvar=idvar
				, type=type
				, coefficients=coefficients
				, method=method
				, ...
			)
			return(est)
		})
	}

	temp <- attr(est[[1]], "ibslabel")
	if (type=="IBS") {
		out <- do.call(rbind, est)
		out <- t(quantile(out, probs = prop, na.rm = TRUE))
		out <- cbind.data.frame(time=temp, out)
		colnames(out) <- c("time", "lower", "estimate", "upper")
		out$model <- modelname
	} else if (type=="BS") {
		out <- do.call(rbind, est)
		out <- aggregate(BScore~Time,out,FUN=function(x){quantile(x, prop, na.rm=TRUE)})
		out <- cbind.data.frame(out$Time, out$BScore)
		colnames(out) <- c("time", "lower", "estimate", "upper")
		out$model <- modelname
	} else {
		out <- do.call("rbind", est)
		IBS <- do.call("rbind", out[, "IBS"])
		IBS <- t(quantile(IBS, probs = prop, na.rm = TRUE))
		IBS <- cbind.data.frame(time=temp, IBS)
		colnames(IBS) <- c("time", "lower", "estimate", "upper")
		IBS$model <- modelname
		
		BS <- do.call("rbind", out[, "BS"])
		BS <- aggregate(BScore~Time,BS,FUN=function(x){quantile(x, prop, na.rm=TRUE)})
		BS <- cbind.data.frame(BS$Time, BS$BScore)
		colnames(BS) <- c("time", "lower", "estimate", "upper")
		BS$model <- modelname
		out <- list(BS=BS, IBS=IBS)
	}
	attr(out, "type") <- type
	class(out) <- c("ibsTDC", class(out))
	return(out)
}

#' Plot brier or integrated brier score estimates
#'
#' @import ggplot2
#' @export

plot.ibsTDC <- function(x, ..., which.plot=c("BS", "IBS", "both"), plotit=TRUE) {
	type <- attr(x, "type")
	which.plot <- match.arg(which.plot)
	if (which.plot!=type && type!="both") {
		which.plot <- type
		message("type in get_ibsTDC is not the same as the specified one in which.plot, sure?")
	}
	if (which.plot=="both") {
		BS <- x$BS
		IBS <- x$IBS
		BS <- (ggplot(BS, aes(x=time, y=estimate, colour=model))
			+ geom_line()
			+ geom_line(aes(y=lower, colour=model), lty=2)
			+ geom_line(aes(y=upper, colour=model), lty=2)
			+ labs(x="Time", y="Brier score")
		)
		IBS <- (ggplot(IBS, aes(x=time, y=estimate, colour=model))
			+ geom_pointrange(aes(ymin=lower, ymax=upper, colour=model))
			+ labs(x="Time", y="Integrated Brier Score")
		)
		if (plotit) {
			print(BS)
			print(IBS)
		}
		invisible(list(BS=BS, IBS=IBS))
	} else if (which.plot=="BS") {
		if (type=="both") x <- x$BS
		BS <- (ggplot(x, aes(x=time, y=estimate, colour=model))
			+ geom_line()
			+ geom_line(aes(y=lower, colour=model), lty=2)
			+ geom_line(aes(y=upper, colour=model), lty=2)
			+ labs(x="Time", y="Brier score")
		)
		if (plotit) {
			print(BS)
		}
		invisible(BS)
	} else {
		if (type=="both") x <- x$IBS
		IBS <- (ggplot(x, aes(x=time, y=estimate, colour=model))
			+ geom_pointrange(aes(ymin=lower, ymax=upper, colour=model))
			+ labs(x="Time", y="Integrated Brier Score")
		)
		if (plotit) {
			print(IBS)
		}
		invisible(IBS)
	}
}
