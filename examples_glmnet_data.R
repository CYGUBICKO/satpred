library(shellpipes)
library(ggplot2)
library(glmnet)
library(pcoxtime); pcoxtheme()
library(survival)
set.seed(8888)

nobs <- 100; nvars <- 15
xvec <- rnorm(nobs * nvars)
xvec[sample.int(nobs * nvars, size = 0.4 * nobs * nvars)] <- 0
x <- matrix(xvec, nrow = nobs)  # non-sparse x
x_sparse <- Matrix::Matrix(xvec, nrow = nobs, sparse = TRUE)  # sparse x

# create start-stop data response
beta <- rnorm(5)
fx <- x_sparse[, 1:5] %*% beta / 3
ty <- rexp(nobs, drop(exp(fx)))
tcens <- rbinom(n = nobs, prob = 0.3, size = 1)
starty <- runif(nobs)
yss <- Surv(starty, starty + ty, tcens)

## Fit models
glmnet_fit <- glmnet(x, yss, family = "cox", lambda = 0)
glmnet_coef <- coef(glmnet_fit)
glmnet_coef <- data.frame(coef=rownames(glmnet_coef), value=as.vector(glmnet_coef), model="glmnet")
coxph_fit <- coxph(yss ~ x, method="breslow")
coxph_coef <- coef(coxph_fit)
coxph_coef <- data.frame(coef=names(coxph_coef), value=as.vector(coxph_coef), model="coxph")
pcox_fit <- pcoxtime(yss ~ x, lambda=0)
pcox_coef <- coef(pcox_fit)
pcox_coef <- data.frame(coef=names(pcox_coef), value=as.vector(pcox_coef), model="pcox")

coef_df <- do.call("rbind", list(glmnet_coef, coxph_coef, pcox_coef))
coef_df$coef <- gsub("V", "x", coef_df$coef)

p1 <- (ggplot(coef_df, aes(x=reorder(coef, -value), y=value, col=model))
	+ geom_point(aes(shape=model), alpha=0.2)
	+ scale_colour_manual(breaks = c("pcox", "coxph", "glmnet")
		, values=c("pcox"="red", "coxph"="black", "glmnet"="blue")
	)
	+ scale_shape(guide = FALSE)
	+ labs(colour="Model")
	+ coord_flip()
)
print(p1)


