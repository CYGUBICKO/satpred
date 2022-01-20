library(shellpipes)
library(ggplot2)
library(glmnetsurv)
library(pcoxtime); pcoxtheme()
library(survival)
set.seed(8888)

## Data
df <- veteran
df$t0 <- 0
## Fit models
glmnet_fit <- glmnetsurv(Surv(t0, time, status)~., data=df, lambda = 0)
glmnet_coef <- coef(glmnet_fit)
glmnet_coef <- data.frame(coef=rownames(glmnet_coef), value=as.vector(glmnet_coef), model="glmnet")
coxph_fit <- coxph(Surv(t0, time, status)~., data=df, method="breslow")
coxph_coef <- coef(coxph_fit)
coxph_coef <- data.frame(coef=names(coxph_coef), value=as.vector(coxph_coef), model="coxph")
pcox_fit <- pcoxtime(Surv(t0, time, status)~., data=df, lambda=0)
pcox_coef <- coef(pcox_fit)
pcox_coef <- data.frame(coef=names(pcox_coef), value=as.vector(pcox_coef), model="pcox")

coef_df <- do.call("rbind", list(glmnet_coef, coxph_coef, pcox_coef))

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

concordance(coxph_fit)[[1]]
pcoxtime::concordScore(pcox_fit)
glmnetsurv::concordScore(glmnet_fit)[[1]]


