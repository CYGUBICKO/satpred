library(shellpipes)
library(survival)
library(reticulate)
library(survivalmodels)
library(pec)
library(satpred)

commandEnvironments()

## Variable importance
vimp_all <- do.call("rbind", list(vimp_rfsrc, vimp_gbm, vimp_deepsurv))
plot(vimp_all)

## Prediction error
prederror_all <- pec(list(rfsrc=fit_rfsrc, gbm=fit_gbm)#, deepsurv=fit_deepsurv)
	, formula = Surv(time, status)~1, data = test_df
)
plotpec(prederror_all)
