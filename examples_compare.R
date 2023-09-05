library(shellpipes)
library(pec)
library(satpred)

commandEnvironments()

## Load deepsurv model
fit_deepsurv$model <- reticulate::py_load_object("fit_deepsurv.pkl")

## Variable importance
vimp_all <- do.call("rbind", list(vimp_coxph, vimp_rfsrc, vimp_gbm, vimp_deepsurv))
plot(vimp_all)

## Prediction error
prederror_all <- pec(list(coxph=fit_coxph, rfsrc=fit_rfsrc, gbm=fit_gbm, deepsurv=fit_deepsurv)
	, formula = Surv(time, status)~1, data = test_df
)

plotpec(prederror_all)


## Score
score_all <- Score(list(coxph=fit_coxph, rfsrc=fit_rfsrc, gbm=fit_gbm, deepsurv=fit_deepsurv)
	, formula = Surv(time, status)~1, data = test_df
	, plots = "roc", times = c(90, 180, 360)
)

### roc
plot(score_all, type = "roc")

### auc
plot(score_all, type = "auc")

### brier
plot(score_all, type = "brier")
