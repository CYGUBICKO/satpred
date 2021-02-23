library(shellpipes)
library(survival)
library(survivalmodels)
library(satpred)

commandEnvironments()

names(fit_deepsurv$model)

quit()
## Does not work
scurves_deepsurv <- get_indivsurv(fit_deepsurv, train_df)
plot(scurves_deepsurv)
