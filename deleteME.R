library(shellpipes)
library(survivalmodels)
library(survival)
library(satpred); satpredtheme()

commandEnvironments()

set.seed(8888)

### Cross-validation
params_deepsurv <- expand.grid(dropout = 0.1
	, learning_rate = 0.1
	, epochs = 1
)
num_nodes <- list(n1 = c(32,32))
tuned_deepsurv <- modtune(Surv(time, status) ~ ., train_df, num_nodes = num_nodes
	, param_grid = params_deepsurv, modfun = deepsurv.satpred
	, parallelize = FALSE, early_stopping = TRUE
)
plot(tuned_deepsurv)

### Fit model
fit_deepsurv <- modfit(tuned_deepsurv, return_data = TRUE)
fit_deepsurv$model


### Individual survival curves
scurves_deepsurv <- get_indivsurv(fit_deepsurv, train_df)
plot(scurves_deepsurv)

### Concordance score
concord_deepsurv <- get_survconcord(fit_deepsurv)
print(concord_deepsurv)

### Permutation variable importance
vimp_deepsurv <- get_varimp(fit_deepsurv, type = "perm", newdata = train_df, nrep = 20, modelname = "deepsurv")
plot(vimp_deepsurv)

saveVars(fit_deepsurv
	, scurves_deepsurv
	, concord_deepsurv
	, vimp_deepsurv
)

quit()

###  we can't pickle the $model component from this object
##  for mysterious reasons ... 
m <- fit_deepsurv$model
try(py_save_object(m, "fit_deepsurv.pkl") )

## 'simple' fits can be saved and restored, but only if
## we separately pickle and restore the $model component
ss <- simsurvdata()
dd <- deepsurv(data=ss,learning_rate = 0.1)
py_save_object(dd$model, "fit_deepsurv.pkl")
save("dd", file="tmp.rda")
rm(dd)
load("tmp.rda")
try(predict(dd))
pp <- py_load_object("fit_deepsurv.pkl")
dd$model <- pp
predict(dd)

## if I pass stuff to modfit.satpred
##  that runs a 'simple' deepsurv call without tuning
##  then I can successfully safe
unlink("checkpoint.pkl")
fit_deepsurv <- modfit.satpred(
    list(modelfun=deepsurv,
         modelargs=list(data=ss)))
file.exists("checkpoint.pkl") ##
m <- py_load_object("checkpoint.pkl")

fit_deepsurv <- modfit.satpred(
    list(modelfun=deepsurv,
         modelargs=list(data=ss)))

## I can also save the $model component directly
py_save_object(fit_deepsurv$model, "fit_deepsurv2.pkl")
