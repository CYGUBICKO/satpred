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
	, parallelize = TRUE, early_stopping = TRUE
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
