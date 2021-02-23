library(shellpipes)
library(survivalmodels)
library(survival)
library(satpred); satpredtheme()

commandEnvironments()

set.seed(8888)

### Cross-validation
params_deepsurv <- expand.grid(dropout = c(0.1, 0.4)
	, learning_rate = c(0.1, 0.01, 0.001, 0.0001)
	, epochs = c(500, 1000)
)
num_nodes <- list(n1 = c(32,32), n2 = c(25, 25), n3 = c(5, 10, 15))
tuned_deepsurv <- modtune(Surv(time, status) ~ ., train_df, num_nodes = num_nodes
	, param_grid = params_deepsurv, modfun = deepsurv.satpred, lr_decay = 0.001
	, parallelize = TRUE, early_stopping = FALSE
)
plot(tuned_deepsurv)

### Fit model
fit_deepsurv <- modfit(tuned_deepsurv, return_data = TRUE)

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
