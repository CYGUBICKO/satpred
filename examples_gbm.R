library(shellpipes)
library(gbm)
library(survival)
library(satpred); satpredtheme()

commandEnvironments()

set.seed(8888)

### Cross-validation
params_gbm <- expand.grid(shrinkage = seq(0.05, 0.1, length.out = 3)
	, n.trees = c(800, 1000), n.minobsinnode = c(10, 15)
	, interaction.depth = c(2, 4, 8)
)
tuned_gbm <- modtune(Surv(time, status) ~ ., train_df, distribution = "coxph", param_grid = params_gbm
	, modfun = gbm.satpred, parallelize = TRUE
)
plot(tuned_gbm)

### Fit model
fit_gbm <- modfit(tuned_gbm, return_data = TRUE)

### Individual survival curves
scurves_gbm <- get_indivsurv(fit_gbm, train_df)
plot(scurves_gbm)

### Concordance score
concord_gbm <- get_survconcord(fit_gbm)
print(concord_gbm)

### Permutation variable importance
vimp_gbm <- get_varimp(fit_gbm, type = "perm", newdata = train_df, nrep = 20, modelname = "gbm")
plot(vimp_gbm)

saveVars(fit_gbm
	, scurves_gbm
	, concord_gbm
	, vimp_gbm
)
