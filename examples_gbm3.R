library(shellpipes)
library(gbm3)
library(survival)
library(satpred); satpredtheme()

commandEnvironments()

set.seed(8888)

### Cross-validation
params_gbm3 <- expand.grid(shrinkage = seq(0.05, 0.1, length.out = 3)
	, n.trees = c(200, 500), n.minobsinnode = 10
	, interaction.depth = 1
)
tuned_gbm3 <- modtune(Surv(Start, Stop, Event) ~ age + alk.phos + ast + chol + edema
	, train_df
	, distribution = "coxph"
	, param_grid = params_gbm3
	, modfun = gbm3.satpred
	, parallelize = TRUE
)
plot(tuned_gbm3)

### Fit model
fit_gbm3 <- modfit(tuned_gbm3, return_data = TRUE)

### Individual survival curves
scurves_gbm3 <- get_indivsurv(fit_gbm3, train_df)
plot(scurves_gbm3)

### Concordance score
concord_gbm3 <- get_survconcord(fit_gbm3)
print(concord_gbm3)

### Permutation variable importance
vimp_gbm3 <- get_varimp(fit_gbm3, type = "perm", newdata = train_df, nrep = 20, modelname = "gbm3")
plot(vimp_gbm3)

saveVars(fit_gbm3
	, scurves_gbm3
	, concord_gbm3
	, vimp_gbm3
)
