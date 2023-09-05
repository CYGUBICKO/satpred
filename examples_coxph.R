
library(shellpipes)
library(randomForestSRC)
library(survival)
library(satpred); satpredtheme()

commandEnvironments()

set.seed(8888)

### Cross-validation
tuned_coxph <- modtune(Surv(time, status) ~ ., train_df, param_grid = NULL,
	, modfun = coxph.satpred, parallelize = TRUE
)

print(tuned_coxph$besTune)
plot(tuned_coxph)

### Fit model
fit_coxph <- modfit(tuned_coxph, return_data = FALSE)

### Individual survival curves
scurves_coxph <- get_indivsurv(fit_coxph, train_df)
plot(scurves_coxph)

### Concordance score
concord_coxph <- get_survconcord(fit_coxph)
print(concord_coxph)

### Permutation variable importance
vimp_coxph <- get_varimp(fit_coxph, type = "perm", newdata = train_df, nrep = 20, modelname = "coxph")
plot(vimp_coxph)

saveVars(fit_coxph
	, scurves_coxph
	, concord_coxph
	, vimp_coxph
)
