library(shellpipes)
library(randomForestSRC)
library(survival)
library(satpred); satpredtheme()

commandEnvironments()

set.seed(8888)

### Cross-validation
params_rfsrc <- expand.grid(mtry = c(2, 3), nodesize = seq(2,20,length.out=5)
#	, nodedepth=c(0, 500, 1000)
)

# ll <- rfsrc.satpred(formula=Surv(time, status) ~ ., train_df=train_df, param_grid=params_rfsrc, forest=TRUE)
# 
# quit()
tuned_rfsrc <- modtune(Surv(time, status) ~ ., train_df, param_grid = params_rfsrc,
	, modfun = rfsrc.satpred, nodedepth=NULL, forest=TRUE, parallelize = TRUE, seed = 8888
)

print(tuned_rfsrc$besTune)
plot(tuned_rfsrc)

### Fit model
fit_rfsrc <- modfit(tuned_rfsrc, return_data = FALSE)

### Individual survival curves
scurves_rfsrc <- get_indivsurv(fit_rfsrc, train_df)
plot(scurves_rfsrc)

### Concordance score
concord_rfsrc <- get_survconcord(fit_rfsrc)
print(concord_rfsrc)

### Permutation variable importance
vimp_rfsrc <- get_varimp(fit_rfsrc, type = "perm", newdata = train_df, nrep = 20, modelname = "rfsrc")
plot(vimp_rfsrc)

saveVars(fit_rfsrc
	, scurves_rfsrc
	, concord_rfsrc
	, vimp_rfsrc
)
