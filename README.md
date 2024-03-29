# Survival Analysis Trianing and PREDiction (satpred)

Installation:

- `remotes::install_github("CYGUBICKO/satpred")`

Example:

```
library(survival)
set.seed(8888)
library(randomForestSRC)
library(pec)
library(riskRegression)
library(satpred); satpredtheme()

## Implemented methods
### rfsrc.satpred: random survival forest
### gbm.satpred: gradient boosted trees
### deepsurv.satpred: survival neural network


df <- veteran
train_prop <- 0.8
index <- runif(nrow(df), 0, 1) <= train_prop

## Train data
train_df <- droplevels(df[index, ])

## Test data
test_df <- droplevels(df[!index, ])

## Survival forest

### Cross-validation
params_rfsrc <- expand.grid(mtry = c(4, 5, 6), nodesize = c(5, 10,15), ntree=c(800, 1000, 1200))
tuned_rfsrc <- modtune(Surv(time, status) ~ ., train_df, param_grid = params_rfsrc
	, modfun = rfsrc.satpred, parallelize = TRUE, seed = 8888
)
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
```
