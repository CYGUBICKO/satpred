library(shellpipes)
library(survival)
set.seed(8888)

df <- veteran
train_prop <- 0.8
index <- runif(nrow(df), 0, 1) <= train_prop

## Train data
train_df <- droplevels(df[index, ])

## Test data
test_df <- droplevels(df[!index, ])

saveVars(train_df, test_df)
