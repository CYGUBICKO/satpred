library(shellpipes)
library(LTRCforests)
library(dplyr)
set.seed(8888)

df <- pbcsample
train_prop <- 0.6
index <- (df
	%>% group_by(ID)
	%>% mutate(index=runif(1, 0, 1) <= train_prop)
	%>% ungroup()
	%>% pull(index)
)

## Train data
train_df <- df #droplevels(df[index, ])

## Test data
test_df <- droplevels(df[!index, ])
str(test_df)

saveVars(train_df, test_df)
