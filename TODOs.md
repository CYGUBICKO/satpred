add return_data = TRUE to modfit
add is.character to pvimp*

Add
if (inherits(finalModel, "gbm")){
	class(finalModel) <- c("gbm.satpred", class(finalModel), "satpred")
} else {
	class(finalModel) <- c(class(finalModel), "satpred")
}

Add colouring
+ scale_color_manual(breaks = nmods
	, values = rainbow(n = length(nmods))
)

