#' Create K-fold indices for TDC data
#'
#' @export

get_tdcfoldids <- function(ids, k=10) {
	folds <- sample(rep_len(1:k, length(unique(ids))))[as.numeric(as.factor(as.character(ids)))]
	return(folds)
}
