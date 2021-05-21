
#' doubledeepms__mann_whitney_U_wrapper_rand
#'
#' Mann-Whitney U test
#'
#' @param vals1 sample 1 (required)
#' @param vals2 sample 2 (required)
#' @param sd1 sample 1 sd (required)
#' @param sd2 sample 2 sd (required)
#' @param n number of randomisations (default:100)
#' @param strictly_positive values are strictly positive (default:TRUE)
#'
#' @return named vector with effect size (AUC) and p-value obtained from coin package
#' @export
#' @import data.table
doubledeepms__mann_whitney_U_wrapper_rand <- function(
  vals1, 
  vals2,
  sd1,
  sd2,
  n=100,
  strictly_positive=TRUE
  ){
	data <- t(mapply(rnorm, mean = c(vals1, vals2), sd = c(sd1, sd2), n = n))
	if(strictly_positive){
		data <- abs(data)
	}
	data_dt <- data.table(data)
	r_list <- list()
	for(i in colnames(data_dt)){
		r_list[[i]] <- doubledeepms__mann_whitney_U_wrapper(
			data_dt[1:length(vals1),.SD[[1]],,.SDcols = i], 
			data_dt[(length(vals1)+1):nrow(data_dt),.SD[[1]],,.SDcols = i])
	}
	result <- c(mean(sapply(r_list, '[', 1)), sd(sapply(r_list, '[', 1)))
	names(result) <- c("effect_size_mean", "effect_size_sd")
	result
}