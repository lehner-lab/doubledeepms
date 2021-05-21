
#' doubledeepms__mann_whitney_U_wrapper
#'
#' Mann-Whitney U test
#'
#' @param vals1 sample 1 (required)
#' @param vals2 sample 2 (required)
#'
#' @return named vector with effect size (AUC) and p-value obtained from coin package
#' @export
doubledeepms__mann_whitney_U_wrapper <- function(
  vals1, 
  vals2
  ){
  if(length(vals1)==0 | length(vals2)==0){
    temp_test<-c(NA, NA)
  }else{
    g = factor(c(rep("GroupA", length(vals1)), rep("GroupB", length(vals2))))
    v = c(vals1, vals2)
    g_wt<-coin::wilcox_test(v ~ g)
    suppressWarnings(temp_test<-c(wilcox.test(vals1, vals2)$statistic/(as.numeric(length(vals1))*as.numeric(length(vals2))),
                 coin::pvalue(g_wt)))
  }
  names(temp_test)<-c('effect_size', 'p_value')
  temp_test    
}