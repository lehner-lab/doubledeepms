
#' doubledeepms__aa_properties_pca_WT_loadings
#'
#' Create results folder for analysis script plots and saved objects.
#'
#' @param input_dt data.table with double mutant codes (required)
#' @param aa_properties_file path to amino acid properties file (required)
#' @param selected_identifiers path to file with selected subset of identifiers
#'
#' @return results from PCA
#' @export
doubledeepms__aa_properties_pca_WT_loadings <- function(
  input_dt,
  aa_properties_file,
  selected_identifiers = NULL
  ){
  #AA properties PCA
  exp_pca <- doubledeepms__aa_properties_pca(aa_properties_file = aa_properties_file, selected_identifiers = selected_identifiers)
  #Mut codes
  mut_code <- input_dt[,paste0(WT_AA, Mut)]
  #List of mutation effects on AA properties
  mut_list <- list()
  for(i in unique(c(mut_code))){
    mut_list[[i]] <- exp_pca$x[substr(i, 1, 1),]
  }
  #Single mutant effects on AA properties
  singles_daaprop <- data.table::as.data.table(do.call("rbind", mut_list[mut_code]))
  names(singles_daaprop) <- paste0(names(singles_daaprop), "_score")
  #Append to input data table
  return(cbind(input_dt, singles_daaprop))
}
