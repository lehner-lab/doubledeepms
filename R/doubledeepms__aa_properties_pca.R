
#' doubledeepms__aa_properties_pca
#'
#' Create results folder for analysis script plots and saved objects.
#'
#' @param aa_properties_file path to amino acid properties file (required)
#' @param selected_identifiers path to file with selected subset of identifiers
#' @param return_evidences return AA properties evidences (default:F)
#' @param return_matrix_only only return AA properties matrix (default:F)
#'
#' @return results from PCA
#' @export
doubledeepms__aa_properties_pca <- function(
  aa_properties_file,
  selected_identifiers = NULL,
  return_evidences = F,
  return_matrix_only = F){
  #Load amino acid properties
  aa_df <- read.table(file=aa_properties_file, sep = "\t", header = T, stringsAsFactors = F, row.names = 1)
  #Subset to selected amino acid properties if supplied
  if(!is.null(selected_identifiers)){
    aa_df <- aa_df[rownames(aa_df) %in% selected_identifiers | grepl("BBSUPPL", rownames(aa_df)),]
  }
  #Evidences list
  aa_evidences <- as.list(aa_df[,1])
  names(aa_evidences) <- rownames(aa_df)
  #Reformat and add AA identities
  aa_df <- aa_df[,-1]
  #Remove properties with NAs
  aa_evidences <- aa_evidences[apply(is.na(aa_df), 1, sum)==0]
  aa_df <- aa_df[apply(is.na(aa_df), 1, sum)==0,]
  #Normalise
  aa_mat <- scale(t(aa_df))
  #PCA on all marks
  if(return_matrix_only){
    return(aa_mat)
  }
  if(return_evidences){
    return(list(
      PCA = prcomp(aa_mat),
      evidences = aa_evidences))
  }else{
    return(prcomp(aa_mat))
  }
}

