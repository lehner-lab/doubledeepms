
#' doubledeepms__mochi__construct_feature_matrix
#'
#' Contruct feature matrix for model fitting.
#'
#' @param input_dt Input data.table (required)
#' @param prefix Prefix for feature names (default:empty string)
#' @param mut_ids Mutation id subset (default:NULL)
#'
#' @return Feature matrix data.table
#' @export
#' @import data.table
doubledeepms__mochi__construct_feature_matrix <- function(
  input_dt,
  prefix = "",
  mut_ids = NULL
  ){

  wt_seq <- input_dt[id=="-0-",aa_seq][1]

  #Mutated residues
  input_dt[, mut1 := sapply(strsplit(id, ","), '[', 1)]
  input_dt[, mut2 := sapply(strsplit(id, ","), '[', 2)]
  input_dt[, mut_pos1 := as.integer(substr(mut1, 2, nchar(mut1)-1))]
  input_dt[, mut_pos2 := as.integer(substr(mut2, 2, nchar(mut2)-1))]

  #All mutated positions
  all_mut_pos <- unique(unlist(as.data.frame(input_dt[id!="-0-",.(mut_pos1, mut_pos2)])))
  all_mut_pos <- all_mut_pos[order(all_mut_pos)]

  #Feature matrix
  feature_matrix <- data.table(model.matrix(~., as.data.frame(do.call("rbind", strsplit(input_dt[,var], "")))))

  #Reformat column names
  col_dt <- data.table(name = names(feature_matrix)[-1], wt_seq = wt_seq)
  col_dt[,pos := as.integer(substr(name, 2, nchar(name)-1))]
  col_dt[,new_pos := all_mut_pos[pos]]
  col_dt[,new_name := paste0(substr(wt_seq, new_pos, new_pos), new_pos, substr(name, nchar(name), nchar(name)))]
  names(feature_matrix)[-1] <- paste0(prefix, "_", col_dt[,new_name])
  names(feature_matrix)[1] <- "WT"

  #Subset to desired ids
  if(!is.null(mut_ids)){
    feature_matrix <- feature_matrix[,.SD,,.SDcols = names(feature_matrix)[names(feature_matrix) %in% c("WT", paste0(prefix, "_", mut_ids))]]
  }

  return(feature_matrix)
}
