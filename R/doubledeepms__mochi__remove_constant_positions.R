
#' doubledeepms__mochi__remove_constant_positions
#'
#' Remove invariable positions from sequence strings.
#'
#' @param input_strings vector of sequence strings (required)
#'
#' @return Vector of sequence strings with invariable positions removed
#' @export
#' @import data.table
doubledeepms__mochi__remove_constant_positions <- function(
  input_strings
  ){
  str_dt <- data.table(seq = input_strings)

  #Add mutant codes
  for(i in 1:nchar(str_dt[1,seq])){
    str_dt[, paste0("Mut_pos", i) := paste0(i, substr(seq, i, i))]
    str_dt[, paste0("Mut_seq", i) := substr(seq, i, i)]
  }

  #Mutated positions
  mut_colnames <- names(str_dt)[grep("^Mut_pos", names(str_dt))]
  var_num <- sapply(lapply(as.list(as.data.frame(str_dt[,.SD,,.SDcols = mut_colnames])), table), length)
  mut_colnames <- mut_colnames[var_num!=1]
  str_dt[,seq_new := apply(.SD, 1, paste, collapse=""),,.SDcols = gsub("_pos", "_seq", mut_colnames)]

  return(str_dt[,seq_new])
}
