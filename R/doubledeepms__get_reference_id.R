
#' doubledeepms__get_reference_id
#'
#' Get reference id with absolute amino acid residue position.
#'
#' @param input_dt data.table with id and mut_order columns (required)
#' @param position_offset residue position offset (default:0)
#'
#' @return Reference id with absolute amino acid residue position
#' @export
#' @import data.table
doubledeepms__get_reference_id <- function(
  input_dt, 
  position_offset=0
  ){
	id_dt <- copy(input_dt[,.(id, mut_order)])
  id_dt[id!="-0-", mut1 := sapply(strsplit(id, ","), '[', 1)]
  id_dt[mut_order==2, mut2 := sapply(strsplit(id, ","), '[', 2)]
  id_dt[id!="-0-", pos1_ref := as.integer(substr(mut1, 2, nchar(mut1)-1))+position_offset]
  id_dt[mut_order==2, pos2_ref := as.integer(substr(mut2, 2, nchar(mut2)-1))+position_offset]
  id_dt[id!="-0-", mut1_ref := paste0(substr(mut1, 1, 1), pos1_ref, substr(mut1, nchar(mut1), nchar(mut1)))]
  id_dt[mut_order==2, mut2_ref := paste0(substr(mut2, 1, 1), pos2_ref, substr(mut2, nchar(mut2), nchar(mut2)))]
  id_dt[id!="-0-", id_ref := mut1_ref]
  id_dt[mut_order==2, id_ref := paste(mut1_ref, mut2_ref, sep = ",")]
  id_dt[id=="-0-", id_ref := id]  
  return(id_dt[,id_ref])
}
