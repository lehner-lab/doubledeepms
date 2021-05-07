
#' doubledeepms__secondary_structure_from_PDB
#'
#' Get secondary structure.
#'
#' @param input_file path to PDB file (required)
#' @param chain chain id (default:A)
#'
#' @return data.table with secondary structure
#' @export
#' @import data.table
doubledeepms__secondary_structure_from_PDB <- function(
	input_file,
  chain = "A"
  ){
  
  #load PDB structure
	sink(file = "/dev/null")
  pdb <- bio3d::read.pdb(input_file, rm.alt = TRUE)
	sink()

  ### Get secondary structure
  ###########################

  #Get secondary structure
  ss_list <- list()
  for(ss_type in c("helix", "sheet")){
    for(i in 1:length(pdb[[ss_type]][["start"]])){
      ss_list[[paste(ss_type, i, sep = "_")]] <- data.table(
        Pos = pdb[[ss_type]][["start"]][i]:pdb[[ss_type]][["end"]][i],
        SS = ss_type,
        pdb_chain = pdb[[ss_type]][["chain"]][i])
    }
  }
  result_dt <- rbindlist(ss_list)

  #Subset to relevant chain
	result_dt <- result_dt[pdb_chain==chain,.(Pos, SS)]

  #Return
  return(result_dt)

}

