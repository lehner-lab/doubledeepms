
#' doubledeepms__relative_SASA_from_PDB
#'
#' Get relative solvent accessible surface area.
#'
#' @param input_file path to PDB file (required)
#' @param chain chain id (default:A)
#'
#' @return data.table with relative solvent accessible surface area
#' @export
#' @import data.table
doubledeepms__relative_SASA_from_PDB <- function(
	input_file,
  chain = "A"
  ){
  
  #load PDB structure
	sink(file = "/dev/null")
  pdb <- bio3d::read.pdb(input_file, rm.alt = TRUE)
	sink()

  ### Atom selections
  ###########################

	#C-alpha atoms
  sele_ca <- bio3d::atom.select(pdb, "calpha", chain = chain, verbose=FALSE)

  ### Get relative SASA
  ###########################

  #Subset to c-alpha atoms
	pdb_sub <- bio3d::trim.pdb(pdb, sele_ca)
  #Result data.table
  result_dt <- data.table(
  	Pos = pdb_sub$atom[,"resno"],
  	RSASA = pdb_sub$atom[,"b"])

  #Return
  return(result_dt)

}

