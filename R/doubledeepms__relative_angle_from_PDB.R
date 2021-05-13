
#' doubledeepms__relative_angle_from_PDB
#'
#' Get relative angle of c-beta to c-alpha vector w.r.t. center of mass to c-alpha vector from PDB file.
#'
#' @param input_file path to PDB file (required)
#' @param chain chain id (default:A)
#'
#' @return data.table with relative angles for all residues 
#' @export
#' @import data.table
doubledeepms__relative_angle_from_PDB <- function(
	input_file,
  chain = "A"
  ){
  
  #load PDB structure
	sink(file = "/dev/null")
  pdb <- bio3d::read.pdb(input_file, rm.alt = TRUE)
	sink()

  ### Atom selections
  ###########################

  #Center of mass
  sele_chain <- bio3d::atom.select(pdb, chain = chain, verbose=FALSE)

	#C-alpha atoms
  sele_ca <- bio3d::atom.select(pdb, "calpha", chain = chain, verbose=FALSE)

  #C-beta atoms
  sele_cb <- bio3d::atom.select(pdb, "cbeta", chain = chain, verbose=FALSE)

  ### Calculate relative angle of c-beta to c-alpha vector w.r.t. centre of mass to c-alpha vector
  ###########################

  #Subset to chain atoms
  pdb_sub <- bio3d::trim.pdb(pdb, sele_chain)

  #Calculate center of mass
  com_xyz <- as.data.table(bio3d::com(pdb_sub$xyz))

  #Subset to c-alpha atoms
	calpha_atoms_xyz <- as.data.table(bio3d::trim.pdb(pdb, sele_ca)[["atom"]])[resid!="GLY",.(resno, x,y,z)]

  #Subset to c-beta atoms
  cbeta_atoms_xyz <- as.data.table(bio3d::trim.pdb(pdb, sele_cb)[["atom"]])[resid!="GLY" & elety=="CB",.(resno, x,y,z)]

  #c-alpha to c-beta vectors at origin
  cbeta_atoms_xyz[, x := x - calpha_atoms_xyz[,x]]
  cbeta_atoms_xyz[, y := y - calpha_atoms_xyz[,y]]
  cbeta_atoms_xyz[, z := z - calpha_atoms_xyz[,z]]

  #c-alpha vectors at origin
  calpha_atoms_xyz[, x := x - com_xyz[,x]]
  calpha_atoms_xyz[, y := y - com_xyz[,y]]
  calpha_atoms_xyz[, z := z - com_xyz[,z]]

  #Calculate relative angle
  for(i in 1:nrow(calpha_atoms_xyz)){
    calpha_atoms_xyz[i, relative_angle := doubledeepms__calculate_angle(
      a = unlist(calpha_atoms_xyz[i,.(x,y,z)]),
      b = unlist(cbeta_atoms_xyz[i,.(x,y,z)]))]
  }

  #Result data.table
  result_dt <- data.table(
  	Pos = calpha_atoms_xyz[,resno],
  	bfactor = calpha_atoms_xyz[,relative_angle])

  #Return
  return(result_dt)

}

