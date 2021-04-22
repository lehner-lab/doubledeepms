
#' doubledeepms_structure_metrics
#'
#' Add 3D structure metrics.
#'
#' @param input_file path to MoCHI thermo model fit results (required)
#' @param outpath output path for plots and saved objects (required)
#' @param pdb_file path to PDB file (required)
#' @param pdb_chain_query query chain id (default:A)
#' @param pdb_chain_target target chain id (default:B)
#' @param execute whether or not to execute the analysis (default: TRUE)
#'
#' @return Nothing
#' @export
#' @import data.table
doubledeepms_structure_metrics <- function(
  input_file,
  outpath,
  pdb_file,
  pdb_chain_query = "A",
  pdb_chain_target = "B",
  execute = TRUE
  ){

  #Return if analysis not executed
  if(!execute){
    return()
  }

  #Display status
  message(paste("\n\n*******", "running stage: doubledeepms_structure_metrics", "*******\n\n"))

  #Create output directory
  doubledeepms__create_dir(doubledeepms_dir = outpath)

  #Load free energies
  dg_dt <- fread(input_file)

  #Calculate minimum ligand distances
  dist_dt <- doubledeepms__minimum_interchain_distances_from_PDB(
    input_file = pdb_file,
    chain_query = pdb_chain_query,
    chain_target = pdb_chain_target)
  dist_dt[, Pos_ref := Pos]

  #Get relative SASA
  sasa_dt <- doubledeepms__relative_SASA_from_PDB(
    input_file = pdb_file,
    chain = pdb_chain_query)
  sasa_dt[, Pos_ref := Pos]

  #Merge with free energies
  dg_dt <- merge(dg_dt, dist_dt[,.SD,,.SDcols = names(dist_dt)[names(dist_dt)!="Pos"]], by = "Pos_ref")
  dg_dt <- merge(dg_dt, sasa_dt[,.SD,,.SDcols = names(sasa_dt)[names(sasa_dt)!="Pos"]], by = "Pos_ref")

  #Define position class
  dg_dt[RSASA<0.5, Pos_class := "core"]
  dg_dt[RSASA>=0.5, Pos_class := "surface"]
  dg_dt[scHAmin_ligand<5, Pos_class := "binding_interface"]

  #Save dGs and ddGs
  write.table(dg_dt, 
    file = file.path(outpath, "dg_singles.txt"), 
    quote = F, sep = "\t", row.names = F)
}