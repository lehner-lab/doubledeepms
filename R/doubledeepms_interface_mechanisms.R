
#' doubledeepms_interface_mechanisms
#'
#' Plot terface structures and binding free energy heatmaps subsets
#'
#' @param input_folder_prefix folder name where free energy heatmap data is (required)
#' @param domain_name domain name (required)
#' @param outpath output path for plots and saved objects (required)
#' @param colour_scheme colour scheme file (required)
#' @param mut_subset_list list of aa mutations to highlight in the heatmap (required)
#' @param pos_subset_list list of domain positions to highlight in the heatmap and 3D structure (required)
#' @param ligand_pos_list list of domain positions to highlight in the 3D structure (required)
#' @param pdb_id name of the pdb file (required)
#' @param pdb_view vector with pdb coordinates to orientate the 3D structure (required)
#' @param zoom_pdb point to zoom into the center of the pdb structure (default: 25)
#' @param plot_width heatmap plot width in inches (default:10)
#' @param plot_height heatmap plot height in inches (default:4)
#' @param execute whether or not to execute the analysis (default: TRUE)
#'
#' @return Nothing
#' @export
#' @import data.table
doubledeepms_interface_mechanisms <- function(
  base_dir,
  input_folder_prefix = "006_doubledeepms_free_energy_heatmaps",
  domain_name,
  outpath,
  colour_scheme,
  mut_subset_list,
  pos_subset_list,
  ligand_pos_list,
  pdb_id,
  pdb_view,
  zoom_pdb,
  plot_width = 10,
  plot_height = 4,
  execute = TRUE
){
  
  #Return if analysis not executed
  if(!execute){
    return()
  }
  
  #Display status
  message(paste("\n\n*******", paste("running stage: doubledeepms_interface_mechanisms for", domain_name), "*******\n\n"))
  
  #Create output directory
  outpath <- paste0(outpath, "_", domain_name)
  doubledeepms__create_dir(doubledeepms_dir = outpath)
  
  #Load free energy heatmap info
  heat_df <- fread(file.path(base_dir, paste0(input_folder_prefix, "_", domain_name), "binding_heatmap_data.txt"))
  heat_df[, Pos := as.numeric(unlist(strsplit(x, "_"))[[1]]), by=x]
  
  ###########################
  ### Binding heatmap subsets
  ###########################
  
  for (x in 1:length(mut_subset_list)){
    mut_subset <- mut_subset_list[[x]]
    pos_subset <- pos_subset_list[[x]]
    
    #Adjust weight and height of the plot for each different subset
    w <- 0.5 + 9*length(pos_subset_list[[x]])/max(heat_df$Pos)
    h <- 0.5 +4*length(mut_subset_list[[x]])/length(unique(heat_df$y))
    
    doubledeepms__plot_heatmap_subset(
      input_file = heat_df,
      output_file = file.path(outpath, paste0("binding_subset_heatmap_pos_", paste(pos_subset, collapse="_"), ".pdf")),
      mut_subset = mut_subset,
      pos_subset = pos_subset,
      width = w,
      height = h,
      plot_title = "",
      colour_clip = 4,
      colour_mid = "lightgrey",
      colour_low = colour_scheme[["shade 0"]][[3]],
      colour_high = colour_scheme[["shade 0"]][[1]],
      xaxis_angle = 0,
      xaxis_hjust = 0.5, 
      xaxis_size = 5, 
      na_colour = "white",
      input_matrix_point = TRUE
    )
  }
  
  
  ###########################
  ### Structure plots
  ###########################
  
  domain_pos <- do.call("c", lapply(pos_subset_list, function(x){x}))
  ligand_pos <- do.call("c", lapply(ligand_pos_list, function(x){x}))
  
  doubledeepms__plot_sturcture_examples_residues(
    pdb_id = pdb_id,
    pdb_view = pdb_view,
    domain_pos = domain_pos,
    ligand_pos = ligand_pos,
    domain_backbone_color = colour_scheme[["shade 5"]][[1]],
    ligand_backbone_color = colour_scheme[["shade 5"]][[2]],
    domain_residues_color = colour_scheme[["shade 0"]][[4]],
    ligand_residues_color = colour_scheme[["shade 5"]][[4]],
    zoom_center = zoom_pdb,
    output_file = file.path(outpath, paste0("pymol_script_highlighted_residues_", paste(domain_pos, collapse="_"), ".txt"))
  )
  
}

