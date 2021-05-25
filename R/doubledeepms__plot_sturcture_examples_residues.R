#' doubledeepms__plot_sturcture_examples_residues
#'
#' generate pymol script to get figures of pdb structures highlighting example residues
#'
#' @param pdb_id name of the pdb structure (required)
#' @param pdb_view vector of the pdb sturcture orientation (required)
#' @param domain_pos vector of domain positions to highlight in the heatmap and 3D structure (required)
#' @param ligand_pos vector of domain positions to highlight in the 3D structure (required)
#' @param domain_backbone_color color of the domain backbone in the 3D structure
#' @param ligand_backbone_color color of the ligand backbone in the 3D structure
#' @param domain_residues_color color of the domain residues in the 3D structure
#' @param ligand_residues_color color of the ligand residues in the 3D structure
#' @param zoom_center points to zoom in (negative) or out (positive) the pdb structure in points (Default: 25)
#' @param output_file output file path
#' 
doubledeepms__plot_sturcture_examples_residues <- function(
  pdb_id,
  pdb_view,
  domain_pos,
  ligand_pos,
  domain_backbone_color,
  ligand_backbone_color,
  domain_residues_color,
  ligand_residues_color,
  zoom_center = 25,
  output_file
) {
  
  pymol_script = "reinitialize"
  pymol_script[length(pymol_script)+1] = paste0("fetch ", pdb_id, ", async=0")
  pymol_script[length(pymol_script)+1] = "hide everything"
  pymol_script[length(pymol_script)+1] = "show cartoon, chain A"
  pymol_script[length(pymol_script)+1] = paste0("set cartoon_color, ", gsub("\\#", "0x", domain_backbone_color))
  pymol_script[length(pymol_script)+1] = "show sticks, chain B"
  pymol_script[length(pymol_script)+1] = paste0("set stick_color, ", gsub("\\#", "0x", ligand_backbone_color),  ", chain B")
  pymol_script[length(pymol_script)+1] = "remove resn hoh"
  pymol_script[length(pymol_script)+1] = "set stick_radius, 0.4"
  pymol_script[length(pymol_script)+1] = "set ray_opaque_background, 0"
  pymol_script[length(pymol_script)+1] = "set ray_shadow, 0"
  pymol_script[length(pymol_script)+1] = "set ray_trace_fog, 0"
  pymol_script[length(pymol_script)+1] = "set antialias, 1"
  pymol_script[length(pymol_script)+1] = "bg_color white"
  pymol_script[length(pymol_script)+1] = "zoom center, 25"
  pymol_script[length(pymol_script)+1] = paste0("set_view ", pdb_view)
  
  for (r in 1:length(domain_pos)) {
    pymol_script[length(pymol_script)+1] = paste0("select a", r, ", (chain A & resi ", domain_pos[r], ")") 
    pymol_script[length(pymol_script)+1] = paste0("show sticks, a", r)
    pymol_script[length(pymol_script)+1] = paste0("color ",  gsub("\\#", "0x", domain_residues_color), ", a", r)
  }
  
  for (r in 1:length(ligand_pos)) {
    pymol_script[length(pymol_script)+1] = paste0("select l", r, ", (chain B & resi ", ligand_pos[r], ")") 
    pymol_script[length(pymol_script)+1] = paste0("set stick_color, ",  gsub("\\#", "0x", ligand_residues_color), ", l", r)
  }
  
  pymol_script[length(pymol_script)+1] = "ray 2400,2400"
  pymol_script[length(pymol_script)+1] = paste0("png pymol_structure_highlighted_residues_", paste(domain_pos, collapse="_"), ".png, dpi=600")
  
  write(x = pymol_script, file = output_file)
  
  
}
  
 