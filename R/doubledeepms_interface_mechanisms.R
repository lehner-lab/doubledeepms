
#' doubledeepms_interface_mechanisms
#'
#' Plot terface structures and binding free energy heatmaps subsets
#'
#' @param input_folder_prefix folder name where free energy heatmap data is (required)
#' @param domain_name domain name (required)
#' @param outpath output path for plots and saved objects (required)
#' @param colour_scheme colour scheme file (required)
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
  
  
  ###########################
  ### Binding heatmap subsets
  ###########################
  
  doubledeepms__plot_heatmap_subset(
    input_file = heat_df,
    output_file = file.path(outpath, "binding_subset_heatmap_pos_.pdf"),
    width = 10,
    height = 4,
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

