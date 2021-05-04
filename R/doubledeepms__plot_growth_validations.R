
#' doubledeepms__plot_fitness_replicates_cor
#'
#' Create results folder for analysis script plots and saved objects.
#'
#' @param fitness_path list of folder paths to fitness data (required)
#' @param outpath output path for plots and saved objects (required)
#' @param val_inpath  path to experimental fitness validations (required)
#' @param execute whether or not the system command will be executed (required)
#'
#' @return Nothing
#' @export
doubledeepms__plot_growth_validations <- function(
  fitness_path, 
  outpath,
  val_inpath,
  execute = TRUE, 
  message = NULL){
  
  #Return if analysis not executed
  if(!execute){
    return()
  }
  
  #Display status
  message(paste("\n\n*******", "running stage: doubledeepms__plot_fitness_replicates_cor", "*******\n\n"))
  
  #Create output directory
  doubledeepms__create_dir(doubledeepms_dir = outpath)
  
  
}






