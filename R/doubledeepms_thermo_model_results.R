
#' doubledeepms_thermo_model_results
#'
#' Evaluate thermo model results.
#'
#' @param mochi_outpath path to MoCHI thermo model fit results (required)
#' @param outpath output path for plots and saved objects (required)
#' @param temperature temperature in degrees celcuis (default:24)
#' @param literature_free_energies path to literature free energies (default:NA)
#' @param position_offset residue position offset (default:0)
#' @param folding_energy_max_sd maximum folding energy standard deviation (default:1/(1.96*2))
#' @param binding_energy_max_sd maximum binding energy standard deviation (default:1/(1.96*2))
#' @param execute whether or not to execute the analysis (default: TRUE)
#'
#' @return Nothing
#' @export
#' @import data.table
doubledeepms_thermo_model_results <- function(
  mochi_outpath,
  outpath,
  temperature = 24,
  literature_free_energies = NA,
  position_offset = 0,
  folding_energy_max_sd = 1/(1.96*2),
  binding_energy_max_sd = 1/(1.96*2),
  execute = TRUE
  ){

  #Return if analysis not executed
  if(!execute){
    return()
  }

  #Display status
  message(paste("\n\n*******", "running stage: doubledeepms_thermo_model_results", "*******\n\n"))

  #Create output directory
  doubledeepms__create_dir(doubledeepms_dir = outpath)

  #Constants
  gas_constant <- 0.001987
  RT <- gas_constant*(273+temperature)

  #Load model data
  fitness_dt <- fread(file.path(mochi_outpath, "model_data.txt"))

  #Load model results
  pred_dt <- doubledeepms__get_model_results(
    input_folder = mochi_outpath, 
    input_dt = fitness_dt, RT = RT)

  #Call confident ddGs
  pred_dt_conf <- doubledeepms__define_confident_free_energies(
    pred_dt, outpath, 
    folding_energy_max_sd = folding_energy_max_sd, 
    binding_energy_max_sd = binding_energy_max_sd)

  #Plot model performance
  doubledeepms__plot_model_performance(pred_dt_conf, outpath)

  #Plot folding additive trait 
  doubledeepms__plot_additive_trait_folding(pred_dt_conf, outpath)

  #Plot correlation with validation data
  doubledeepms__plot_validation_scatter(
    pred_dt_conf, 
    literature_free_energies, 
    outpath, 
    RT = RT, position_offset = position_offset)

  #Add id with reference amino acid position
  pred_dt_conf[, id_ref := doubledeepms__get_reference_id(pred_dt_conf[,.(id, mut_order)], position_offset)]

  #Save dGs and ddGs
  write.table(pred_dt_conf, 
    file = file.path(outpath, "model_results.txt"), 
    quote = F, sep = "\t", row.names = F)
  write.table(pred_dt_conf[mut_order<=1 & !duplicated(id),.SD,,.SDcols = c(
    "id", 
    "id_ref",
    "mut_order", 
    "f_dg_pred", 
    "f_ddg_pred", 
    "f_ddg_pred_sd", 
    "b_dg_pred", 
    "b_ddg_pred", 
    "b_ddg_pred_sd", 
    "f_ddg_pred_conf", 
    "b_ddg_pred_conf")], 
  file = file.path(outpath, "dg_singles.txt"), 
  quote = F, sep = "\t", row.names = F)
}