
#' doubledeepms_thermo_model_results
#'
#' Evaluate thermo model results.
#'
#' @param mochi_outpath path to MoCHI thermo model fit results (required)
#' @param temperature temperature in degrees celcuis (default:24)
#' @param literature_free_energies path to literature free energies (default:NA)
#' @param position_offset residue position offset (default:0)
#' @param folding_energy_max_sd maximum folding energy standard deviation (default:1/(1.96*2))
#' @param binding_energy_max_sd maximum binding energy standard deviation (default:1/(1.96*2))
#' @param outpath output path for plots and saved objects (required)
#' @param colour_scheme colour scheme file (required)
#' @param execute whether or not to execute the analysis (default: TRUE)
#'
#' @return Nothing
#' @export
#' @import data.table
doubledeepms_thermo_model_results <- function(
  mochi_outpath,
  temperature = 24,
  literature_free_energies = NA,
  position_offset = 0,
  folding_energy_max_sd = 1/(1.96*2),
  binding_energy_max_sd = 1/(1.96*2),
  outpath,
  colour_scheme,
  execute = TRUE
  ){

  #Return if analysis not executed
  if(!execute){
    return()
  }

  #Domain name
  domain_name <- rev(unlist(strsplit(basename(outpath), "_")))[1]

  #Display status
  message(paste("\n\n*******", paste("running stage: doubledeepms_thermo_model_results for", domain_name), "*******\n\n"))

  #Create output directory
  doubledeepms__create_dir(doubledeepms_dir = outpath)

  #Constants
  gas_constant <- 0.001987
  RT <- gas_constant*(273+temperature)

  #Load model data
  fitness_dt <- fread(file.path(mochi_outpath, "model_data.txt"))

  #Plot mutation distributions
  doubledeepms__plot_mutation_distributions(
    input_dt = fitness_dt,
    report_outpath = outpath,
    colour_scheme = colour_scheme)

  #Load model results
  pred_dt <- doubledeepms__get_model_results(
    input_folder = mochi_outpath, 
    input_dt = fitness_dt, 
    RT = RT)

  #Call confident ddGs
  pred_dt_conf <- doubledeepms__define_confident_free_energies(
    input_dt = pred_dt, 
    report_outpath = outpath, 
    highlight_colour = colour_scheme[["shade 0"]][[1]],
    folding_energy_max_sd = folding_energy_max_sd, 
    binding_energy_max_sd = binding_energy_max_sd)

  #Plot model performance
  doubledeepms__plot_model_performance(
    input_dt = pred_dt_conf, 
    report_outpath = outpath, 
    highlight_colour = colour_scheme[["shade 0"]][[1]])

  #Plot 2D global epistasis (folding energy vs. folding fitness) 
  doubledeepms__plot_additive_trait_folding(
    mochi_outpath = mochi_outpath,
    input_dt = pred_dt_conf, 
    RT = RT,
    report_outpath = outpath,
    colour_scheme = colour_scheme)

  #Plot 3D global epistasis (folding+binding energies vs. binding fitness) 
  doubledeepms__plot_additive_trait_binding(
    mochi_outpath = mochi_outpath,
    input_dt = pred_dt_conf, 
    RT = RT,
    report_outpath = outpath,
    colour_scheme = colour_scheme)

  #Plot folding versus binding energy coloured by fraction bound showing isochores for arbitrary double mutant
  doubledeepms__plot_isochore_fraction_bound(
    mochi_outpath = mochi_outpath,
    input_dt = pred_dt_conf, 
    RT = RT,
    report_outpath = outpath,
    colour_scheme = colour_scheme)

  #Plot correlation with validation data
  doubledeepms__plot_validation_scatter(
    input_dt = pred_dt_conf, 
    lit_inpath = literature_free_energies, 
    report_outpath = outpath, 
    highlight_colour = colour_scheme[["shade 0"]][[1]],
    RT = RT, 
    position_offset = position_offset)

  #Add id with reference amino acid position
  pred_dt_conf[, id_ref := doubledeepms__get_reference_id(pred_dt_conf[,.(id, mut_order)], position_offset)]

  #Add residue position for singles
  pred_dt_conf[mut_order==1, Pos := as.integer(substr(id, 2, nchar(id)-1))]
  pred_dt_conf[mut_order==1, Pos_ref := as.integer(substr(id_ref, 2, nchar(id_ref)-1))]

  #Save dGs and ddGs
  write.table(pred_dt_conf, 
    file = file.path(outpath, "model_results.txt"), 
    quote = F, sep = "\t", row.names = F)
  write.table(pred_dt_conf[mut_order<=1 & !duplicated(id),.SD,,.SDcols = c(
    "id", 
    "id_ref",
    "Pos",
    "Pos_ref",
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

  #Add ddG significance
  pred_dt_conf[mut_order==1 & !duplicated(id), b_ddg_pred_sig := p.adjust(doubledeepms__pvalue(b_ddg_pred, b_ddg_pred_sd), method = "BH")<0.05]
  pred_dt_conf[mut_order==1 & !duplicated(id), f_ddg_pred_sig := p.adjust(doubledeepms__pvalue(f_ddg_pred, f_ddg_pred_sd), method = "BH")<0.05]

  #Save ids of singles with significant ddGs
  write.table(pred_dt_conf[mut_order==1 & !duplicated(id) & b_ddg_pred_sig==T,id], 
    file = file.path(outpath, "b_ddg_sig_singles_id.txt"), 
    quote = F, sep = "\n", row.names = F, col.names = "id")
  write.table(pred_dt_conf[mut_order==1 & !duplicated(id) & f_ddg_pred_sig==T,id], 
    file = file.path(outpath, "f_ddg_sig_singles_id.txt"), 
    quote = F, sep = "\n", row.names = F, col.names = "id")
}