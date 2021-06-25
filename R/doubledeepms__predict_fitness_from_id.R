
#' doubledeepms__predict_fitness_from_id
#'
#' Predict fitness from id.
#'
#' @param mochi_outpath path to MoCHI thermo model fit results (required)
#' @param mut_id mutation id (required)
#' @param RT constant (default:0.001987*(273+24))
#' @param fitness_type fitness type, either folding or binding (default:folding)
#'
#' @return fitness predictions
#' @export
#' @import data.table
doubledeepms__predict_fitness_from_id <- function(
  mochi_outpath,
  mut_id,
  RT = 0.001987*(273+24),
  fitness_type = "folding"
  ){
  #Load coefficients
  coef_dt <- fread(file.path(mochi_outpath, "model_weights_0.txt"))

  #Dictionaries
  folding_coef_list <- coef_dt[,folding_coefficient]
  names(folding_coef_list) <- coef_dt[,id]
  binding_coef_list <- coef_dt[,binding_coefficient]
  names(binding_coef_list) <- coef_dt[,id]

  #Mutation ids
  mut_ids <- unlist(strsplit(mut_id, ","))

  #Predict fitness
  pred_list <- doubledeepms__predict_fitness(
    mochi_outpath = mochi_outpath, 
    folding_energy = folding_coef_list[["WT"]] + sum(folding_coef_list[mut_ids], na.rm = T),
    binding_energy = binding_coef_list[["WT"]] + sum(binding_coef_list[mut_ids], na.rm = T),
    RT = RT)

  #Return
  return(pred_list[[paste0("fitness_", fitness_type)]])
}