
#' doubledeepms__predict_fitness
#'
#' Predict fitness from folding and binding energies for 3-state model.
#'
#' @param mochi_outpath path to MoCHI thermo model fit results (required)
#' @param folding_energy folding energies for 3-state model (required)
#' @param binding_energy binding energies for 3-state model (required for binding fitness prediction)
#' @param RT constant (default:0.001987*(273+24))
#'
#' @return A list with folding and binding fitness predictions
#' @export
#' @import data.table
doubledeepms__predict_fitness <- function(
  mochi_outpath,
  folding_energy,
  binding_energy,
  RT = 0.001987*(273+24)
  ){
  #Load model
  modpar <- fread(file.path(mochi_outpath, "model_parameters_0.txt"), header = F)[,V1]
  modpar_list <- as.list(as.numeric(modpar[seq(2, length(modpar), 2)]))
  names(modpar_list) <- modpar[seq(1, length(modpar), 2)]

  #Predicted fitness
  pred_list <- list()

  #Folding fitness
  if(!is.null(folding_energy)){
    pred_list[["fraction_folded"]] <- doubledeepms__fraction_folded(folding_energy, RT)
    pred_list[["fitness_folding"]] <- pred_list[["fraction_folded"]] * modpar_list[["folding_linear_kernel"]] + modpar_list[["folding_linear_bias"]]
    if("fitness_scaler_bias" %in% names(modpar_list) & "fitness_scaler_kernel" %in% names(modpar_list)){
      pred_list[["fitness_folding"]] <- (pred_list[["fitness_folding"]] - modpar_list[["fitness_scaler_bias"]]) / modpar_list[["fitness_scaler_kernel"]]
    }
  }

  #Binding fitness
  if(!is.null(folding_energy) & !is.null(binding_energy) & length(folding_energy) == length(binding_energy)){
    pred_list[["fraction_bound"]] <- doubledeepms__fraction_bound(folding_energy, binding_energy, RT)
    pred_list[["fitness_binding"]] <- pred_list[["fraction_bound"]] * modpar_list[["binding_linear_kernel"]] + modpar_list[["binding_linear_bias"]]
    if("fitness_scaler_bias" %in% names(modpar_list) & "fitness_scaler_kernel" %in% names(modpar_list)){
      pred_list[["fitness_binding"]] <- (pred_list[["fitness_binding"]] - modpar_list[["fitness_scaler_bias"]]) / modpar_list[["fitness_scaler_kernel"]]
    }
  }

  #Return
  return(pred_list)
}