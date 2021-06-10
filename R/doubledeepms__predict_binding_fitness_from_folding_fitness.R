
#' doubledeepms__predict_binding_fitness_from_folding_fitness
#'
#' Predict binding fitness from folding fitness (for no binding affinity model).
#'
#' @param folding_fitness Folding fitness (required)
#' @param folding_linear_kernel Folding fitness kernel (required)
#' @param folding_linear_bias Folding fitness bias (required)
#' @param binding_linear_kernel Binding fitness kernel (required)
#' @param binding_linear_bias Binding fitness bias (required)
#' @param b_dg_wt Wild-type binding delta G (required)
#' @param RT constant (default:0.001987*(273+24))
#'
#' @return Predicted binding fitness
#' @export
#' @import data.table
doubledeepms__predict_binding_fitness_from_folding_fitness <- function(
  folding_fitness, 
  folding_linear_kernel,
  folding_linear_bias,
  binding_linear_kernel,
  binding_linear_bias,
  b_dg_wt,
  RT = 0.001987*(273+24)
  ){
  Cb <- exp(b_dg_wt/RT)
  return(binding_linear_kernel/(Cb*folding_linear_kernel/(folding_fitness-folding_linear_bias)+1)+binding_linear_bias)
}