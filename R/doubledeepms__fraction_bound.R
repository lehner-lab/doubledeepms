
#' doubledeepms__fraction_bound
#'
#' Calculate fraction bound in 3-state model.
#'
#' @param folding_energy folding energies for 3-state model (required)
#' @param binding_energy binding energies for 3-state model (required)
#' @param RT constant (default:0.001987*(273+24))
#'
#' @return Fraction bound
#' @export
doubledeepms__fraction_bound <- function(
  folding_energy,
  binding_energy,
  RT = 0.001987*(273+24)){
  return(1/(1+exp(binding_energy/RT)*(1+exp(folding_energy/RT))))
}
