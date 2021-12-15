
#' doubledeepms__plot_isochore_fraction_bound
#'
#' Plot folding versus binding energy coloured by fraction bound showing isochores for arbitrary double mutant.
#'
#' @param mochi_outpath path to MoCHI thermo model fit results (required)
#' @param input_dt data.table with model free energy estimates (required)
#' @param RT constant (default:0.001987*(273+24))
#' @param report_outpath output path for scatterplots (required)
#' @param colour_scheme colour scheme file (default:ggplot colours)
#'
#' @return Nothing
#' @export
#' @import data.table
doubledeepms__plot_isochore_fraction_bound <- function(
  mochi_outpath,
  input_dt,
  RT = 0.001987*(273+24),
  report_outpath,
  colour_scheme
  ){

  #Number of grid points
  num_grid <- 50

  #Restrict folding and binding energies to [-5, 5]
  folding_energy_range <- c(-5, 5)
  binding_energy_range <- c(-5, 5)

  #Fraction bound grid
  folding_energy_grid <- seq(folding_energy_range[1], folding_energy_range[2], (folding_energy_range[2]-folding_energy_range[1])/num_grid)
  binding_energy_grid <- seq(binding_energy_range[1], binding_energy_range[2], (binding_energy_range[2]-binding_energy_range[1])/num_grid)
  energy_grid_dt <- as.data.table(expand.grid(folding_energy_grid = folding_energy_grid, binding_energy_grid = binding_energy_grid))

  #Predicted fitness
  pred_fitness_list <- doubledeepms__predict_fitness(
    mochi_outpath = mochi_outpath,
    folding_energy = energy_grid_dt[,folding_energy_grid],
    binding_energy = energy_grid_dt[,binding_energy_grid],
    RT = RT)
  #Predicted fitness data table
  pred_fitness_dt <- data.table(
    folding_energy = energy_grid_dt[,folding_energy_grid],
    binding_energy = energy_grid_dt[,binding_energy_grid],
    fraction_bound = pred_fitness_list[["fraction_bound"]])

  isochore_dt_list <- list()
  wt_p_bound <- doubledeepms__fraction_bound(input_dt[id=="-0-",f_dg_pred][1], input_dt[id=="-0-",b_dg_pred][1], RT = RT)
  single1_p_bound <- doubledeepms__fraction_bound(input_dt[id=="-0-",f_dg_pred+1][1], input_dt[id=="-0-",b_dg_pred+1][1], RT = RT)
  single2_p_bound <- doubledeepms__fraction_bound(input_dt[id=="-0-",f_dg_pred+1.5][1], input_dt[id=="-0-",b_dg_pred+2.5][1], RT = RT)
  double_p_bound <- doubledeepms__fraction_bound(input_dt[id=="-0-",f_dg_pred+2.5][1], input_dt[id=="-0-",b_dg_pred+3.5][1], RT = RT)
  for(p_bound in c(wt_p_bound, single1_p_bound, single2_p_bound, double_p_bound)){
    isochore_dt_list[[as.character(p_bound)]] <- data.table(
      folding_energy = folding_energy_grid,
      binding_energy = RT*log((p_bound^-1-1)/(1+exp(folding_energy_grid/RT))),
      fraction_bound = p_bound)
  }
  isochore_dt <- rbindlist(isochore_dt_list)
  isochore_dt[, fraction_bound_factor := as.factor(as.integer(fraction_bound*100))]
  isochore_point_dt <- data.table(
    folding_energy = c(
      input_dt[id=="-0-",(f_dg_pred)][1], 
      input_dt[id=="-0-",(f_dg_pred+1)][1], 
      input_dt[id=="-0-",(f_dg_pred+1.5)][1], 
      input_dt[id=="-0-",(f_dg_pred+2.5)][1]), 
    binding_energy = c(
      input_dt[id=="-0-",(b_dg_pred)][1], 
      input_dt[id=="-0-",(b_dg_pred+1)][1], 
      input_dt[id=="-0-",(b_dg_pred+2.5)][1], 
      input_dt[id=="-0-",(b_dg_pred+3.5)][1]), 
    fraction_bound = c(
      wt_p_bound,
      single1_p_bound,
      single2_p_bound,
      double_p_bound))
  isochore_point_dt[, fraction_bound_factor := as.factor(as.integer(fraction_bound*100))]
  plot_cols <- c("#FFFFFF", "#000000", "#000000", "#000000")
  names(plot_cols) <- as.integer(unique(as.character(isochore_dt[, fraction_bound_factor])))
  d <- ggplot2::ggplot(pred_fitness_dt,ggplot2::aes(folding_energy, binding_energy, fill = fraction_bound)) +
    ggplot2::geom_tile() + 
    ggplot2::geom_line(data = isochore_dt, ggplot2::aes(color = fraction_bound_factor, linetype = fraction_bound_factor)) + 
    ggplot2::geom_point(data = isochore_point_dt, size = 2) + 
    ggplot2::xlab(expression(Delta*"G Folding")) +
    ggplot2::ylab(expression(Delta*"G Binding")) +
    ggplot2::theme_bw() +
    ggplot2::scale_fill_gradientn(
      # colours = c(colour_scheme[["shade 0"]][[1]], colour_scheme[["shade 0"]][[2]], colour_scheme[["shade 0"]][[4]], colour_scheme[["shade 0"]][[3]]),
      colours = c("white", "lightgrey", "darkgrey", "black"),
      aesthetics = "fill") +
    ggplot2::scale_colour_manual(
      values = plot_cols)
  ggplot2::ggsave(file.path(report_outpath, paste0("dG_binding_dG_folding_fraction_bound_isochores.pdf")), d, width = 4.5, height = 3, useDingbats=FALSE)


}
