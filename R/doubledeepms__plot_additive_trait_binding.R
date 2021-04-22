
#' doubledeepms__plot_additive_trait_binding
#'
#' Plot folding additive trait.
#'
#' @param mochi_outpath path to MoCHI thermo model fit results (required)
#' @param input_dt data.table with model free energy estimates (required)
#' @param report_outpath output path for scatterplots (required)
#' @param colour_scheme colour scheme file (default:ggplot colours)
#'
#' @return Nothing
#' @export
#' @import data.table
doubledeepms__plot_additive_trait_binding <- function(
  mochi_outpath,
  input_dt,
  report_outpath,
  colour_scheme
  ){

  #Check foldiong data exists
  if(input_dt[dataset_binding==1,.N]==0){return()}

  #Plot data.table
  plot_dt <- input_dt[dataset_binding==1][!duplicated(id)]

  ### All data
  ###########################

  #Model data.table (for geom_line)
  folding_energy_range <- plot_dt[mut_order>0,range(additive_trait_folding)]
  binding_energy_range <- plot_dt[mut_order>0,range(additive_trait_binding)]
  folding_energy_grid <- seq(folding_energy_range[1], folding_energy_range[2], (folding_energy_range[2]-folding_energy_range[1])/30)
  binding_energy_grid <- seq(binding_energy_range[1], binding_energy_range[2], (binding_energy_range[2]-binding_energy_range[1])/30)
  
  energy_grid_dt <- as.data.table(expand.grid(folding_energy_grid = folding_energy_grid, binding_energy_grid = binding_energy_grid))

  #Predicted fitness
  pred_fitness_list <- doubledeepms__predict_fitness(
    mochi_outpath = mochi_outpath,
    folding_energy = energy_grid_dt[,folding_energy_grid],
    binding_energy = energy_grid_dt[,binding_energy_grid])
  #Predicted fitness data table
  pred_fitness_dt <- data.table(
    additive_trait_folding = rep(energy_grid_dt[,folding_energy_grid], 2),
    additive_trait_binding = rep(energy_grid_dt[,binding_energy_grid], 2),
    observed_fitness = pred_fitness_list[["fitness_binding"]],
    mut_order = rep(c(1, 2), each = energy_grid_dt[,.N]))

  Cairo::CairoPDF(file = file.path(report_outpath, "dG_observed_binding_scatter.pdf"))
  plot3D::persp3D(
    x = folding_energy_grid, 
    y = binding_energy_grid, 
    z = matrix(data=pred_fitness_dt[,observed_fitness], nrow=length(folding_energy_grid), ncol=length(folding_energy_grid)), 
    r=2, shade=0.4, axes=TRUE,scale=TRUE, box=TRUE, nticks=5, ticktype="detailed", colvar=F, col="white", alpha = 0, border=colour_scheme[["shade 0"]][[1]], lwd=0.2,
    xlab = "dG Folding",
    ylab = "dG Binding",
    zlab = "Fitness (Binding)")
  plot3D::scatter3D(
    x = plot_dt[,additive_trait_folding], 
    y = plot_dt[,additive_trait_binding], 
    z = plot_dt[,observed_fitness], 
    add = T, col = "black", alpha = 0.2, cex = 0.2)
  dev.off()

  ### Restrict folding and binding energies to [-5, 10]
  ###########################

  plot_xylim <- c(-5, 10)
  plot_dt <- plot_dt[additive_trait_folding>plot_xylim[1] & additive_trait_folding<plot_xylim[2]]
  plot_dt <- plot_dt[additive_trait_binding>plot_xylim[1] & additive_trait_binding<plot_xylim[2]]

  #Model data.table (for geom_line)
  folding_energy_range <- plot_xylim
  binding_energy_range <- plot_xylim
  folding_energy_grid <- seq(folding_energy_range[1], folding_energy_range[2], (folding_energy_range[2]-folding_energy_range[1])/30)
  binding_energy_grid <- seq(binding_energy_range[1], binding_energy_range[2], (binding_energy_range[2]-binding_energy_range[1])/30)
  
  energy_grid_dt <- as.data.table(expand.grid(folding_energy_grid = folding_energy_grid, binding_energy_grid = binding_energy_grid))

  #Predicted fitness
  pred_fitness_list <- doubledeepms__predict_fitness(
    mochi_outpath = mochi_outpath,
    folding_energy = energy_grid_dt[,folding_energy_grid],
    binding_energy = energy_grid_dt[,binding_energy_grid])
  #Predicted fitness data table
  pred_fitness_dt <- data.table(
    additive_trait_folding = rep(energy_grid_dt[,folding_energy_grid], 2),
    additive_trait_binding = rep(energy_grid_dt[,binding_energy_grid], 2),
    observed_fitness = pred_fitness_list[["fitness_binding"]],
    mut_order = rep(c(1, 2), each = energy_grid_dt[,.N]))

  Cairo::CairoPDF(file = file.path(report_outpath, "dG_observed_binding_scatter_xylim.pdf"))
  plot3D::persp3D(
    x = folding_energy_grid, 
    y = binding_energy_grid, 
    z = matrix(data=pred_fitness_dt[,observed_fitness], nrow=length(folding_energy_grid), ncol=length(folding_energy_grid)), 
    r=2, shade=0.4, axes=TRUE,scale=TRUE, box=TRUE, nticks=5, ticktype="detailed", colvar=F, col="white", alpha = 0, border=colour_scheme[["shade 0"]][[1]], lwd=0.2,
    xlab = "dG Folding",
    ylab = "dG Binding",
    zlab = "Fitness (Binding)")
  plot3D::scatter3D(
    x = plot_dt[,additive_trait_folding], 
    y = plot_dt[,additive_trait_binding], 
    z = plot_dt[,observed_fitness], 
    add = T, col = "black", alpha = 0.2, cex = 0.2)
  dev.off()

}
