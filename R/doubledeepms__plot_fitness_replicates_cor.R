
#' doubledeepms__plot_fitness_replicates_cor
#'
#' Plot fitness replicate correlations.
#'
#' @param input_dt list of folder paths to fitness data (required)
#' @param outpath output path for plots and saved objects (required)
#' @param colour_scheme colour scheme file (required)
#'
#' @return Nothing
#' @export
#' @import data.table
doubledeepms__plot_fitness_replicates_cor <- function(
  input_dt, 
  outpath,
  colour_scheme
  ){
  
  fitness_dt <- input_dt

  #Plot replicate toxicity correlation
  num_rep = 3
  for(i in c("GRB2-SH3", "PSD95-PDZ3")){
    for(pca in c("Abundance", "Binding")){
    doubledeepms__ggpairs_binhex(
      input_dt = fitness_dt[Nham_aa %in% c(1,2) & STOP==F & STOP_readthrough==F & protein==i & pca_type==pca,.SD,.SDcols = paste0('fitness', c(1:num_rep), "_uncorr")], 
      output_file = file.path(outpath, paste0(i, "_replicate_scatter_binhex_", pca, ".pdf")),
      xlab = "Fitness",
      ylab = "Fitness",
      width = 3,
      height = 3,
      label_size = 2,
      bins = 20,
      plot_colours = c("lightgrey", "black"))
    }
  }

}
