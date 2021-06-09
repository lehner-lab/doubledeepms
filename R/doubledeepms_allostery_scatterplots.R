
#' doubledeepms_allostery_scatterplots
#'
#' Plot free energy scatterplots for allosteric mutations
#'
#' @param input_file path to input file (required)
#' @param fitness_list list of folder paths to fitness data (required)
#' @param outpath output path for plots and saved objects (required)
#' @param colour_scheme colour scheme file (required)
#' @param execute whether or not to execute the analysis (default: TRUE)
#'
#' @return Nothing
#' @export
#' @import data.table
doubledeepms_allostery_scatterplots <- function(
  input_file,
  fitness_list,
  outpath,
  colour_scheme,
  execute = TRUE
  ){

  #Return if analysis not executed
  if(!execute){
    return()
  }

  #Display status
  message(paste("\n\n*******", "running stage: doubledeepms_allostery_scatterplots", "*******\n\n"))

  #Create output directory
  doubledeepms__create_dir(doubledeepms_dir = outpath)

  ### Load single mutant free energies
  ###########################

  #Load dg data
  dg_dt <- fread(input_file)

  ###########################
  ### Scatterplots for allosteric sites
  ###########################

  #Free energy scatterplots by protein - all - conf
  for(i in dg_dt[,unique(protein)]){
    plot_dt <- copy(dg_dt)[protein==i][,.(f_dg_pred, f_ddg_pred_sd, b_dg_pred, b_ddg_pred_sd, f_ddg_pred_conf, b_ddg_pred_conf, Pos_class, allosteric, id, Pos_ref)]
    plot_dt <- plot_dt[f_ddg_pred_conf==T & b_ddg_pred_conf==T,]
    plot_dt[, Pos_ref_plot := factor(Pos_ref)]
    #Plot
    d <- ggplot2::ggplot(plot_dt[id!="-0-"],ggplot2::aes(f_dg_pred, b_dg_pred)) +
      ggplot2::geom_point(alpha = 0.5, size = 1, color = "lightgrey") +
      ggplot2::geom_point(data = plot_dt[allosteric==T & Pos_class !="binding_interface",], ggplot2::aes(color = Pos_ref_plot), size = 2) +
      ggplot2::geom_linerange(data = plot_dt[allosteric==T & Pos_class !="binding_interface",], ggplot2::aes(ymin = b_dg_pred-b_ddg_pred_sd*1.96, ymax = b_dg_pred+b_ddg_pred_sd*1.96, color = Pos_ref_plot)) +
      ggplot2::geom_linerange(data = plot_dt[allosteric==T & Pos_class !="binding_interface",], ggplot2::aes(xmin = f_dg_pred-f_ddg_pred_sd*1.96, xmax = f_dg_pred+f_ddg_pred_sd*1.96, color = Pos_ref_plot)) +
      ggplot2::geom_vline(xintercept = 0) +
      ggplot2::geom_vline(data = plot_dt[id=="-0-",], ggplot2::aes(xintercept = f_dg_pred), linetype = 2) +
      ggplot2::geom_hline(yintercept = 0) +
      ggplot2::geom_hline(data = plot_dt[id=="-0-",], ggplot2::aes(yintercept = b_dg_pred), linetype = 2) +
      ggplot2::xlab(expression("Folding "*Delta*"G")) +
      ggplot2::ylab(expression("Binding "*Delta*"G")) +
      ggplot2::labs(color = "Residue\nposition") +
      ggplot2::theme_classic()
    ggplot2::ggsave(file.path(outpath, paste0("ddG_scatter_", i, ".pdf")), d, width = 4, height = 3, useDingbats=FALSE)
    #Plot
    d <- ggplot2::ggplot(plot_dt[id!="-0-"],ggplot2::aes(f_dg_pred, b_dg_pred)) +
      ggplot2::geom_point(alpha = 0.5, size = 1, color = "lightgrey") +
      ggplot2::geom_point(data = plot_dt[allosteric==T & Pos_class !="binding_interface",], ggplot2::aes(color = Pos_ref_plot), size = 2) +
      ggplot2::geom_linerange(data = plot_dt[allosteric==T & Pos_class !="binding_interface",], ggplot2::aes(ymin = b_dg_pred-b_ddg_pred_sd*1.96, ymax = b_dg_pred+b_ddg_pred_sd*1.96, color = Pos_ref_plot)) +
      ggplot2::geom_linerange(data = plot_dt[allosteric==T & Pos_class !="binding_interface",], ggplot2::aes(xmin = f_dg_pred-f_ddg_pred_sd*1.96, xmax = f_dg_pred+f_ddg_pred_sd*1.96, color = Pos_ref_plot)) +
      ggplot2::geom_vline(xintercept = 0) +
      ggplot2::geom_vline(data = plot_dt[id=="-0-",], ggplot2::aes(xintercept = f_dg_pred), linetype = 2) +
      ggplot2::geom_hline(yintercept = 0) +
      ggplot2::geom_hline(data = plot_dt[id=="-0-",], ggplot2::aes(yintercept = b_dg_pred), linetype = 2) +
      ggplot2::xlab(expression("Folding "*Delta*"G")) +
      ggplot2::ylab(expression("Binding "*Delta*"G")) +
      ggplot2::labs(color = "Residue\nposition") +
      ggplot2::theme_classic()
    if(i!="GB1"){
      d <- d + ggplot2::coord_cartesian(ylim = c(-2.5, 2.5), xlim = c(-2.5,2.5))
    }
    ggplot2::ggsave(file.path(outpath, paste0("ddG_scatter_", i, "_xylim.pdf")), d, width = 4, height = 3, useDingbats=FALSE)
  }



}

