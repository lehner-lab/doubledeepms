
#' doubledeepms_free_energy_scatterplots
#'
#' Plot free energy heatmaps.
#'
#' @param input_list path to MoCHI thermo model fit results (required)
#' @param outpath output path for plots and saved objects (required)
#' @param colour_scheme colour scheme file (required)
#' @param execute whether or not to execute the analysis (default: TRUE)
#'
#' @return Nothing
#' @export
#' @import data.table
doubledeepms_free_energy_scatterplots <- function(
  input_list,
  outpath,
  colour_scheme,
  execute = TRUE
  ){

  #Return if analysis not executed
  if(!execute){
    return()
  }

  #Display status
  message(paste("\n\n*******", "running stage: doubledeepms_free_energy_scatterplots", "*******\n\n"))

  #Create output directory
  doubledeepms__create_dir(doubledeepms_dir = outpath)

  #dG max absolute value
  ddg_max <- 5

  #Load dg data
  dg_list <- list()
  for(protein in names(input_list)){
    temp_dt <- fread(input_list[[protein]])
    temp_dt[, protein := protein]
    temp_dt[abs(f_ddg_pred) > ddg_max, f_ddg_pred_conf := F]
    temp_dt[abs(b_ddg_pred) > ddg_max, b_ddg_pred_conf := F]
    dg_list[[protein]] <- temp_dt
  }
  dg_dt <- rbindlist(dg_list)

  ###########################
  ### Free energy distributions
  ###########################

  #Free energy distributions by protein - all - conf
  plot_dt <- reshape2::melt(copy(dg_dt)[,.(protein, f_ddg_pred, b_ddg_pred, f_ddg_pred_conf, b_ddg_pred_conf)], id = c("protein", "f_ddg_pred_conf", "b_ddg_pred_conf"))
  plot_dt <- plot_dt[(variable=="f_ddg_pred" & f_ddg_pred_conf) | (variable=="b_ddg_pred" & b_ddg_pred_conf),]
  plot_dt[variable=="f_ddg_pred", variable_plot := "Folding"]
  plot_dt[variable=="b_ddg_pred", variable_plot := "Binding"]
  #Plot
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(value, fill = variable_plot)) +
    ggplot2::geom_density(alpha = 0.5) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2) +
    # ggplot2::geom_vline(data = plot_dt[STOP_detrimental==T,.(fitness = median(fitness)),.(pca_type, protein)], ggplot2::aes(xintercept = fitness), linetype = 2) +
    ggplot2::facet_wrap(~protein, scales = "free") +
    ggplot2::xlab(expression(Delta*Delta*"G")) +
    ggplot2::ylab("Density") +
    # ggplot2::coord_cartesian(xlim = c(-3, 6)) +
    ggplot2::labs(fill = "Trait") +
    ggplot2::theme_classic()
  if(!is.null(colour_scheme)){
    d <- d + ggplot2::scale_fill_manual(values = unlist(colour_scheme[["shade 0"]][c(1, 3)]))
  }
  ggplot2::ggsave(file.path(outpath, "ddG_densities_all.pdf"), d, width = 9, height = 2, useDingbats=FALSE)

  #Free energy distributions by protein - GRB2 and PSD95 - conf
  plot_dt <- reshape2::melt(copy(dg_dt)[protein != "GB1",.(protein, f_ddg_pred, b_ddg_pred, f_ddg_pred_conf, b_ddg_pred_conf)], id = c("protein", "f_ddg_pred_conf", "b_ddg_pred_conf"))
  plot_dt <- plot_dt[(variable=="f_ddg_pred" & f_ddg_pred_conf) | (variable=="b_ddg_pred" & b_ddg_pred_conf),]
  plot_dt[variable=="f_ddg_pred", variable_plot := "Folding"]
  plot_dt[variable=="b_ddg_pred", variable_plot := "Binding"]
  #Plot
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(value, fill = variable_plot)) +
    ggplot2::geom_density(alpha = 0.5) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2) +
    # ggplot2::geom_vline(data = plot_dt[STOP_detrimental==T,.(fitness = median(fitness)),.(pca_type, protein)], ggplot2::aes(xintercept = fitness), linetype = 2) +
    ggplot2::facet_wrap(~protein, scales = "free") +
    ggplot2::xlab(expression(Delta*Delta*"G")) +
    ggplot2::ylab("Density") +
    # ggplot2::coord_cartesian(xlim = c(-3, 6)) +
    ggplot2::labs(fill = "Trait") +
    ggplot2::theme_classic()
  if(!is.null(colour_scheme)){
    d <- d + ggplot2::scale_fill_manual(values = unlist(colour_scheme[["shade 0"]][c(1, 3)]))
  }
  ggplot2::ggsave(file.path(outpath, "ddG_densities.pdf"), d, width = 7, height = 2, useDingbats=FALSE)

  ###########################
  ### Free energy scatterplots
  ###########################

  #Free energy scatterplots by protein - all - conf
  plot_dt <- copy(dg_dt)[,.(protein, f_dg_pred, b_dg_pred, f_ddg_pred_conf, b_ddg_pred_conf, Pos_class, id)]
  plot_dt <- plot_dt[f_ddg_pred_conf & b_ddg_pred_conf,]
  plot_dt[, Pos_class_plot := stringr::str_to_title(Pos_class)]
  plot_dt[Pos_class=="binding_interface", Pos_class_plot := "Binding\ninterface"]
  #Plot
  d <- ggplot2::ggplot(plot_dt[id!="-0-"],ggplot2::aes(f_dg_pred, b_dg_pred, colour = Pos_class_plot)) +
    ggplot2::geom_point(alpha = 0.5) +
    ggplot2::geom_vline(xintercept = 0) +
    ggplot2::geom_vline(data = plot_dt[id=="-0-",], ggplot2::aes(xintercept = f_dg_pred), linetype = 2) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::geom_hline(data = plot_dt[id=="-0-",], ggplot2::aes(yintercept = b_dg_pred), linetype = 2) +
    # ggplot2::geom_vline(data = plot_dt[STOP_detrimental==T,.(fitness = median(fitness)),.(pca_type, protein)], ggplot2::aes(xintercept = fitness), linetype = 2) +
    ggplot2::facet_wrap(~protein, scales = "free") +
    ggplot2::xlab(expression("Folding "*Delta*"G")) +
    ggplot2::ylab(expression("Binding "*Delta*"G")) +
    ggplot2::labs(color = "Residue\nposition") +
    ggplot2::theme_classic()
  if(!is.null(colour_scheme)){
    d <- d + ggplot2::scale_colour_manual(values = unlist(colour_scheme[["shade 0"]][c(1, 3, 4)]))
  }
  ggplot2::ggsave(file.path(outpath, "ddG_scatter_all.pdf"), d, width = 9, height = 3, useDingbats=FALSE)

  #Free energy scatterplots by protein - GRB2 and PSD95 - conf
  plot_dt <- copy(dg_dt)[protein != "GB1",.(protein, f_dg_pred, b_dg_pred, f_ddg_pred_conf, b_ddg_pred_conf, Pos_class, id)]
  plot_dt <- plot_dt[f_ddg_pred_conf & b_ddg_pred_conf,]
  plot_dt[, Pos_class_plot := stringr::str_to_title(Pos_class)]
  plot_dt[Pos_class=="binding_interface", Pos_class_plot := "Binding\ninterface"]
  #Plot
  d <- ggplot2::ggplot(plot_dt[id!="-0-"],ggplot2::aes(f_dg_pred, b_dg_pred, colour = Pos_class_plot)) +
    ggplot2::geom_point(alpha = 0.5) +
    ggplot2::geom_vline(xintercept = 0) +
    ggplot2::geom_vline(data = plot_dt[id=="-0-",], ggplot2::aes(xintercept = f_dg_pred), linetype = 2) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::geom_hline(data = plot_dt[id=="-0-",], ggplot2::aes(yintercept = b_dg_pred), linetype = 2) +
    # ggplot2::geom_vline(data = plot_dt[STOP_detrimental==T,.(fitness = median(fitness)),.(pca_type, protein)], ggplot2::aes(xintercept = fitness), linetype = 2) +
    ggplot2::facet_wrap(~protein, scales = "free") +
    ggplot2::xlab(expression("Folding "*Delta*"G")) +
    ggplot2::ylab(expression("Binding "*Delta*"G")) +
    ggplot2::labs(color = "Residue\nposition") +
    ggplot2::theme_classic()
  if(!is.null(colour_scheme)){
    d <- d + ggplot2::scale_colour_manual(values = unlist(colour_scheme[["shade 0"]][c(1, 3, 4)]))
  }
  ggplot2::ggsave(file.path(outpath, "ddG_scatter.pdf"), d, width = 7, height = 3, useDingbats=FALSE)

}

