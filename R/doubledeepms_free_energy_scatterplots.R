
#' doubledeepms_free_energy_scatterplots
#'
#' Plot free energy heatmaps.
#'
#' @param input_list path to MoCHI thermo model fit results (required)
#' @param input_MSA_list path to MSA frequencies data (required)
#' @param outpath output path for plots and saved objects (required)
#' @param colour_scheme colour scheme file (required)
#' @param execute whether or not to execute the analysis (default: TRUE)
#'
#' @return Nothing
#' @export
#' @import data.table
doubledeepms_free_energy_scatterplots <- function(
  input_list,
  input_MSA_list,
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

  #Load dg data
  dg_list <- list()
  for(protein in names(input_list)){
    temp_dt <- fread(input_list[[protein]])
    temp_dt[, protein := protein]
    #Add WT and mutant AAs
    temp_dt[, WT_AA := substr(id, 1, 1)]
    temp_dt[, Mut := substr(id, nchar(id), nchar(id))]
    #Per residue metrics
    for(i in c("f_ddg", "b_ddg")){
      temp_dt[,paste0(i, "_wposmean") := sum(.SD[[1]]/.SD[[2]]^2, na.rm = T)/sum(1/.SD[[2]]^2, na.rm = T),Pos_ref,.SDcols = paste0(i, c("_pred", "_pred_sd"))]
      temp_dt[,paste0(i, "_wposse") := sqrt(1/sum(1/.SD[[2]]^2, na.rm = T)),Pos_ref,.SDcols = paste0(i, c("_pred", "_pred_sd"))]
    }
    if(is.null(input_MSA_list[[protein]])){
      temp_dt[, conservation := NA]
    } else {
      MSA_dt <- fread(input_MSA_list[[protein]])
      temp_dt <- data.table::merge.data.table(temp_dt, MSA_dt[, .(i, conservation)], by.x = "Pos_ref", by.y="i", all.x = T)
    }
    dg_list[[protein]] <- temp_dt
  }
  dg_dt <- rbindlist(dg_list)

  ###########################
  ### Free energy distributions
  ###########################

  #Free energy distributions by protein - all - conf
  plot_dt <- data.table::melt.data.table(copy(dg_dt)[,.(protein, f_ddg_pred, b_ddg_pred, f_ddg_pred_conf, b_ddg_pred_conf)], id = c("protein", "f_ddg_pred_conf", "b_ddg_pred_conf"))
  plot_dt <- plot_dt[(variable=="f_ddg_pred" & f_ddg_pred_conf) | (variable=="b_ddg_pred" & b_ddg_pred_conf),]
  plot_dt[variable=="f_ddg_pred", variable_plot := "Folding"]
  plot_dt[variable=="b_ddg_pred", variable_plot := "Binding"]
  #Plot
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(value, fill = variable_plot)) +
    ggplot2::geom_density(alpha = 0.5) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2) +
    ggplot2::facet_wrap(~protein, scales = "free") +
    ggplot2::xlab(expression(Delta*Delta*"G")) +
    ggplot2::ylab("Density") +
    # ggplot2::coord_cartesian(xlim = c(-3, 6)) +
    ggplot2::labs(fill = "Trait") +
    ggplot2::theme_classic()
  if(!is.null(colour_scheme)){
    d <- d + ggplot2::scale_fill_manual(values = unlist(colour_scheme[["shade 0"]][c(1, 3)]))
  }
  ggplot2::ggsave(file.path(outpath, "ddG_densities_all.pdf"), d, width = 9, height = 3, useDingbats=FALSE)

  #Free energy distributions by protein - GRB2 and PSD95 - conf
  plot_dt <- data.table::melt.data.table(copy(dg_dt)[protein != "GB1",.(protein, f_ddg_pred, b_ddg_pred, f_ddg_pred_conf, b_ddg_pred_conf)], id = c("protein", "f_ddg_pred_conf", "b_ddg_pred_conf"))
  plot_dt <- plot_dt[(variable=="f_ddg_pred" & f_ddg_pred_conf) | (variable=="b_ddg_pred" & b_ddg_pred_conf),]
  plot_dt[variable=="f_ddg_pred", variable_plot := "Folding"]
  plot_dt[variable=="b_ddg_pred", variable_plot := "Binding"]
  #Plot
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(value, fill = variable_plot)) +
    ggplot2::geom_density(alpha = 0.5) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2) +
    ggplot2::facet_wrap(~protein, scales = "free") +
    ggplot2::xlab(expression(Delta*Delta*"G")) +
    ggplot2::ylab("Density") +
    # ggplot2::coord_cartesian(xlim = c(-3, 6)) +
    ggplot2::labs(fill = "Trait") +
    ggplot2::theme_classic()
  if(!is.null(colour_scheme)){
    d <- d + ggplot2::scale_fill_manual(values = unlist(colour_scheme[["shade 0"]][c(1, 3)]))
  }
  ggplot2::ggsave(file.path(outpath, "ddG_densities.pdf"), d, width = 6, height = 2, useDingbats=FALSE)

  #Free energy distributions by protein - GRB2 and PSD95 - conf - xlim
  plot_dt <- data.table::melt.data.table(copy(dg_dt)[protein != "GB1",.(protein, f_ddg_pred, b_ddg_pred, f_ddg_pred_conf, b_ddg_pred_conf)], id = c("protein", "f_ddg_pred_conf", "b_ddg_pred_conf"))
  plot_dt <- plot_dt[(variable=="f_ddg_pred" & f_ddg_pred_conf) | (variable=="b_ddg_pred" & b_ddg_pred_conf),]
  plot_dt[variable=="f_ddg_pred", variable_plot := "Folding"]
  plot_dt[variable=="b_ddg_pred", variable_plot := "Binding"]
  #Plot
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(value, fill = variable_plot)) +
    ggplot2::geom_density(alpha = 0.5) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2) +
    ggplot2::facet_wrap(~protein, scales = "free") +
    ggplot2::xlab(expression(Delta*Delta*"G")) +
    ggplot2::ylab("Density") +
    ggplot2::coord_cartesian(xlim = c(-1.5,3.5)) +
    ggplot2::labs(fill = "Trait") +
    ggplot2::theme_classic()
  if(!is.null(colour_scheme)){
    d <- d + ggplot2::scale_fill_manual(values = unlist(colour_scheme[["shade 0"]][c(1, 3)]))
  }
  ggplot2::ggsave(file.path(outpath, "ddG_densities_xlim.pdf"), d, width = 6, height = 2, useDingbats=FALSE)

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
    ggplot2::facet_wrap(~protein, scales = "free") +
    ggplot2::xlab(expression("Folding "*Delta*"G")) +
    ggplot2::ylab(expression("Binding "*Delta*"G")) +
    ggplot2::labs(color = "Residue\nposition") +
    ggplot2::theme_classic()
  if(!is.null(colour_scheme)){
    d <- d + ggplot2::scale_colour_manual(values = unlist(colour_scheme[["shade 0"]][c(1, 3, 4)]))
  }
  ggplot2::ggsave(file.path(outpath, "ddG_scatter.pdf"), d, width = 6, height = 3, useDingbats=FALSE)

  #Free energy scatterplots by protein - GRB2 and PSD95 - conf - xylim
  plot_dt <- copy(dg_dt)[protein != "GB1",.(protein, f_dg_pred, b_dg_pred, f_ddg_pred_conf, b_ddg_pred_conf, Pos_class, id)]
  plot_dt <- plot_dt[f_ddg_pred_conf & b_ddg_pred_conf,]
  plot_dt[, Pos_class_plot := stringr::str_to_title(Pos_class)]
  plot_dt[Pos_class=="binding_interface", Pos_class_plot := "Binding\ninterface"]
  #Plot
  d <- ggplot2::ggplot(plot_dt[id!="-0-" & protein=="PSD95-PDZ3"],ggplot2::aes(f_dg_pred, b_dg_pred, colour = Pos_class_plot)) +
    ggplot2::geom_point(alpha = 0.5, size = 1) +
    ggplot2::coord_cartesian(ylim = c(-2.5, 2.5), xlim = c(-2.5,2.5)) +
    ggplot2::geom_density_2d(contour_var = "ndensity", bins = 6) +
    ggplot2::geom_vline(xintercept = 0) +
    ggplot2::geom_vline(data = plot_dt[id=="-0-",], ggplot2::aes(xintercept = f_dg_pred), linetype = 2) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::geom_hline(data = plot_dt[id=="-0-",], ggplot2::aes(yintercept = b_dg_pred), linetype = 2) +
    ggplot2::facet_wrap(~protein, scales = "free") +
    ggplot2::xlab(expression("Folding "*Delta*"G")) +
    ggplot2::ylab(expression("Binding "*Delta*"G")) +
    ggplot2::labs(color = "Residue\nposition") +
    ggplot2::theme_classic()
  if(!is.null(colour_scheme)){
    d <- d + ggplot2::scale_colour_manual(values = unlist(colour_scheme[["shade 0"]][c(1, 3, 4)]))
  }
  ggplot2::ggsave(file.path(outpath, "ddG_scatter_contour_xylim.pdf"), d, width = 4, height = 3, useDingbats=FALSE)

  #Free energy scatterplots by protein - GB1 - conf - xylim
  plot_dt <- copy(dg_dt)[,.(protein, f_dg_pred, b_dg_pred, f_ddg_pred, b_ddg_pred, f_ddg_pred_conf, b_ddg_pred_conf, Pos_class, id)]
  plot_dt <- plot_dt[f_ddg_pred_conf & b_ddg_pred_conf,]
  plot_dt[, Pos_class_plot := stringr::str_to_title(Pos_class)]
  plot_dt[Pos_class=="binding_interface", Pos_class_plot := "Binding\ninterface"]
  #Plot
  d <- ggplot2::ggplot(plot_dt[id!="-0-" & protein=="GB1"],ggplot2::aes(f_dg_pred, b_dg_pred, colour = Pos_class_plot)) +
    ggplot2::geom_point(alpha = 0.5, size = 1) +
    ggplot2::geom_density_2d(contour_var = "ndensity", bins = 6) +
    ggplot2::geom_vline(xintercept = 0) +
    ggplot2::geom_vline(data = plot_dt[id=="-0-" & protein=="GB1",], ggplot2::aes(xintercept = f_dg_pred), linetype = 2) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::geom_hline(data = plot_dt[id=="-0-" & protein=="GB1",], ggplot2::aes(yintercept = b_dg_pred), linetype = 2) +
    ggplot2::xlab(expression("Folding "*Delta*"G")) +
    ggplot2::ylab(expression("Binding "*Delta*"G")) +
    ggplot2::labs(color = "Residue\nposition") +
    ggplot2::coord_cartesian(ylim = c(-1.5, 2), xlim = c(-6,1)) +
    ggplot2::theme_classic()
  if(!is.null(colour_scheme)){
    d <- d + ggplot2::scale_colour_manual(values = unlist(colour_scheme[["shade 0"]][c(1, 3, 4)]))
  }
  ggplot2::ggsave(file.path(outpath, "dG_scatter_contour_GB1_xylim.pdf"), d, width = 4, height = 3, useDingbats=FALSE)

  #Free energy change scatterplots by protein - GB1 - conf - xylim
  #Plot
  d <- ggplot2::ggplot(plot_dt[id!="-0-" & protein=="GB1"],ggplot2::aes(f_ddg_pred, b_ddg_pred, colour = Pos_class_plot)) +
    ggplot2::geom_point(alpha = 0.5, size = 1) +
    ggplot2::geom_density_2d(contour_var = "ndensity", bins = 6) +
    ggplot2::geom_vline(xintercept = 0) +
    # ggplot2::geom_vline(data = plot_dt[id=="-0-" & protein=="GB1",], ggplot2::aes(xintercept = f_dg_pred), linetype = 2) +
    ggplot2::geom_hline(yintercept = 0) +
    # ggplot2::geom_hline(data = plot_dt[id=="-0-" & protein=="GB1",], ggplot2::aes(yintercept = b_dg_pred), linetype = 2) +
    ggplot2::xlab(expression("Folding "*Delta*Delta*"G")) +
    ggplot2::ylab(expression("Binding "*Delta*Delta*"G")) +
    ggplot2::labs(color = "Residue\nposition") +
    ggplot2::coord_cartesian(ylim = c(-1, 2.5), xlim = c(-2,5)) +
    ggplot2::theme_classic()
  if(!is.null(colour_scheme)){
    d <- d + ggplot2::scale_colour_manual(values = unlist(colour_scheme[["shade 0"]][c(1, 3, 4)]))
  }
  ggplot2::ggsave(file.path(outpath, "ddG_scatter_contour_GB1_xylim.pdf"), d, width = 4, height = 3, useDingbats=FALSE)

  #Free energy scatterplots by protein - GRB2 - conf - xylim
  plot_dt <- copy(dg_dt)[protein != "GB1",.(protein, f_dg_pred, b_dg_pred, f_ddg_pred, b_ddg_pred, f_ddg_pred_conf, b_ddg_pred_conf, Pos_class, id)]
  plot_dt <- plot_dt[f_ddg_pred_conf & b_ddg_pred_conf,]
  plot_dt[, Pos_class_plot := stringr::str_to_title(Pos_class)]
  plot_dt[Pos_class=="binding_interface", Pos_class_plot := "Binding\ninterface"]
  #Plot
  d <- ggplot2::ggplot(plot_dt[id!="-0-" & protein=="GRB2-SH3"],ggplot2::aes(f_dg_pred, b_dg_pred, colour = Pos_class_plot)) +
    ggplot2::geom_point(alpha = 0.5, size = 1) +
    ggplot2::geom_density_2d(contour_var = "ndensity", bins = 6) +
    ggplot2::geom_vline(xintercept = 0) +
    ggplot2::geom_vline(data = plot_dt[id=="-0-" & protein=="GRB2-SH3",], ggplot2::aes(xintercept = f_dg_pred), linetype = 2) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::geom_hline(data = plot_dt[id=="-0-" & protein=="GRB2-SH3",], ggplot2::aes(yintercept = b_dg_pred), linetype = 2) +
    ggplot2::xlab(expression("Folding "*Delta*"G")) +
    ggplot2::ylab(expression("Binding "*Delta*"G")) +
    ggplot2::labs(color = "Residue\nposition") +
    ggplot2::coord_cartesian(ylim = c(-2.5, 2), xlim = c(-2,2.5)) +
    ggplot2::theme_classic()
  if(!is.null(colour_scheme)){
    d <- d + ggplot2::scale_colour_manual(values = unlist(colour_scheme[["shade 0"]][c(1, 3, 4)]))
  }
  ggplot2::ggsave(file.path(outpath, "dG_scatter_contour_GRB2-SH3_xylim.pdf"), d, width = 4, height = 3, useDingbats=FALSE)

  #Free energy scatterplots by protein - GRB2 - conf - xylim
  #Plot
  d <- ggplot2::ggplot(plot_dt[id!="-0-" & protein=="GRB2-SH3"],ggplot2::aes(f_ddg_pred, b_ddg_pred, colour = Pos_class_plot)) +
    ggplot2::geom_point(alpha = 0.5, size = 1) +
    ggplot2::geom_density_2d(contour_var = "ndensity", bins = 6) +
    ggplot2::geom_vline(xintercept = 0) +
    # ggplot2::geom_vline(data = plot_dt[id=="-0-" & protein=="GRB2-SH3",], ggplot2::aes(xintercept = f_dg_pred), linetype = 2) +
    ggplot2::geom_hline(yintercept = 0) +
    # ggplot2::geom_hline(data = plot_dt[id=="-0-" & protein=="GRB2-SH3",], ggplot2::aes(yintercept = b_dg_pred), linetype = 2) +
    ggplot2::xlab(expression("Folding "*Delta*Delta*"G")) +
    ggplot2::ylab(expression("Binding "*Delta*Delta*"G")) +
    ggplot2::labs(color = "Residue\nposition") +
    ggplot2::coord_cartesian(ylim = c(-1.5, 3), xlim = c(-1,3.5)) +
    ggplot2::theme_classic()
  if(!is.null(colour_scheme)){
    d <- d + ggplot2::scale_colour_manual(values = unlist(colour_scheme[["shade 0"]][c(1, 3, 4)]))
  }
  ggplot2::ggsave(file.path(outpath, "ddG_scatter_contour_GRB2-SH3_xylim.pdf"), d, width = 4, height = 3, useDingbats=FALSE)

  #Free energy scatterplots by protein - PSD95 - conf - xylim
  plot_dt <- copy(dg_dt)[protein != "GB1",.(protein, f_dg_pred, b_dg_pred, f_ddg_pred, b_ddg_pred, f_ddg_pred_conf, b_ddg_pred_conf, Pos_class, id)]
  plot_dt <- plot_dt[f_ddg_pred_conf & b_ddg_pred_conf,]
  plot_dt[, Pos_class_plot := stringr::str_to_title(Pos_class)]
  plot_dt[Pos_class=="binding_interface", Pos_class_plot := "Binding\ninterface"]
  #Plot
  d <- ggplot2::ggplot(plot_dt[id!="-0-" & protein=="PSD95-PDZ3"],ggplot2::aes(f_dg_pred, b_dg_pred, colour = Pos_class_plot)) +
    ggplot2::geom_point(alpha = 0.5, size = 1) +
    ggplot2::geom_density_2d(contour_var = "ndensity", bins = 6) +
    ggplot2::geom_vline(xintercept = 0) +
    ggplot2::geom_vline(data = plot_dt[id=="-0-" & protein=="PSD95-PDZ3",], ggplot2::aes(xintercept = f_dg_pred), linetype = 2) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::geom_hline(data = plot_dt[id=="-0-" & protein=="PSD95-PDZ3",], ggplot2::aes(yintercept = b_dg_pred), linetype = 2) +
    ggplot2::xlab(expression("Folding "*Delta*"G")) +
    ggplot2::ylab(expression("Binding "*Delta*"G")) +
    ggplot2::labs(color = "Residue\nposition") +
    ggplot2::coord_cartesian(ylim = c(-1.5, 1.2), xlim = c(-2.5,1.5)) +
    ggplot2::theme_classic()
  if(!is.null(colour_scheme)){
    d <- d + ggplot2::scale_colour_manual(values = unlist(colour_scheme[["shade 0"]][c(1, 3, 4)]))
  }
  ggplot2::ggsave(file.path(outpath, "dG_scatter_contour_PSD95-PDZ3_xylim.pdf"), d, width = 4, height = 3, useDingbats=FALSE)

  #Free energy scatterplots by protein - PSD95 - conf - xylim
  #Plot
  d <- ggplot2::ggplot(plot_dt[id!="-0-" & protein=="PSD95-PDZ3"],ggplot2::aes(f_ddg_pred, b_ddg_pred, colour = Pos_class_plot)) +
    ggplot2::geom_point(alpha = 0.5, size = 1) +
    ggplot2::geom_density_2d(contour_var = "ndensity", bins = 6) +
    ggplot2::geom_vline(xintercept = 0) +
    # ggplot2::geom_vline(data = plot_dt[id=="-0-" & protein=="PSD95-PDZ3",], ggplot2::aes(xintercept = f_dg_pred), linetype = 2) +
    ggplot2::geom_hline(yintercept = 0) +
    # ggplot2::geom_hline(data = plot_dt[id=="-0-" & protein=="PSD95-PDZ3",], ggplot2::aes(yintercept = b_dg_pred), linetype = 2) +
    ggplot2::xlab(expression("Folding "*Delta*Delta*"G")) +
    ggplot2::ylab(expression("Binding "*Delta*Delta*"G")) +
    ggplot2::labs(color = "Residue\nposition") +
    ggplot2::coord_cartesian(ylim = c(-1, 1.7), xlim = c(-1,3)) +
    ggplot2::theme_classic()
  if(!is.null(colour_scheme)){
    d <- d + ggplot2::scale_colour_manual(values = unlist(colour_scheme[["shade 0"]][c(1, 3, 4)]))
  }
  ggplot2::ggsave(file.path(outpath, "ddG_scatter_contour_PSD95-PDZ3_xylim.pdf"), d, width = 4, height = 3, useDingbats=FALSE)

  ###########################
  ### Free energy vs conservation scatterplots
  ###########################
  
  proteins_MSA <- names(input_MSA_list)[(do.call("c", lapply(input_MSA_list, function(x){!is.null(x)})))]
  
  plot_dt <- dg_dt[id!="-0-" & protein %in% proteins_MSA][!duplicated(paste(Pos_ref, protein, sep = ":"))]
  cor_dt <- plot_dt[,.(cor = round(cor(conservation, f_ddg_wposmean, use = "pairwise.complete"), 2), sig = format(cor.test(conservation, f_ddg_wposmean, use = "pairwise.complete")$p.value), scientific = T, digits = 2),.(protein, Pos_class)]
  cor_dt[, conservation := 0.5]
  cor_dt[, f_ddg_wposmean := seq(0, plot_dt[,max(f_ddg_wposmean)], plot_dt[,max(f_ddg_wposmean)]/(3+1))[2:(3+1)], protein]
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(conservation, f_ddg_wposmean)) +
    ggplot2::geom_point(ggplot2::aes(color = Pos_class)) +
    ggplot2::geom_linerange(ggplot2::aes(ymin = f_ddg_wposmean-f_ddg_wposse*1.96, ymax = f_ddg_wposmean+f_ddg_wposse*1.96, color = Pos_class)) +
    ggplot2::geom_smooth(method = "lm", formula = 'y~x', ggplot2::aes(color = Pos_class), se = F, show.legend = F) +
    ggplot2::xlab("Residue conservation") +
    ggplot2::ylab(expression("Weighted mean |Folding "*Delta*Delta*"G|")) +
    ggplot2::geom_text(data = cor_dt, ggplot2::aes(label=paste("r = ", cor, ", P = ", sig, sep=""), color = Pos_class)) +
    ggplot2::facet_grid(~protein, scales = "free") + 
    ggplot2::theme_bw()
  if(!is.null(colour_scheme)){
    d <- d + ggplot2::scale_color_manual(values = unlist(colour_scheme[["shade 0"]][c(1, 3, 4)]))
  }
  suppressWarnings(ggplot2::ggsave(file.path(outpath, "mean_ddG_folding_conservation_scatter.pdf"), d, width = 7, height = 3))
  
  cor_dt <- plot_dt[,.(cor = round(cor(conservation, b_ddg_wposmean, use = "pairwise.complete"), 2), sig = format(cor.test(conservation, b_ddg_wposmean, use = "pairwise.complete")$p.value), scientific = T, digits = 2),.(protein, Pos_class)]
  cor_dt[, conservation := 0.5]
  cor_dt[, b_ddg_wposmean := seq(0, plot_dt[,max(b_ddg_wposmean)], plot_dt[,max(b_ddg_wposmean)]/(3+1))[2:(3+1)], protein]
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(conservation, b_ddg_wposmean)) +
    ggplot2::geom_point(ggplot2::aes(color = Pos_class)) +
    ggplot2::geom_linerange(ggplot2::aes(ymin = b_ddg_wposmean-b_ddg_wposse*1.96, ymax = b_ddg_wposmean+f_ddg_wposse*1.96, color = Pos_class)) +
    ggplot2::geom_smooth(method = "lm", formula = 'y~x', ggplot2::aes(color = Pos_class), se = F, show.legend = F) +
    ggplot2::xlab("Residue conservation") +
    ggplot2::ylab(expression("Weighted mean |Binding "*Delta*Delta*"G|")) +
    ggplot2::geom_text(data = cor_dt, ggplot2::aes(label=paste("r = ", cor, ", P = ", sig, sep=""), color = Pos_class)) +
    ggplot2::facet_grid(~protein, scales = "free") + 
    ggplot2::theme_bw()
  if(!is.null(colour_scheme)){
    d <- d + ggplot2::scale_color_manual(values = unlist(colour_scheme[["shade 0"]][c(1, 3, 4)]))
  }
  suppressWarnings(ggplot2::ggsave(file.path(outpath, "mean_ddG_binding_conservation_scatter.pdf"), d, width = 7, height = 3))
  
}

