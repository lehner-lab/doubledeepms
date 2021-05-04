
#' doubledeepms_free_energy_heatmaps
#'
#' Plot free energy heatmaps.
#'
#' @param input_file path to MoCHI thermo model fit results (required)
#' @param domain_name domain name (required)
#' @param outpath output path for plots and saved objects (required)
#' @param colour_scheme colour scheme file (required)
#' @param plot_width heatmap plot width in inches (default:10)
#' @param plot_height heatmap plot height in inches (default:4)
#' @param execute whether or not to execute the analysis (default: TRUE)
#'
#' @return Nothing
#' @export
#' @import data.table
doubledeepms_free_energy_heatmaps <- function(
  input_file,
  domain_name,
  outpath,
  colour_scheme,
  plot_width = 10,
  plot_height = 4,
  execute = TRUE
  ){

  #Return if analysis not executed
  if(!execute){
    return()
  }

  #Display status
  message(paste("\n\n*******", paste("running stage: doubledeepms_free_energy_heatmaps for", domain_name), "*******\n\n"))

  #Create output directory
  doubledeepms__create_dir(doubledeepms_dir = outpath)

  #Load free energies
  dg_dt <- fread(input_file)[id!="-0-"]

  #Add WT and mutant AAs
  dg_dt[, WT_AA := substr(id, 1, 1)]
  dg_dt[, Mut := substr(id, nchar(id), nchar(id))]

  #Per residue metrics
  for(i in c("f_ddg", "b_ddg")){
    dg_dt[,paste0(i, "_posmeanabs") := mean(abs(.SD[[1]]), na.rm = T),Pos_ref,.SDcols = paste0(i, c("_pred"))]
    dg_dt[,paste0(i, "_posse") := sd(abs(.SD[[1]]), na.rm = T)/sqrt(sum(!is.na(.SD[[1]]))),Pos_ref,.SDcols = paste0(i, c("_pred"))]
    dg_dt[,paste0(i, "_wposmeanabs") := sum(abs(.SD[[1]])/.SD[[2]]^2, na.rm = T)/sum(1/.SD[[2]]^2, na.rm = T),Pos_ref,.SDcols = paste0(i, c("_pred", "_pred_sd"))]
    dg_dt[,paste0(i, "_wposse") := sqrt(1/sum(1/.SD[[2]]^2, na.rm = T)),Pos_ref,.SDcols = paste0(i, c("_pred", "_pred_sd"))]
  }

  #Per residue metrics - confident only
  for(i in c("f_ddg", "b_ddg")){
    dg_dt[get(paste0(i, "_pred_conf"))==TRUE,paste0(i, "_pred_filtered") := .SD[[1]],,.SDcols = paste0(i, "_pred")]
    dg_dt[get(paste0(i, "_pred_conf"))==TRUE,paste0(i, "_pred_sd_filtered") := .SD[[1]],,.SDcols = paste0(i, "_pred_sd")]
    dg_dt[,paste0(i, "_posmeanabs_conf") := mean(abs(.SD[[1]]), na.rm = T),Pos_ref,.SDcols = paste0(i, c("_pred_filtered"))]
    dg_dt[,paste0(i, "_posse_conf") := sd(abs(.SD[[1]]), na.rm = T)/sqrt(sum(!is.na(.SD[[1]]))),Pos_ref,.SDcols = paste0(i, c("_pred_filtered"))]
    dg_dt[,paste0(i, "_wposmeanabs_conf") := sum(abs(.SD[[1]])/.SD[[2]]^2, na.rm = T)/sum(1/.SD[[2]]^2, na.rm = T),Pos_ref,.SDcols = paste0(i, c("_pred_filtered", "_pred_sd_filtered"))]
    dg_dt[,paste0(i, "_wposse_conf") := sqrt(1/sum(1/.SD[[2]]^2, na.rm = T)),Pos_ref,.SDcols = paste0(i, c("_pred_filtered", "_pred_sd_filtered"))]
  }

  ###########################
  ### Folding heatmap
  ###########################

  doubledeepms__plot_heatmap(
    input_dt = dg_dt,
    variable_name = "f_ddg_pred",
    metric_name = "f_ddg_posmeanabs_conf",
    output_file = file.path(outpath, "folding_heatmap.pdf"),
    width = plot_width,
    height = plot_height,
    plot_title = paste0(domain_name, " amino acid position"),
    colour_clip = 2.5,
    colour_low = colour_scheme[["shade 0"]][[3]],
    colour_high = colour_scheme[["shade 0"]][[1]])

  ###########################
  ### Binding heatmap
  ###########################

  doubledeepms__plot_heatmap(
    input_dt = dg_dt,
    variable_name = "b_ddg_pred",
    metric_name = "b_ddg_posmeanabs_conf",
    output_file = file.path(outpath, "binding_heatmap.pdf"),
    width = plot_width,
    height = plot_height,
    plot_title = paste0(domain_name, " amino acid position"),
    colour_clip = 2.5,
    colour_low = colour_scheme[["shade 0"]][[3]],
    colour_high = colour_scheme[["shade 0"]][[1]])

  ###########################
  ### Folding correlation plot
  ###########################

  #Metrics
  metric_name <- "f_ddg_posmeanabs_conf"
  metric_name_plot <- expression("Mean |Folding "*Delta*Delta*"G|")
  metric_name_se <- "f_ddg_posse_conf"

  #Subset
  plot_dt <- dg_dt[order(Pos_ref)][!duplicated(Pos_ref)]

  #Metric
  plot_dt[, plot_metric := .SD[[1]],,.SDcols = metric_name]
  plot_dt[, plot_metric_se := .SD[[1]],,.SDcols = metric_name_se]

  #Allostery linear model
  plot_dt[, class := "Remainder"]
  plot_dt[Pos_class=="binding_interface", class := "Binding\ninterface"]
  plot_dt[, class := factor(class, levels = c("Distal", "Binding\ninterface", "Remainder"))]
  plot_cols <- c(colour_scheme[["shade 0"]][[1]], colour_scheme[["shade 0"]][[4]], colour_scheme[["shade 0"]][[3]])
  names(plot_cols) <- c("Distal", "Binding\ninterface", "Remainder")

  #Plot
  d <- ggplot2::ggplot(plot_dt[!is.na(plot_metric)],ggplot2::aes(scHAmin_ligand, plot_metric, colour = class)) +
    ggplot2::geom_vline(xintercept = 5, linetype = 2) +
    ggplot2::geom_point(alpha = 3/4, size = 3) +
    # ggplot2::geom_linerange(ggplot2::aes(ymin = plot_metric-plot_metric_se*1.96, ymax = plot_metric+plot_metric_se*1.96)) +
    ggplot2::xlab(expression("Distance to ligand ("*ring(A)*")")) +
    ggplot2::ylab(metric_name_plot) +
    ggplot2::geom_text(ggplot2::aes(label = Pos_ref), colour = "black", size = 1) +
    # ggplot2::annotate("text", label=paste("Pearson's r = ", plot_dt[,round(cor(scHAmin_ligand, plot_metric, use = "pairwise.complete"), 2)], sep=""), x=-Inf, y=Inf, hjust = 0, vjust = 1) +
    ggplot2::theme_bw() +
    ggplot2::scale_colour_manual(values=plot_cols) +
    ggplot2::labs(color = "Residue\nposition")   
  ggplot2::ggsave(file.path(outpath, "persite_folding_energy_vs_distance_scatter.pdf"), d, width = 4, height = 3, useDingbats=FALSE)

  ###########################
  ### Binding correlation plot
  ###########################

  #Metrics
  metric_name <- "b_ddg_posmeanabs_conf"
  metric_name_plot <- expression("Mean |Binding "*Delta*Delta*"G|")
  metric_name_se <- "b_ddg_posse_conf"

  #Subset
  plot_dt <- dg_dt[order(Pos_ref)][!duplicated(Pos_ref)]

  #Metric
  plot_dt[, plot_metric := .SD[[1]],,.SDcols = metric_name]
  plot_dt[, plot_metric_se := .SD[[1]],,.SDcols = metric_name_se]

  #Allostery linear model
  plot_dt[, class := "remainder"]
  plot_dt[Pos_class=="binding_interface", class := "bind"]

  #Threshold for binding interface modulating residues
  reg_threshold <- plot_dt[Pos_class=="binding_interface"][order(plot_metric, decreasing = T)][5][,plot_metric]

  #Allosteric sites
  allostery_pos <- plot_dt[plot_metric>=reg_threshold,Pos_ref]
  plot_dt[, class := "Remainder"]
  plot_dt[Pos_ref %in% allostery_pos, class := "Distal"]
  plot_dt[Pos_class=="binding_interface", class := "Binding\ninterface"]
  plot_dt[, class := factor(class, levels = c("Distal", "Binding\ninterface", "Remainder"))]
  plot_cols <- c(colour_scheme[["shade 0"]][[1]], colour_scheme[["shade 0"]][[4]], colour_scheme[["shade 0"]][[3]])
  names(plot_cols) <- c("Distal", "Binding\ninterface", "Remainder")

  #Plot
  d <- ggplot2::ggplot(plot_dt[!is.na(plot_metric)],ggplot2::aes(scHAmin_ligand, plot_metric)) +
    ggplot2::geom_smooth(color = "grey", method = "glm", formula = y ~ exp(-x)) +
    ggplot2::geom_vline(xintercept = 5, linetype = 2) +
    ggplot2::geom_hline(yintercept = reg_threshold, linetype = 2) +
    ggplot2::geom_point(alpha = 3/4, size = 2.5, ggplot2::aes(colour = class)) +
    # ggplot2::geom_linerange(ggplot2::aes(ymin = plot_metric-plot_metric_se*1.96, ymax = plot_metric+plot_metric_se*1.96)) +
    ggplot2::xlab(expression("Distance to ligand ("*ring(A)*")")) +
    ggplot2::ylab(metric_name_plot) +
    ggplot2::geom_text(ggplot2::aes(label = Pos_ref), colour = "black", size = 1) +
    ggplot2::annotate("text", label=paste("R-squared = ", plot_dt[,round(cor(exp(-scHAmin_ligand), plot_metric, use = "pairwise.complete")^2, 2)], sep=""), x=Inf, y=Inf, hjust = 1, vjust = 1) +
      # "R-squared = ", plot_dt[,round(cor(scHAmin_ligand, plot_metric, use = "pairwise.complete")^2, 2)], sep=""), x=Inf, y=Inf, hjust = 1, vjust = 1) +
    ggplot2::theme_bw() + 
    ggplot2::scale_colour_manual(values=plot_cols) +
    ggplot2::labs(color = "Residue\nposition")   
  ggplot2::ggsave(file.path(outpath, "persite_binding_energy_vs_distance_scatter.pdf"), d, width = 4, height = 3, useDingbats=FALSE)

  ###########################
  ### Binding ROC plot
  ###########################

  metric_names <- c("b_ddg_wposmeanabs", "b_ddg_wposmeanabs_conf", "b_ddg_posmeanabs", "b_ddg_posmeanabs_conf")
  metric_names_plot <- c(
    "Weighted mean |Binding ddG|",
    "Weighted mean |Binding ddG| (conf.)",
    "Mean |Binding ddG|",
    "Mean |Binding ddG| (conf.)")
  names(metric_names_plot) <- metric_names

  #Subset
  subset_dt <- dg_dt[order(Pos_ref)][!duplicated(Pos_ref)]

  #Performance of all metrics
  perf_list <- list()
  for(metric_name in c("b_ddg_wposmeanabs", "b_ddg_wposmeanabs_conf", "b_ddg_posmeanabs", "b_ddg_posmeanabs_conf")){
    #Metric
    subset_dt[, plot_metric := .SD[[1]],,.SDcols = metric_name]
    #ROC curve data
    roc_df <- data.frame(
      predictions = subset_dt[,plot_metric], 
      labels = subset_dt[,as.numeric(Pos_class=="binding_interface")])
    pred <- ROCR::prediction(roc_df$predictions, roc_df$labels)
    perf <- ROCR::performance(pred,"tpr","fpr")
    auc <- round(ROCR::performance(pred, measure = "auc")@'y.values'[[1]], 2)
    #Save
    perf_list[[metric_name]] <- data.table(
      FPR = perf@'x.values'[[1]],
      TPR = perf@'y.values'[[1]],
      measure = metric_name,
      auc = auc)
  }
  plot_dt <- rbindlist(perf_list)
  plot_dt[, measure_plot := metric_names_plot[measure]]
  plot_cols <- c(colour_scheme[["shade 0"]][[1]], colour_scheme[["shade 0"]][[2]], colour_scheme[["shade 0"]][[3]], colour_scheme[["shade 0"]][[4]])
  names(plot_cols) <- metric_names_plot

  #Plot
  auc_dt <- plot_dt[!duplicated(measure_plot)][order(measure_plot, decreasing = T)]
  auc_dt[, FPR := 0.5]
  auc_dt[, TPR := seq(0, 1, 1/(.N+1))[2:(.N+1)]]
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(FPR, TPR, color = measure_plot)) +
    ggplot2::geom_line() +
    ggplot2::xlab("FPR") +
    ggplot2::ylab("TPR") +
    ggplot2::geom_text(data = auc_dt, ggplot2::aes(label=paste("AUC = ", auc, sep=""))) +
    ggplot2::theme_bw() +
    ggplot2::scale_colour_manual(values=plot_cols) +
    ggplot2::labs(color = "Binding interface\nprediction metric")   
  ggplot2::ggsave(file.path(outpath, "binding_ROC.pdf"), d, width = 6, height = 3, useDingbats=FALSE)

}

