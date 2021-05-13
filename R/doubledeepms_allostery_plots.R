
#' doubledeepms_allostery_plots
#'
#' Plot free energy heatmaps.
#'
#' @param input_list path to MoCHI thermo model fit results (required)
#' @param aaprop_file path to amino acid properties file (required)
#' @param aaprop_file_selected path to file with selected subset of identifiers
#' @param outpath output path for plots and saved objects (required)
#' @param colour_scheme colour scheme file (required)
#' @param execute whether or not to execute the analysis (default: TRUE)
#'
#' @return Nothing
#' @export
#' @import data.table
doubledeepms_allostery_plots <- function(
  input_list,
  aaprop_file,
  aaprop_file_selected,
  outpath,
  colour_scheme,
  execute = TRUE
  ){

  #Return if analysis not executed
  if(!execute){
    return()
  }

  #Display status
  message(paste("\n\n*******", "running stage: doubledeepms_allostery_plots", "*******\n\n"))

  #Create output directory
  doubledeepms__create_dir(doubledeepms_dir = outpath)

  ### Load single mutant free energies
  ###########################

  #Load dg data
  dg_list <- list()
  for(protein in names(input_list)){
    temp_dt <- fread(input_list[[protein]])[id!="-0-"]
    temp_dt[, protein := protein]

    #Add WT and mutant AAs
    temp_dt[, WT_AA := substr(id, 1, 1)]
    temp_dt[, Mut := substr(id, nchar(id), nchar(id))]

    #Per residue metrics
    for(i in c("f_ddg", "b_ddg")){
      temp_dt[,paste0(i, "_posmeanabs") := mean(abs(.SD[[1]]), na.rm = T),Pos_ref,.SDcols = paste0(i, c("_pred"))]
      temp_dt[,paste0(i, "_posse") := sd(abs(.SD[[1]]), na.rm = T)/sqrt(sum(!is.na(.SD[[1]]))),Pos_ref,.SDcols = paste0(i, c("_pred"))]
      temp_dt[,paste0(i, "_wposmeanabs") := sum(abs(.SD[[1]])/.SD[[2]]^2, na.rm = T)/sum(1/.SD[[2]]^2, na.rm = T),Pos_ref,.SDcols = paste0(i, c("_pred", "_pred_sd"))]
      temp_dt[,paste0(i, "_wposse") := sqrt(1/sum(1/.SD[[2]]^2, na.rm = T)),Pos_ref,.SDcols = paste0(i, c("_pred", "_pred_sd"))]
    }

    #Per residue metrics - confident only
    for(i in c("f_ddg", "b_ddg")){
      temp_dt[get(paste0(i, "_pred_conf"))==TRUE,paste0(i, "_pred_filtered") := .SD[[1]],,.SDcols = paste0(i, "_pred")]
      temp_dt[get(paste0(i, "_pred_conf"))==TRUE,paste0(i, "_pred_sd_filtered") := .SD[[1]],,.SDcols = paste0(i, "_pred_sd")]
      temp_dt[,paste0(i, "_posmeanabs_conf") := mean(abs(.SD[[1]]), na.rm = T),Pos_ref,.SDcols = paste0(i, c("_pred_filtered"))]
      temp_dt[,paste0(i, "_posse_conf") := sd(abs(.SD[[1]]), na.rm = T)/sqrt(sum(!is.na(.SD[[1]]))),Pos_ref,.SDcols = paste0(i, c("_pred_filtered"))]
      temp_dt[,paste0(i, "_wposmeanabs_conf") := sum(abs(.SD[[1]])/.SD[[2]]^2, na.rm = T)/sum(1/.SD[[2]]^2, na.rm = T),Pos_ref,.SDcols = paste0(i, c("_pred_filtered", "_pred_sd_filtered"))]
      temp_dt[,paste0(i, "_wposse_conf") := sqrt(1/sum(1/.SD[[2]]^2, na.rm = T)),Pos_ref,.SDcols = paste0(i, c("_pred_filtered", "_pred_sd_filtered"))]
    }
    dg_list[[protein]] <- temp_dt
  }
  dg_dt <- rbindlist(dg_list)

  ###########################
  ### Folding energy distance correlation plots
  ###########################

  for(i in dg_dt[,unique(protein)]){
    doubledeepms__persite_energy_vs_distance_plot(
      input_dt = copy(dg_dt)[protein==i],
      outpath = file.path(outpath, paste0(i, "_persite_folding_energy_vs_distance_scatter.pdf")),
      colour_scheme = colour_scheme,
      trait_name = "folding")
  }

  ###########################
  ### Binding energy distance correlation plots
  ###########################

  for(i in dg_dt[,unique(protein)]){
    doubledeepms__persite_energy_vs_distance_plot(
      input_dt = copy(dg_dt)[protein==i],
      outpath = file.path(outpath, paste0(i, "_persite_binding_energy_vs_distance_scatter.pdf")),
      colour_scheme = colour_scheme,
      trait_name = "binding")
  }

  # ###########################
  # ### Binding ROC plot
  # ###########################

  # metric_names <- c("b_ddg_wposmeanabs", "b_ddg_wposmeanabs_conf", "b_ddg_posmeanabs", "b_ddg_posmeanabs_conf")
  # metric_names_plot <- c(
  #   "Weighted mean |Binding ddG|",
  #   "Weighted mean |Binding ddG| (conf.)",
  #   "Mean |Binding ddG|",
  #   "Mean |Binding ddG| (conf.)")
  # names(metric_names_plot) <- metric_names

  # #Subset
  # subset_dt <- dg_dt[order(Pos_ref)][!duplicated(Pos_ref)]

  # #Performance of all metrics
  # perf_list <- list()
  # for(metric_name in c("b_ddg_wposmeanabs", "b_ddg_wposmeanabs_conf", "b_ddg_posmeanabs", "b_ddg_posmeanabs_conf")){
  #   #Metric
  #   subset_dt[, plot_metric := .SD[[1]],,.SDcols = metric_name]
  #   #ROC curve data
  #   roc_df <- data.frame(
  #     predictions = subset_dt[,plot_metric], 
  #     labels = subset_dt[,as.numeric(Pos_class=="binding_interface")])
  #   pred <- ROCR::prediction(roc_df$predictions, roc_df$labels)
  #   perf <- ROCR::performance(pred,"tpr","fpr")
  #   auc <- round(ROCR::performance(pred, measure = "auc")@'y.values'[[1]], 2)
  #   #Save
  #   perf_list[[metric_name]] <- data.table(
  #     FPR = perf@'x.values'[[1]],
  #     TPR = perf@'y.values'[[1]],
  #     measure = metric_name,
  #     auc = auc)
  # }
  # plot_dt <- rbindlist(perf_list)
  # plot_dt[, measure_plot := metric_names_plot[measure]]
  # plot_cols <- c(colour_scheme[["shade 0"]][[1]], colour_scheme[["shade 0"]][[2]], colour_scheme[["shade 0"]][[3]], colour_scheme[["shade 0"]][[4]])
  # names(plot_cols) <- metric_names_plot

  # #Plot
  # auc_dt <- plot_dt[!duplicated(measure_plot)][order(measure_plot, decreasing = T)]
  # auc_dt[, FPR := 0.5]
  # auc_dt[, TPR := seq(0, 1, 1/(.N+1))[2:(.N+1)]]
  # d <- ggplot2::ggplot(plot_dt,ggplot2::aes(FPR, TPR, color = measure_plot)) +
  #   ggplot2::geom_line() +
  #   ggplot2::geom_abline(linetype = 2) +
  #   ggplot2::xlab("FPR") +
  #   ggplot2::ylab("TPR") +
  #   ggplot2::geom_text(data = auc_dt, ggplot2::aes(label=paste("AUC = ", auc, sep=""))) +
  #   ggplot2::theme_bw() +
  #   ggplot2::scale_colour_manual(values=plot_cols) +
  #   ggplot2::labs(color = "Binding interface\nprediction metric")   
  # ggplot2::ggsave(file.path(outpath, "binding_ROC.pdf"), d, width = 6, height = 3, useDingbats=FALSE)

}

