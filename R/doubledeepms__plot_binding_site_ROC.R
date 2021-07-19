
#' doubledeepms__plot_binding_site_ROC
#'
#' Plot binding site prediction ROC plot.
#'
#' @param input_dt input data table (required)
#' @param outpath plot output path (required)
#' @param colour_scheme colour scheme file (required)
#' @param metric_names vector of variable names of the input data with the metric to compute the ROC curve 
#' @param metric_names_plot vector of characters with the names of each metric to appear in the ROC plot 
#' 
#'
#' @return Nothing
#' @export
doubledeepms__plot_binding_site_ROC <- function(
  input_dt, 
  outpath,
  colour_scheme,
  metric_names,
  metric_names_plot
  ){

  #Metric names
  if(is.null(metric_names)){
    metric_names <- c(
    # "b_ddg_posmaxabs", 
    # "b_ddg_posmaxabs_conf",
    # "effect_size_mean",
    "b_ddg_posmeanabs", 
    "b_ddg_posmeanabs_conf",
    "b_ddg_wposmeanabs", 
    "b_ddg_wposmeanabs_conf")
    }
  if(is.null(metric_names_plot)){
    metric_names_plot <- c(
      # "Max |Binding ddG|",
      # "Max |Binding ddG| (conf.)",
      # "Effect size |Binding ddG| (conf.)",
      "Mean |Binding ddG|",
      "Mean |Binding ddG| (conf.)",
      "Weighted mean |Binding ddG|",
      "Weighted mean |Binding ddG| (conf.)")
    }
  names(metric_names_plot) <- metric_names
  
  #Subset
  subset_dt <- input_dt[order(Pos_ref)][!duplicated(Pos_ref)]

  #Performance of all metrics
  perf_list <- list()
  for(metric_name in metric_names){
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
    ggplot2::geom_abline(linetype = 2) +
    ggplot2::xlab("FPR") +
    ggplot2::ylab("TPR") +
    ggplot2::geom_text(data = auc_dt, ggplot2::aes(label=paste("AUC = ", auc, sep=""))) +
    ggplot2::theme_bw() +
    ggplot2::scale_colour_manual(values=plot_cols) +
    ggplot2::labs(color = "Binding interface\nprediction metric")   
  ggplot2::ggsave(outpath, d, width = 6, height = 3, useDingbats=FALSE)
}




