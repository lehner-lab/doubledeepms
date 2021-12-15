
#' doubledeepms__persite_energy_vs_distance_plot
#'
#' Create results folder for analysis script plots and saved objects.
#'
#' @param input_dt input data table (required)
#' @param outpath plot output path (required)
#' @param colour_scheme colour scheme file (required)
#' @param literature_sites literature allosteric sites (default:empty vector)
#' @param trait_name trait name (default:folding)
#' @param absolute_value absolute value of the trait (default:T)
#' @param threshold_residues residues used to determine regulatory threshold: either "proximal" or "distal" (default:"proximal")
#'
#' @return Allosteric positions
#' @export
doubledeepms__persite_energy_vs_distance_plot <- function(
  input_dt, 
  outpath,
  colour_scheme,
  literature_sites=c(), 
  trait_name="folding",
  absolute_value=T,
  threshold_residues="proximal"
  ){

  #Trait
  trait_prefix <- ""
  metric_name_plot <- ""
  if(trait_name=="folding"){
    trait_prefix <- "f_ddg"
    if(absolute_value){
      metric_name_plot <- expression("Weighted mean |Folding "*Delta*Delta*"G|")
    }else{
      metric_name_plot <- expression("Weighted mean Folding "*Delta*Delta*"G")
    }
  }
  if(trait_name=="binding"){
    trait_prefix <- "b_ddg"
    if(absolute_value){
      metric_name_plot <- expression("Weighted mean |Binding "*Delta*Delta*"G|")
    }else{
      metric_name_plot <- expression("Weighted mean Binding "*Delta*Delta*"G")
    }
  }

  #Metrics
  if(absolute_value){
    metric_name <- paste0(trait_prefix, "_wposmeanabs")
  }else{
    metric_name <- paste0(trait_prefix, "_wposmean")
  }
  metric_name_se <- paste0(trait_prefix, "_wposse")

  #Subset
  plot_dt <- input_dt[order(Pos_ref)][!duplicated(Pos_ref)]

  #Metric
  plot_dt[, plot_metric := .SD[[1]],,.SDcols = metric_name]
  plot_dt[, plot_metric_se := .SD[[1]],,.SDcols = metric_name_se]

  #Threshold for binding interface modulating residues
  if(threshold_residues=="proximal"){
    if(absolute_value){
      reg_threshold <- input_dt[Pos_class=="binding_interface",sum(abs(.SD[[1]])/.SD[[2]]^2, na.rm = T)/sum(1/.SD[[2]]^2, na.rm = T),.SDcols = paste0(trait_prefix, c("_pred", "_pred_sd"))]
    }else{
      reg_threshold <- input_dt[Pos_class=="binding_interface",sum(.SD[[1]]/.SD[[2]]^2, na.rm = T)/sum(1/.SD[[2]]^2, na.rm = T),.SDcols = paste0(trait_prefix, c("_pred", "_pred_sd"))]
    }
  }
  if(threshold_residues=="distal"){
    reg_threshold <- input_dt[scHAmin_ligand>=5][!duplicated(Pos_ref),mean(.SD[[1]]) + sd(.SD[[1]]),,.SDcols = metric_name]
  }

  #Residue class
  allostery_pos <- plot_dt[plot_metric>=reg_threshold,Pos_ref]
  plot_dt[, class := "Remainder"]
  if(trait_name=="binding"){
    plot_dt[Pos_ref %in% allostery_pos, class := "Allosteric"]
  }
  plot_dt[Pos_class=="binding_interface", class := "Binding\ninterface"]
  plot_dt[, class := factor(class, levels = c("Allosteric", "Binding\ninterface", "Remainder"))]
  plot_cols <- c(colour_scheme[["shade 0"]][[2]], colour_scheme[["shade 0"]][[1]], "darkgrey")
  names(plot_cols) <- c("Allosteric", "Binding\ninterface", "Remainder")

  plot_dt[, switch := F]
  plot_dt[Pos_ref %in% literature_sites, switch := T]

  #Plot
  d <- ggplot2::ggplot(plot_dt[!is.na(plot_metric)],ggplot2::aes(scHAmin_ligand, plot_metric, colour = class)) +
    ggplot2::geom_vline(xintercept = 5, linetype = 2)
  if(trait_name=="binding"){
    d <- d + 
      ggplot2::geom_hline(yintercept = reg_threshold, linetype = 2) +
      ggplot2::geom_point(alpha = 3/4, size = 3, ggplot2::aes(shape = switch))# +
  }else{
    d <- d + ggplot2::geom_point(alpha = 3/4, size = 3) 
  }
  d <- d +
    ggplot2::geom_linerange(ggplot2::aes(ymin = plot_metric-plot_metric_se*1.96, ymax = plot_metric+plot_metric_se*1.96)) +
    ggrepel::geom_text_repel(ggplot2::aes(label = Pos_ref), show.legend = F, 
                             max.overlaps = 10) +
    ggplot2::xlab(expression("Distance to ligand ("*ring(A)*")")) +
    ggplot2::ylab(metric_name_plot) +
    ggplot2::geom_text(ggplot2::aes(label = Pos_ref), colour = "black", size = 1) +
    # ggplot2::annotate("text", label=paste("Pearson's r = ", plot_dt[,round(cor(scHAmin_ligand, plot_metric, use = "pairwise.complete"), 2)], sep=""), x=-Inf, y=Inf, hjust = 0, vjust = 1) +
    ggplot2::theme_classic() +
    ggplot2::scale_colour_manual(values=plot_cols) +
    ggplot2::labs(color = "Residue\nposition")   
  suppressWarnings(ggplot2::ggsave(outpath, d, width = 4, height = 3, useDingbats=FALSE))

  #Return allosteric positions
  return(allostery_pos)

}




