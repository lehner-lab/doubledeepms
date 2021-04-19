
#' doubledeepms__plot_additive_trait_folding
#'
#' Plot folding additive trait.
#'
#' @param input_dt data.table with model free energy estimates (required)
#' @param report_outpath output path for scatterplots (required)
#'
#' @return Nothing
#' @export
#' @import data.table
doubledeepms__plot_additive_trait_folding <- function(
  input_dt,
  report_outpath
  ){

  #Check foldiong data exists
  if(input_dt[dataset_binding==0,.N]==0){return()}

  plot_dt <- input_dt[dataset_binding==0][!duplicated(id)]
  temp_dt <- plot_dt[,.(additive_trait_folding, observed_fitness = predicted_fitness, mut_order=1, f_ddg_pred_conf)]
  temp_dt2 <- copy(temp_dt)
  temp_dt2[, mut_order := 2]
  temp_dt <- rbind(temp_dt, temp_dt2)
  d <- ggplot2::ggplot(plot_dt[mut_order>0],ggplot2::aes(additive_trait_folding, observed_fitness)) +
    ggplot2::stat_binhex(bins = 100, size = 0, color = "lightgrey") +
    ggplot2::scale_fill_gradientn(colours = c("white", "black"), trans = "log10") +
    ggplot2::xlab(expression(Delta*"G Folding (ddPCA)")) +
    ggplot2::ylab("Fitness (Abundance)") +
    ggplot2::geom_hline(yintercept = 0, linetype = 2) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2) +
    ggplot2::geom_vline(xintercept = plot_dt[mut_order==0,additive_trait_folding][1], color = "grey", linetype = 2) +
    ggplot2::geom_line(data = temp_dt, color = "red") +
    ggplot2::theme_bw() +
    ggplot2::facet_wrap(~mut_order, nrow = 1)
  ggplot2::ggsave(file.path(report_outpath, "dG_observed_folding_scatter.pdf"), d, width = 8, height = 4, useDingbats=FALSE)
  d <- ggplot2::ggplot(plot_dt[mut_order>0],ggplot2::aes(additive_trait_folding, observed_fitness, colour = f_ddg_pred_conf)) +
    ggplot2::geom_point(data = plot_dt[f_ddg_pred_conf==F & mut_order>0], size = 1, alpha = 1/4) +
    ggplot2::geom_point(data = plot_dt[f_ddg_pred_conf==T & mut_order>0], size = 1, alpha = 1/4) +
    ggplot2::xlab(expression(Delta*"G Folding (ddPCA)")) +
    ggplot2::ylab("Fitness (Abundance)") +
    ggplot2::geom_hline(yintercept = 0, linetype = 2) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2) +
    ggplot2::geom_vline(xintercept = plot_dt[mut_order==0,additive_trait_folding][1], color = "grey", linetype = 2) +
    ggplot2::geom_line(data = temp_dt, color = "red") +
    ggplot2::theme_bw() +
    ggplot2::facet_wrap(~mut_order, nrow = 1)
  ggplot2::ggsave(file.path(report_outpath, "dG_observed_folding_scatter_conf.pdf"), d, width = 8, height = 4, useDingbats=FALSE)
}
