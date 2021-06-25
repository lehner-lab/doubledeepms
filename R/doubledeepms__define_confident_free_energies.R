
#' doubledeepms__define_confident_free_energies
#'
#' Call confident ddGs.
#'
#' @param input_dt data.table with model free energy estimates (required)
#' @param report_outpath output path for scatterplots (required)
#' @param highlight_colour colour for highlights (default:red)
#' @param folding_energy_max_sd maximum folding energy standard deviation (default:1/(1.96*2))
#' @param binding_energy_max_sd maximum binding energy standard deviation (default:1/(1.96*2))
#'
#' @return data.table with columns indicating confident free energies
#' @export
#' @import data.table
doubledeepms__define_confident_free_energies <- function(
  input_dt,
  report_outpath,
  highlight_colour = "red",
  folding_energy_max_sd = 1/(1.96*2),
  binding_energy_max_sd = 1/(1.96*2)
  ){

  #ddG error distributions
  plot_dt <- reshape2::melt(copy(input_dt)[mut_order==1 & !duplicated(id),.(id, mut_order, f_ddg_pred_sd, b_ddg_pred_sd)], id = c("id", "mut_order"))
  plot_dt[, type := "Folding"]
  plot_dt[variable=="b_ddg_pred_sd", type := "Binding"]
  d <- ggplot2::ggplot(plot_dt[!is.na(value)],ggplot2::aes(log10(value))) +
    ggplot2::geom_density() +
    ggplot2::geom_vline(data = data.table(value = 1/(1.96*2), type = "Folding"), ggplot2::aes(xintercept = log10(value)), linetype = 2) +
    ggplot2::geom_vline(data = data.table(value = folding_energy_max_sd, type = "Folding"), ggplot2::aes(xintercept = log10(value)), linetype = 2, color = highlight_colour) +
    ggplot2::geom_vline(data = data.table(value = 1/(1.96*2), type = "Binding"), ggplot2::aes(xintercept = log10(value)), linetype = 2) +
    ggplot2::geom_vline(data = data.table(value = binding_energy_max_sd, type = "Binding"), ggplot2::aes(xintercept = log10(value)), linetype = 2, color = highlight_colour) +
    ggplot2::xlab(expression("log10 SD("*Delta*Delta*"G)")) +
    ggplot2::facet_grid(~type) + 
    ggplot2::theme_bw()
  ggplot2::ggsave(file.path(report_outpath, "ddG_sd_densities_singles.pdf"), d, width = 5, height = 3, useDingbats=FALSE)

  #Confident ddGs
  f_ids <- input_dt[f_ddg_pred_sd<folding_energy_max_sd,unique(id)]
  b_ids <- input_dt[b_ddg_pred_sd<binding_energy_max_sd,unique(id)]
  input_dt[, f_ddg_pred_conf := FALSE]
  input_dt[id %in% f_ids, f_ddg_pred_conf := TRUE]
  input_dt[, b_ddg_pred_conf := FALSE]
  input_dt[id %in% b_ids, b_ddg_pred_conf := TRUE]

  # #Upper and lower limits on ddG
  # b_ddg_pred_lower <- input_dt[b_ddg_pred_conf==T,quantile(b_ddg_pred, 0.025),]
  # b_ddg_pred_upper <- input_dt[b_ddg_pred_conf==T,quantile(b_ddg_pred, 0.975),]
  # f_ddg_pred_lower <- input_dt[f_ddg_pred_conf==T,quantile(f_ddg_pred, 0.025),]
  # f_ddg_pred_upper <- input_dt[f_ddg_pred_conf==T,quantile(f_ddg_pred, 0.975),]
  # input_dt[f_ddg_pred_conf==T & (f_ddg_pred > f_ddg_pred_upper | f_ddg_pred < f_ddg_pred_lower),.N,mut_order]
  # input_dt[b_ddg_pred_conf==T & (b_ddg_pred > b_ddg_pred_upper | b_ddg_pred < b_ddg_pred_lower),.N,mut_order]
  # input_dt[(f_ddg_pred > f_ddg_pred_upper | f_ddg_pred < f_ddg_pred_lower),f_ddg_pred_conf := F]
  # input_dt[(b_ddg_pred > b_ddg_pred_upper | b_ddg_pred < b_ddg_pred_lower),b_ddg_pred_conf := F]

  return(input_dt)
}
