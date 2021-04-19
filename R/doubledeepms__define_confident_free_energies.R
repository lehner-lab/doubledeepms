
#' doubledeepms__define_confident_free_energies
#'
#' Call confident ddGs.
#'
#' @param input_dt data.table with model free energy estimates (required)
#' @param report_outpath output path for scatterplots (required)
#' @param folding_energy_max_sd maximum folding energy standard deviation (default:1/(1.96*2))
#' @param binding_energy_max_sd maximum binding energy standard deviation (default:1/(1.96*2))
#'
#' @return data.table with columns indicating confident free energies
#' @export
#' @import data.table
doubledeepms__define_confident_free_energies <- function(
  input_dt,
  report_outpath,
  folding_energy_max_sd = 1/(1.96*2),
  binding_energy_max_sd = 1/(1.96*2)
  ){

  #ddG error distributions
  plot_dt <- reshape2::melt(copy(input_dt)[mut_order==1 & !duplicated(id),.(id, mut_order, f_ddg_pred_sd, b_ddg_pred_sd)], id = c("id", "mut_order"))
  plot_dt[, type := "Folding"]
  plot_dt[variable=="b_ddg_pred_sd", type := "Binding"]
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(log10(value))) +
    ggplot2::geom_density() +
    ggplot2::geom_vline(data = data.table(value = 1/(1.96*2), type = "Folding"), ggplot2::aes(xintercept = log10(value)), linetype = 2) +
    ggplot2::geom_vline(data = data.table(value = folding_energy_max_sd, type = "Folding"), ggplot2::aes(xintercept = log10(value)), linetype = 2, color = "red") +
    ggplot2::geom_vline(data = data.table(value = 1/(1.96*2), type = "Binding"), ggplot2::aes(xintercept = log10(value)), linetype = 2) +
    ggplot2::geom_vline(data = data.table(value = binding_energy_max_sd, type = "Binding"), ggplot2::aes(xintercept = log10(value)), linetype = 2, color = "red") +
    ggplot2::xlab(expression("log10 SD("*Delta*Delta*"G)")) +
    ggplot2::facet_grid(~type) + 
    ggplot2::theme_bw()
  ggplot2::ggsave(file.path(report_outpath, "ddG_sd_densities_singles.pdf"), d, width = 5, height = 3, useDingbats=FALSE)

  # #Residual fitness distributions
  # input_dt[, fitness_resid := predicted_fitness - observed_fitness]
  # plot_dt <- copy(input_dt)
  # #Outlier residuals thresholds
  # f_resid_thresh_upper <- plot_dt[dataset_binding==0 & mut_order==1][!duplicated(id), quantile(fitness_resid, 0.75)+1.5*(quantile(fitness_resid, 0.75)-quantile(fitness_resid, 0.25))]
  # f_resid_thresh_lower <- plot_dt[dataset_binding==0 & mut_order==1][!duplicated(id), quantile(fitness_resid, 0.25)-1.5*(quantile(fitness_resid, 0.75)-quantile(fitness_resid, 0.25))]
  # b_resid_thresh_upper <- plot_dt[dataset_binding==1 & mut_order==1][!duplicated(id), quantile(fitness_resid, 0.75)+1.5*(quantile(fitness_resid, 0.75)-quantile(fitness_resid, 0.25))]
  # b_resid_thresh_lower <- plot_dt[dataset_binding==1 & mut_order==1][!duplicated(id), quantile(fitness_resid, 0.25)-1.5*(quantile(fitness_resid, 0.75)-quantile(fitness_resid, 0.25))]

  # plot_dt <- reshape2::melt(plot_dt[mut_order==1 & !is.na(dataset_binding),.(id, dataset_binding, fitness_resid)], id = c("id", "dataset_binding", "fitness_resid"))
  # plot_dt <- plot_dt[!duplicated(plot_dt[,.(id, dataset_binding)])]
  # plot_dt[, dataset := "Folding"]
  # plot_dt[dataset_binding==1, dataset := "Binding"]
  # d <- ggplot2::ggplot(plot_dt,ggplot2::aes(fitness_resid)) +
  #   ggplot2::geom_density() +
  #   ggplot2::geom_vline(data = data.table(fitness_resid = c(f_resid_thresh_lower, f_resid_thresh_upper), dataset = "Folding"), ggplot2::aes(xintercept = fitness_resid), linetype = 2, color = "red") +
  #   ggplot2::geom_vline(data = data.table(fitness_resid = c(b_resid_thresh_lower, b_resid_thresh_upper), dataset = "Binding"), ggplot2::aes(xintercept = fitness_resid), linetype = 2, color = "red") +
  #   ggplot2::xlab("Residual Fitness (Predicted - Observed)") +
  #   ggplot2::facet_grid(~dataset) + 
  #   ggplot2::theme_bw()
  # ggplot2::ggsave(file.path(report_outpath, "fitness_resid_densities_singles.pdf"), d, width = 5, height = 3, useDingbats=FALSE)

  #Confident ddGs
  # nf_ids <- input_dt[dataset_binding==0 & (fitness_resid>=f_resid_thresh_upper | fitness_resid<=f_resid_thresh_lower),unique(id)]
  # nb_ids <- input_dt[dataset_binding==1 & (fitness_resid>=b_resid_thresh_upper | fitness_resid<=b_resid_thresh_lower),unique(id)]
  # f_ids <- input_dt[f_ddg_pred_sd<folding_energy_max_sd & !id %in% nf_ids,unique(id)]
  # b_ids <- input_dt[b_ddg_pred_sd<binding_energy_max_sd & !id %in% nb_ids,unique(id)]
  f_ids <- input_dt[f_ddg_pred_sd<folding_energy_max_sd,unique(id)]
  b_ids <- input_dt[b_ddg_pred_sd<binding_energy_max_sd,unique(id)]
  input_dt[, f_ddg_pred_conf := FALSE]
  input_dt[id %in% f_ids, f_ddg_pred_conf := TRUE]
  input_dt[, b_ddg_pred_conf := FALSE]
  input_dt[id %in% b_ids, b_ddg_pred_conf := TRUE]
  return(input_dt)
}
