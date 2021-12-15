
#' doubledeepms__plot_mutation_distributions
#'
#' Plot mutation distributions.
#'
#' @param input_dt data.table with fitness estimates (required)
#' @param report_outpath output path for scatterplots (required)
#' @param colour_scheme colour scheme file (required)
#'
#' @return Nothing
#' @export
#' @import data.table
doubledeepms__plot_mutation_distributions <- function(
  input_dt,
  report_outpath,
  colour_scheme
  ){

  #Plot mutation distributions (number of doubles for each single)
  mut_tab_abundance <- table(unlist(input_dt[dataset_binding==0 & mut_order==2,.(mut1, mut2)]))
  mut_tab_binding <- table(unlist(input_dt[dataset_binding==1 & mut_order==2,.(mut1, mut2)]))
  plot_dt <- data.table(
    count = c(as.numeric(mut_tab_abundance), as.numeric(mut_tab_binding)),
    phenotype = c(rep("abundance", length(mut_tab_abundance)), rep("binding", length(mut_tab_binding))))
  median_dt <- plot_dt[, .(median_count = median(count)),phenotype]
  plot_cols <- unlist(colour_scheme[["shade 0"]][c(1, 3)])
  names(plot_cols) <- c("binding", "abundance")
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(count, fill = phenotype, color = phenotype)) +
    # ggplot2::geom_density() +
    ggplot2::geom_histogram(bins = 30, position = "dodge", color = NA) +
    ggplot2::geom_vline(data = median_dt, ggplot2::aes(xintercept = median_count, color = phenotype), linetype = 2) + 
    ggplot2::geom_text(data = median_dt[,.(label = paste("median = ", median_count, sep=""), count = 50, y = 0.5),.(phenotype)], ggplot2::aes(x = count, y = y, color = phenotype, label = label)) +
    ggplot2::scale_x_continuous(trans = "log10") +
    ggplot2::xlab("Number of backgrounds (doubles per single)") +
    ggplot2::ylab("Number of single mutants") +
    ggplot2::theme_bw()
  if(!is.null(colour_scheme)){
    d <- d + ggplot2::scale_fill_manual(values = plot_cols)
    d <- d + ggplot2::scale_color_manual(values = plot_cols)
  }
  ggplot2::ggsave(file.path(report_outpath, "double_count_density.pdf"), d, width = 4, height = 3, useDingbats=FALSE)
  #Save median counts
  write.table(median_dt, file = file.path(report_outpath, "double_count_density.txt"), quote = F, row.names = F)
}
