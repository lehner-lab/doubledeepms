
#' doubledeepms__plot_model_performance
#'
#' Plot model performance.
#'
#' @param input_dt data.table with model free energy estimates (required)
#' @param report_outpath output path for scatterplots (required)
#' @param highlight_colour colour for highlights (default:red)
#'
#' @return Nothing
#' @export
#' @import data.table
doubledeepms__plot_model_performance <- function(
  input_dt,
  report_outpath,
  highlight_colour = "red"
  ){

  #Plot observed versus predicted fitness
  plot_height <- 5
  if(input_dt[dataset_binding==0,.N]==0 | input_dt[dataset_binding==0,.N]==0){plot_height <- 4}
  plot_dt <- input_dt[mut_order!=0 & !is.na(dataset_binding),.(id, dataset_binding, mut_order, observed_fitness, predicted_fitness)]
  plot_dt <- plot_dt[!duplicated(plot_dt[,.(id, dataset_binding)])]
  plot_dt[, dataset := "Folding"]
  plot_dt[dataset_binding==1, dataset := "Binding"]
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(observed_fitness, predicted_fitness)) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2) +
    ggplot2::stat_binhex(bins = 50, size = 0.2, color = "grey") +
    ggplot2::scale_fill_gradientn(colours = c("white", "black"), trans = "log10") +
    # ggplot2::geom_point(alpha = 1/10) +
    ggplot2::xlab("Observed fitness") +
    ggplot2::ylab("Predicted fitness") +
    ggplot2::geom_text(data = plot_dt[,.(label = paste("R-squared = ", round(cor(observed_fitness, predicted_fitness, use = "pairwise.complete")^2, 2), sep="")),.(mut_order, dataset)], ggplot2::aes(label=label, x=-Inf, y=Inf, hjust = 0, vjust = 1)) +
    ggplot2::facet_grid(dataset~mut_order, scales = "free") + 
    ggplot2::theme_classic()
  if(length(input_dt[!is.na(dataset_binding),unique(dataset_binding)])==length(input_dt[training_set==1,unique(dataset_binding)])){
    #Training data includes all available phenotypes
    d <- d + ggplot2::geom_abline(color = highlight_colour, linetype = 2, size = 1)
  }else{
    #Training data doesn't include all available phenotypes (predictions will be in different units)
    d <- d + ggplot2::geom_smooth(method = "lm", formula = "y~x", color = highlight_colour, linetype = 2, size = 1, se = F)
  }
  ggplot2::ggsave(file.path(report_outpath, "fitness_observed_predicted_scatter_mutorderfacet.pdf"), d, width = 7, height = plot_height, useDingbats=FALSE)

  #Plot observed versus predicted fitness - not facetted on mutation order
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(observed_fitness, predicted_fitness)) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2) +
    ggplot2::stat_binhex(bins = 50, size = 0.2, color = "grey") +
    ggplot2::scale_fill_gradientn(colours = c("white", "black"), trans = "log10") +
    # ggplot2::geom_point(alpha = 1/10) +
    ggplot2::xlab("Observed fitness") +
    ggplot2::ylab("Predicted fitness") +
    ggplot2::geom_text(data = plot_dt[,.(label = paste("R-squared = ", round(cor(observed_fitness, predicted_fitness, use = "pairwise.complete")^2, 2), sep="")),.(dataset)], ggplot2::aes(label=label, x=-Inf, y=Inf, hjust = 0, vjust = 1)) +
    ggplot2::facet_wrap(dataset~., scales = "free", nrow = 2) + 
    ggplot2::theme_classic()
  if(length(input_dt[!is.na(dataset_binding),unique(dataset_binding)])==length(input_dt[training_set==1,unique(dataset_binding)])){
    #Training data includes all available phenotypes
    d <- d + ggplot2::geom_abline(color = highlight_colour, linetype = 2, size = 1)
  }else{
    #Training data doesn't include all available phenotypes (predictions will be in different units)
    d <- d + ggplot2::geom_smooth(method = "lm", formula = "y~x", color = highlight_colour, linetype = 2, size = 1, se = F)
  }
  ggplot2::ggsave(file.path(report_outpath, "fitness_observed_predicted_scatter.pdf"), d, width = 4, height = plot_height, useDingbats=FALSE)

  #Plot observed versus predicted fitness - not facetted on mutation order - only validation data
  #Check validation data exists
  if(input_dt[training_set==0,.N]!=0){
    plot_dt <- input_dt[training_set==0 & mut_order!=0 & !is.na(dataset_binding),.(id, dataset_binding, mut_order, observed_fitness, predicted_fitness)]
    plot_dt <- plot_dt[!duplicated(plot_dt[,.(id, dataset_binding)])]
    plot_dt[, dataset := "Folding"]
    plot_dt[dataset_binding==1, dataset := "Binding"]
    d <- ggplot2::ggplot(plot_dt,ggplot2::aes(observed_fitness, predicted_fitness)) +
      ggplot2::geom_hline(yintercept = 0, linetype = 2) +
      ggplot2::geom_vline(xintercept = 0, linetype = 2) +
      ggplot2::stat_binhex(bins = 50, size = 0.2, color = "grey") +
      ggplot2::scale_fill_gradientn(colours = c("white", "black"), trans = "log10") +
      # ggplot2::geom_point(alpha = 1/10) +
      ggplot2::xlab("Observed fitness") +
      ggplot2::ylab("Predicted fitness") +
      ggplot2::geom_text(data = plot_dt[,.(label = paste("R-squared = ", round(cor(observed_fitness, predicted_fitness, use = "pairwise.complete")^2, 2), sep="")),.(dataset)], ggplot2::aes(label=label, x=-Inf, y=Inf, hjust = 0, vjust = 1)) +
      ggplot2::facet_wrap(dataset~., scales = "free", nrow = 2) + 
      ggplot2::theme_classic()
    if(length(input_dt[!is.na(dataset_binding),unique(dataset_binding)])==length(input_dt[training_set==1,unique(dataset_binding)])){
      #Training data includes all available phenotypes
      d <- d + ggplot2::geom_abline(color = highlight_colour, linetype = 2, size = 1)
    }else{
      #Training data doesn't include all available phenotypes (predictions will be in different units)
      d <- d + ggplot2::geom_smooth(method = "lm", formula = "y~x", color = highlight_colour, linetype = 2, size = 1, se = F)
    }
    ggplot2::ggsave(file.path(report_outpath, "fitness_observed_predicted_scatter_val.pdf"), d, width = 4, height = plot_height, useDingbats=FALSE)
  }

}
