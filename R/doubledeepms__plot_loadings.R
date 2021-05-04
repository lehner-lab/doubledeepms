
#' doubledeepms__plot_loadings
#'
#' PCA biplot scatterplot matrix.
#'
#' @param pca_obj result of call to princomp (required)
#' @param output_file path to output plot (required)
#' @param plot_categories a character vector of category names (required)
#' @param plot_colours a character vector of category colours (required)
#' @param mono_cols a character vector of monochrome colours (to plot in background)
#' @param comps a numeric vector of principle component indices (default:1:5)
#' @param plot_symbol plot symbol (default:19)
#' @param points_cex points cex (default:1)
#' @param width plot width (default:10)
#' @param height plot height (default:10)
#'
#' @return Nothing
#' @export
doubledeepms__plot_loadings <- function(
  pca_obj,
  output_file,
  plot_categories, 
  plot_colours, 
  mono_cols=c("white", "grey", "lightgrey", "darkgrey"),
  comps=1:5, 
  plot_symbol=19,
  points_cex=1,
  width=10,
  height=10
  ){
  #Variable loadings
  plot_df <- as.data.frame(pca_obj[["rotation"]])
  plot_df[,"category"] <- factor(plot_categories, levels = names(plot_colours))
  #Plot
  if(length(comps)>2){
    d <- GGally::ggpairs(plot_df,
      columns = comps,
      upper="blank",
      lower="blank", axisLabels = "internal")
    for (x in comps){
      for (y in comps){
        temp_plot <- NULL
        if (y<x) {
          temp_plot <- ggplot2::ggplot(plot_df, ggplot2::aes_string(x=colnames(plot_df)[x],y=colnames(plot_df)[y], color = "category")) + 
            ggplot2::geom_point(data = plot_df[plot_colours[plot_categories] %in% mono_cols,], shape = 1) + 
            ggplot2::geom_point(data = plot_df[!plot_colours[plot_categories] %in% mono_cols,], shape = plot_symbol) + 
            ggplot2::scale_colour_manual(values = plot_colours) +
            ggplot2::theme_bw() +
            ggplot2::geom_vline(xintercept = 0, linetype=2) +
            ggplot2::geom_hline(yintercept = 0, linetype=2)
          d <- GGally::putPlot(d, temp_plot, y, x)
        }
      }
    }
    #Save
    suppressWarnings(suppressMessages(ggplot2::ggsave(file=output_file, d, width=width, height=height, useDingbats=FALSE)))
  }else{
    d <- ggplot2::ggplot(plot_df, ggplot2::aes_string(colnames(plot_df)[comps[1]], colnames(plot_df)[comps[2]], color = "category")) +
      ggplot2::geom_point(data = plot_df[plot_colours[plot_categories] %in% mono_cols,], shape = 1) + 
      ggplot2::geom_point(data = plot_df[!plot_colours[plot_categories] %in% mono_cols,], shape = plot_symbol) + 
      ggplot2::scale_colour_manual(values = plot_colours) +
      ggplot2::theme_bw() +
      ggplot2::geom_vline(xintercept = 0, linetype=2) +
      ggplot2::geom_hline(yintercept = 0, linetype=2)
    #Save
    ggplot2::ggsave(file=output_file, d, width=width, height=height, useDingbats=FALSE)
  }
}
