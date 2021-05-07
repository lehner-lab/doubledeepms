
#' doubledeepms__tile_heatmap_wrapper
#'
#' ggplot tile heatmap wrapper.
#'
#' @param input_matrix matrix of heatmap values (required)
#' @param output_file plot output file path
#' @param width plot width in "units"
#' @param height plot height in "units"
#' @param units plot size units ("in", "cm", or "mm")
#' @param colour_clip maximum absolute value of colour scale
#' @param cluster heirarchically cluster ("none", "row", "column", "none")
#' @param xlab x-axis label
#' @param ylab y-axis label
#' @param xtick_labels display labels for x ticks
#' @param ytick_labels display labels for y ticks
#' @param colour_type colour scale type ("continuous", "categorical")
#' @param colour_low colour scale lower limit colour -- passed to scale_colour_gradient2
#' @param colour_high colour scale upper limit colour -- passed to scale_colour_gradient2 
#' @param colour_mid colour scale zero colour -- passed to scale_colour_gradient2
#' @param colour_midpoint midpoint for colour scale -- passed to scale_colour_gradient2
#' @param colour_limits upper and lower value limits of colour scale -- passed to scale_colour_gradient2
#' @param mono use monotype font (True, False)
#' @param na_colour colour to use for NA values
#' @param xaxis_angle rotation angle for x tick labels -- passed to element_text
#' @param xaxis_hjust horizontal justification of x tick labels (in [0, 1]) -- passed to element_text
#' @param xaxis_vjust vertical justification of x tick labels (in [0, 1]) -- passed to element_text
#' @param xaxis_size text size of x tick labels (in pts) -- passed to element_text
#' @param omit_xtext omit x tick labels (True, False)
#' @param omit_xticks omit x ticks (True, False)
#' @param yaxis_angle rotation angle for y tick labels -- passed to element_text
#' @param yaxis_hjust horizontal justification of y tick labels (in [0, 1]) -- passed to element_text
#' @param yaxis_vjust vertical justification of y tick labels (in [0, 1]) -- passed to element_text
#' @param yaxis_size text size of y tick labels (in pts) -- passed to element_text
#' @param omit_ytext omit y tick labels (True, False)
#' @param omit_yticks omit y ticks (True, False)
#' @param plot_title main title for plot
#' @param input_matrix_text matrix of heatmap text
#' @param input_matrix_point matrix of heatmap points (plotted instead of tiles)
#' @param text_size size of heatmap text
#' @param text_colour colour of heatmap text
#' @param highlight_regions list of highlighted regions of form: list("red" = list("region1" = c(_min_, _max_), "region2" = c(_min_, _max_), ...), "blue" = list("region3" = c(_min_, _max_), ...), ...)
#' @param x_breaks x-axis breaks (for displaing xtick_labels)
#' @param y_breaks y-axis breaks (for displaing ytick_labels)
#' @param plot whether to plot the heatmap (True, False)
#'
#' @return Nothing
#' @export
#' @import data.table
doubledeepms__tile_heatmap_wrapper <- function(
  input_matrix, 
  output_file, 
  width=10, 
  height=4, 
  units="in",
  colour_clip=4, 
  cluster='both', 
  xlab='x', 
  ylab='y', 
  xtick_labels=NULL, 
  ytick_labels=NULL,
  colour_type='continuous', 
  colour_low='blue', 
  colour_high='red', 
  colour_mid='white', 
  colour_midpoint=0, 
  colour_limits=NULL, 
  mono=F, 
  na_colour="grey50",
  xaxis_angle=330, 
  xaxis_hjust=0, 
  xaxis_vjust=NULL, 
  xaxis_size=5, 
  omit_xtext=F, 
  omit_xticks=F,
  yaxis_angle=NULL, 
  yaxis_hjust=NULL, 
  yaxis_vjust=NULL, 
  yaxis_size=NULL, 
  omit_ytext=F, 
  omit_yticks=F, 
  plot_title='', 
  input_matrix_text=NULL, 
  input_matrix_point=NULL, 
  text_size=0.25, 
  text_colour="black", 
  highlight_regions=NULL, 
  x_breaks=ggplot2::waiver(), 
  y_breaks=ggplot2::waiver(), 
  plot = T){

  order_row <- rev(1:dim(input_matrix)[1])
  order_col <- 1:dim(input_matrix)[2]
  if(cluster %in% c('both', 'row')){
    d <- dist(input_matrix, method = "euclidean") # distance matrix
    order_row <- hclust(d, method="ward.D")$order
  }
  if(cluster %in% c('both', 'column')){
    d <- dist(t(input_matrix), method = "euclidean") # distance matrix
    order_col <- hclust(d, method="ward.D")$order    
  }
  plot_df <- reshape2::melt(input_matrix[order_row,order_col])
  colnames(plot_df) <- c('y', 'x', 'value')
  plot_df$label <- ""
  if(!is.null(input_matrix_text)){
    plot_df_text <- reshape2::melt(input_matrix_text[order_row,order_col])
    colnames(plot_df_text) <- c('y', 'x', 'label')
    plot_df$label <- plot_df_text$label
  }
  if(colour_type=='continuous' & colour_clip){
    plot_df$value[plot_df$value>colour_clip] <- colour_clip
    plot_df$value[plot_df$value<(-colour_clip)] <- (-colour_clip)    
  }
  if(!is.null(input_matrix_point)){
    plot_df_point <- reshape2::melt(input_matrix_point[order_row,order_col])
    colnames(plot_df_point) <- c('y', 'x', 'point')
    plot_df$point <- plot_df_point$point
    #Plot points on top of 0 tiles
    plot_df$value_point <- plot_df$value
    plot_df[plot_df$point==T,]$value <- 0
  }
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x, y)) + 
    ggplot2::geom_tile(ggplot2::aes(fill = value)) + 
    ggplot2::geom_text(ggplot2::aes(label = label), size=text_size, colour=text_colour) +
    ggplot2::theme(
      axis.text.x=list(
        ggplot2::element_text(angle = xaxis_angle, hjust = xaxis_hjust, vjust = xaxis_vjust, size = xaxis_size, family=c('', 'mono')[as.numeric(mono)+1]), 
        ggplot2::element_blank())[[as.numeric(omit_xtext)+1]],
      axis.text.y=list(
        ggplot2::element_text(angle = yaxis_angle, hjust = yaxis_hjust, vjust = yaxis_vjust, size = yaxis_size, family=c('', 'mono')[as.numeric(mono)+1]), 
        ggplot2::element_blank())[[as.numeric(omit_ytext)+1]],
      axis.ticks.x=list(
        ggplot2::element_line(), 
        ggplot2::element_blank())[[as.numeric(omit_xticks)+1]],
      axis.ticks.y=list(
        ggplot2::element_line(), 
        ggplot2::element_blank())[[as.numeric(omit_yticks)+1]],
      panel.grid.major = ggplot2::element_blank(), 
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank()) + 
    # ggplot2::theme_classic() + 
    ggplot2::xlab(xlab) + 
    ggplot2::ylab(ylab) + 
    ggplot2::labs(title = plot_title)
  if(!is.null(highlight_regions)){
    aes_now <- function(...) {
      structure(list(...),  class = "uneval")
    }
    for(i in names(highlight_regions)){
      for(j in names(highlight_regions[[i]])){
        p <- p + 
          ggplot2::geom_rect(
            data = NULL, 
            mapping = aes_now(
              xmin=highlight_regions[[i]][[j]][1]-0.5, 
              xmax=highlight_regions[[i]][[j]][2]+0.5, 
              ymin=highlight_regions[[i]][[j]][1]-0.5, 
              ymax=highlight_regions[[i]][[j]][2]+0.5), 
            fill = NA, 
            colour = i)
      }
    }
  }
  if(!is.null(input_matrix_point)){
    p <- p + ggplot2::geom_point(data = plot_df[plot_df$point == T,], ggplot2::aes(colour = value_point), size = 2)
  }

  #xtick labels specified
  if(!is.null(xtick_labels)){
    if(is.numeric(plot_df$x)){
      #X is numeric
      p <- p + ggplot2::scale_x_continuous(breaks=x_breaks, labels=xtick_labels)
    }else{
      #X is discrete
      p <- p + ggplot2::scale_x_discrete(breaks=x_breaks, labels=xtick_labels)
    }
  }else{
    if(is.numeric(plot_df$x)){
      p <- p + ggplot2::scale_x_continuous(breaks=x_breaks)
    }
  }
  #ytick labels specified
  if(!is.null(ytick_labels)){
    if(is.numeric(plot_df$y)){
      #Y is numeric
      p <- p + ggplot2::scale_y_continuous(breaks=y_breaks, labels = ytick_labels)
    }else{
      #Y is discrete
      p <- p + ggplot2::scale_y_discrete(breaks=y_breaks, labels=ytick_labels)
    }
  }else{
    if(is.numeric(plot_df$y)){
      p <- p + ggplot2::scale_y_continuous(breaks=y_breaks)
    }
  }
  if(colour_type=='continuous'){
    p <- p + ggplot2::scale_fill_gradient2(low = colour_low, high = colour_high, mid = colour_mid, midpoint = colour_midpoint, limits=colour_limits, na.value=na_colour)    
    p <- p + ggplot2::scale_colour_gradient2(low = colour_low, high = colour_high, mid = colour_mid, midpoint = colour_midpoint, limits=colour_limits, na.value=na_colour)    
  }
  if(colour_type=='categorical'){
    p <- p + ggplot2::scale_fill_brewer(palette='Set1')    
    p <- p + ggplot2::scale_colour_brewer(palette='Set1')    
  }
  if(plot){
    ggplot2::ggsave(file=output_file, width=width, units=units, height=height, useDingbats=FALSE)
  }
  return(p)
}
