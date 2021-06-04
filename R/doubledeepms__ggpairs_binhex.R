
#' doubledeepms__ggpairs_binhex
#'
#' GGpairs plot for all variables in input data table with correlation in upper triangle and 2d binned hexagons in lower triangle.
#'
#' @param input_dt input data table (required)
#' @param output_file plot output path (required)
#' @param width plot width (default: 12)
#' @param height plot height (default: 12)
#' @param bins number of hexagon bins (default: 50)
#' @param xlab the title of the x-axis
#' @param ylab the title of the y-axis
#' @param title the plot title
#' @param label_size the size of the label text
#' @param plot_colours the lower and upper gradient limit colours
#'
#' @return Nothing
#' @export
doubledeepms__ggpairs_binhex <- function(
  input_dt, 
  output_file, 
  input_dt_upper=NULL,
  width = 12, 
  height = 12, 
  bins = 50,
  xlab = "x",
  ylab = "y",
  title = "",
  label_size = 4,
  plot_colours = c("lightgrey", "black")
  ){
  #Check if something to plot
  if(dim(input_dt)[1]==0){
    warning("doubledeepms__ggpairs_binhex.R: No data to plot (empty data.table 'input_dt').", call. = FALSE, immediate. = TRUE, noBreaks. = TRUE)
    return(NULL)
  }
  upper_plot <- list(continuous = "cor")
  if(!is.null(input_dt_upper)){
    upper_plot <- "blank"
  }
  d <- GGally::ggpairs(input_dt,
    columns = 1:dim(input_dt)[2],
    upper=upper_plot,
    lower="blank", xlab = xlab, ylab = ylab, title = title)
  for (x in 1:dim(input_dt)[2]){
    for (y in 1:dim(input_dt)[2]){
      temp_plot <- NULL
      if (y>x) {
        axis_lims <- range(input_dt, na.rm = T)
        plot_xy <- quantile(input_dt, probs = c(0.01), na.rm = T)
        temp_cor <- input_dt[,.(cor = cor(.SD, use = "pairwise.complete.obs")[1,2], nx = sum(!is.na(.SD[[1]])), ny = sum(!is.na(.SD[[2]]))),,.SDcols = c(colnames(input_dt)[x],colnames(input_dt)[y])]
        temp_plot <- ggplot2::ggplot(input_dt, ggplot2::aes_string(x=colnames(input_dt)[x],y=colnames(input_dt)[y])) + 
          ggplot2::coord_cartesian(xlim = axis_lims, ylim = axis_lims) +
          ggplot2::stat_binhex(bins=bins, size = 0.2) +
          ggplot2::geom_abline(linetype = 2) + 
          # ggplot2::annotate("text", label = paste0("R = ", round(temp_cor[,"cor"], 2), " (nx=", temp_cor[,"nx"], ", ny=", temp_cor[,"ny"], ")") , x = plot_xy, y = plot_xy, size = label_size) +
          ggplot2::annotate("text", label = round(temp_cor[,"cor"], 2) , x = plot_xy, y = plot_xy, size = label_size) +
          ggplot2::scale_fill_gradientn(colours=plot_colours,name = "Frequency", na.value=plot_colours[length(plot_colours)])
        d <- GGally::putPlot(d, temp_plot, y, x)
      }
      if (y<x & !is.null(input_dt_upper)){
        axis_lims <- range(input_dt_upper, na.rm = T)
        plot_xy <- quantile(input_dt_upper, probs = c(0.01), na.rm = T)
        temp_cor <- input_dt_upper[,.(cor = cor(.SD, use = "pairwise.complete.obs")[1,2], nx = sum(!is.na(.SD[[1]])), ny = sum(!is.na(.SD[[2]]))),,.SDcols = c(colnames(input_dt_upper)[x],colnames(input_dt_upper)[y])]
        temp_plot <- ggplot2::ggplot(input_dt_upper, ggplot2::aes_string(x=colnames(input_dt_upper)[x],y=colnames(input_dt_upper)[y])) + 
          ggplot2::coord_cartesian(xlim = axis_lims, ylim = axis_lims) +
          ggplot2::stat_binhex(bins=bins, size = 0.2) +
          ggplot2::geom_abline(linetype = 2) + 
          # ggplot2::annotate("text", label = paste0("R = ", round(temp_cor[,"cor"], 2), " (nx=", temp_cor[,"nx"], ", ny=", temp_cor[,"ny"], ")") , x = plot_xy, y = plot_xy, size = label_size) +
          ggplot2::annotate("text", label = round(temp_cor[,"cor"], 2) , x = plot_xy, y = plot_xy, size = label_size) +
          ggplot2::scale_fill_gradientn(colours=plot_colours,name = "Frequency", na.value=plot_colours[length(plot_colours)])
        d <- GGally::putPlot(d, temp_plot, y, x)
      }
    }
  }
  d <- d + ggplot2::theme_classic()
  suppressWarnings(suppressMessages(ggplot2::ggsave(output_file, d, width = width, height = height, useDingbats=FALSE)))
}
