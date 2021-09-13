#' doubledeepms__plot_fitness_replicates_cor
#'
#' Create results folder for analysis script plots and saved objects.
#'
#' @param input_dt processed data table with DiMSum fitness scores (required)
#' @param outpath output path for plots and saved objects (required)
#' @param val_inpath  path to experimental fitness validations (required)
#'
#' @return Nothing
#' @export
doubledeepms__plot_growth_validations <- function(
  input_dt, 
  outpath,
  val_inpath){
  
  #Check validation data exists
  if(is.na(val_inpath)){return()}
  
  singles_dt <- input_dt[Nham_aa <= 1 & protein == "GRB2-SH3",]
  singles_dt[, id := paste0(WT_AA, Pos, Mut)]
  singles_dt[WT == T, id := "GRB2"]

  val_dt <- data.table::fread(input = val_inpath, header = T)
  gr <- do.call("rbind", lapply(unique(val_dt$date_well), function(x){
    ds = val_dt[date_well == x,]
    time_od50 = max(ds[od600 <= max(ds$od600)/2, time])
    lmds = data.frame(t = ds[time <= time_od50+2 & time >= time_od50-2, time],
                      ods = log(ds[time <= time_od50+2 & time >= time_od50-2, od600]))
    lmfit = lm(formula = ods ~ t, data = lmds)
    data.frame(genotype = unique(ds$genotype),
               pca_type = unique(ds$pca_type),
               growth_rate_slope = lmfit$coefficients[2],
               date_well = x)
  }))

  gr_mean <- as.data.table(aggregate(growth_rate_slope ~ genotype+pca_type, data=gr[, c("growth_rate_slope", "genotype", "pca_type")], FUN = mean))
  
  plot_dt <- merge(gr_mean, singles_dt[, .(id, growthrate, pca_type)], by.x = c("genotype", "pca_type"), by.y=c("id", "pca_type"), all=F)
  plot_shapes <- c(21, 19)
  names(plot_shapes) <- c("Abundance", "Binding")
  
  p <- ggplot2::ggplot(plot_dt, ggplot2::aes(x=growth_rate_slope, y=growthrate)) +
    ggplot2::geom_point(size=2, ggplot2::aes(shape=pca_type)) +
    ggplot2::geom_smooth(method = "lm", formula = 'y~x', linetype = 2, color = "grey", se = F) +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "bottom", legend.direction = "vertical", plot.margin = ggplot2::unit(c(30, 30, 30, 30), "points")) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::scale_shape_manual("PCA assay",values = plot_shapes) +
    ggrepel::geom_text_repel(ggplot2::aes(label = genotype), show.legend = F, 
                             max.overlaps = Inf, xlim = c(-10, 10), ylim = c(-10, 10)) +
    ggplot2::geom_text(data = plot_dt[,.(label = paste(" r = ", round(cor(growth_rate_slope, growthrate, use = "pairwise.complete", method = "pearson"), 2), sep=""))], 
                       ggplot2::aes(label=label, x=-Inf, y=Inf, hjust = 0, vjust = 1)) +
    ggplot2::labs(x="growth rate (individual measurements)",
                  y="growth rate (deep sequencing)")
    
  ggplot2::ggsave(file.path(outpath, "growthrate_validations_scatter.pdf"), width = 3.2, height = 4)

  #Correlation P-value
  temp_cor <- plot_dt[!is.na(growth_rate_slope) & !is.na(growthrate),cor.test(growth_rate_slope, growthrate)]
  print(paste0("Correlation with growth validations (pooled): p-value=", format(temp_cor[["p.value"]], scientific = T, digits = 2)))

}




