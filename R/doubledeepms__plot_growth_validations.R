#' doubledeepms__plot_fitness_replicates_cor
#'
#' Create results folder for analysis script plots and saved objects.
#'
#' @param input_dt processed data table with DiMSum fitness scores (required)
#' @param outpath output path for plots and saved objects (required)
#' @param val_inpath  path to experimental fitness validations (required)
#' @param colour_scheme colour scheme file (required)
#'
#' @return Nothing
#' @export
doubledeepms__plot_growth_validations <- function(
  input_dt, 
  outpath,
  val_inpath,
  colour_scheme){
  
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
  # gr_mean[pca_type == "Abundance", growth_rate_slope_rel2wt := growth_rate_slope-gr_mean[pca_type == "Abundance" & genotype == "GRB2", growth_rate_slope]]
  # gr_mean[pca_type == "Binding", growth_rate_slope_rel2wt := growth_rate_slope-gr_mean[pca_type == "Abundance" & genotype == "GRB2", growth_rate_slope]]
  
  plot_dt <- merge(gr_mean, singles_dt[, .(id, fitness, pca_type)], by.x = c("genotype", "pca_type"), by.y=c("id", "pca_type"), all=F)
  
  p <- ggplot2::ggplot(plot_dt, ggplot2::aes(x=growth_rate_slope, y=fitness)) +
    ggplot2::geom_point(size=2, ggplot2::aes(color=pca_type)) +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "bottom", legend.direction = "vertical") +
    ggplot2::scale_color_manual("PCA assay",values = c("#9061A7", "#F7941E")) +
    ggrepel::geom_text_repel(ggplot2::aes(label = genotype, color=pca_type), show.legend = F, max.overlaps = 20) +
    ggplot2::geom_text(data = plot_dt[,.(label = paste(" r = ", round(cor(growth_rate_slope, fitness, use = "pairwise.complete"), 2), sep=""))], 
                       ggplot2::aes(label=label, x=-Inf, y=Inf, hjust = 0, vjust = 1)) +
    ggplot2::labs(x="growth rate (individual measurements)",
                  y="fitness (deep sequencing)")
    
    
  ggplot2::ggsave(file.path(outpath, "fitness_growth_validations_scatter.pdf"), width = 3.2, height = 4)
}




