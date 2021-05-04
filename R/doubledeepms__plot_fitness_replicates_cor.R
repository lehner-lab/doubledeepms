
#' doubledeepms__plot_fitness_replicates_cor
#'
#' Create results folder for analysis script plots and saved objects.
#'
#' @param fitness_list list of folder paths to fitness data (required)
#' @param outpath output path for plots and saved objects (required)
#' @param scatter_examples_list list of two vectors containing as elements the protein, assay, replicate X and replicate Y, to plot the fitness comparison (required)
#' @param execute whether or not the system command will be executed (required)
#'
#' @return Nothing
#' @export
doubledeepms__plot_fitness_replicates_cor <- function(
  fitness_list, 
  outpath,
  scatter_examples_list,
  execute = TRUE, 
  message = NULL){
  
  #Return if analysis not executed
  if(!execute){
    return()
  }
  
  #Display status
  message(paste("\n\n*******", "running stage: doubledeepms__plot_fitness_replicates_cor", "*******\n\n"))
  
  #Create output directory
  doubledeepms__create_dir(doubledeepms_dir = outpath)
  
  fit_list <- list()
  yrep_list <- list()
  for(protein in names(fitness_list)){
    fit_list[[protein]] <- list()
    # modify y-axis rep name to follow correct orientation in heatmap
    if (scatter_examples_list[[protein]][4] == "rep1"){
      yrep_list[[protein]] = 3
    } else if (scatter_examples_list[[protein]][4] == "rep3") {
      yrep_list[[protein]] = 1
    } else{
      yrep_list[[protein]] = 2
    }
    for(pca_type in c("Abundance", "Binding")){
      rdata_file <- list.files(file.path(fitness_list[[protein]], pca_type))
      load(file.path(fitness_list[[protein]], pca_type, rdata_file))
      all_variants[, protein := protein]
      all_variants[, pca_type := pca_type]
      #Single mutant position
      wt_seq <- unlist(strsplit(all_variants[WT==T,aa_seq], ""))
      #Single mutant position class
      all_variants[Nham_aa==1, Pos := which(unlist(strsplit(aa_seq, "")) != wt_seq), aa_seq]
      fit_list[[protein]][[pca_type]] <- all_variants
    }
  }
  fitness_dt <- rbindlist(unlist(fit_list, recursive = FALSE), fill = T)
  
  cor_df_list <- list()
  fitness_reps_names <- colnames(fitness_dt)[grepl("^fitness[1-9]_", colnames(fitness_dt))]
  cor_df_list <- apply(unique(fitness_dt[, .(protein, pca_type)]), 1, function(prot_pca){
    corm1 <- stats::cor(fitness_dt[protein == prot_pca[1] & pca_type == prot_pca[2] & Nham_aa == 1, ..fitness_reps_names], use = "pairwise.complete.obs")
    corm2 <- stats::cor(fitness_dt[protein == prot_pca[1] & pca_type == prot_pca[2] & Nham_aa == 2, ..fitness_reps_names], use = "pairwise.complete.obs")
    corm <- matrix(NA, nrow=length(fitness_reps_names), ncol = length(fitness_reps_names))
    repnames <- gsub("^fitness", "rep", gsub("_uncorr", "", fitness_reps_names))
    corm[upper.tri(corm, diag = T)] <- corm1[upper.tri(corm1, diag = T)]
    corm[lower.tri(corm)] <- corm2[lower.tri(corm2)]
    cor_df <- setnames(data.frame(corm), repnames)
    cor_df$y <- repnames
    cor_df$protein <- prot_pca[1]
    cor_df$pca_type <- prot_pca[2]
    reshape2::melt(cor_df, id.vars = c("y", "protein", "pca_type"), variable.name = "x", value.name = "r")
  })
  cor_df <- data.table::rbindlist(cor_df_list)
  cor_df$y <- factor(cor_df$y, levels = rev(unique(cor_df$y)))
  col_rect <- list("singles" = colour_scheme[["shade 0"]][[3]],
                   "doubles" = colour_scheme[["shade 0"]][[1]])

  p1 <- ggplot2::ggplot(cor_df, ggplot2::aes(x=x, y=y, fill=r)) + geom_tile() +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.title = element_blank()) +
    ggplot2::scale_fill_gradient(low = "#F7941E", high = "#9061A7", breaks = seq(0.85, 1, 0.05), labels = format(seq(0.85, 1, 0.05))) +
    ggplot2::theme(legend.position = "bottom", legend.direction = "horizontal") +
    ggplot2::guides( fill = guide_colourbar(barwidth = 10, barheight = 0.5)) +
    ggplot2::geom_rect(inherit.aes = F, 
                       data = data.frame(protein = names(scatter_examples_list)[1],
                                         pca_type =  scatter_examples_list[[1]][1],
                                         xmin = as.numeric(gsub("rep", "", scatter_examples_list[[1]][3]))-0.5,
                                         xmax = as.numeric(gsub("rep", "", scatter_examples_list[[1]][3]))+0.5,
                                         ymin = yrep_list[[names(scatter_examples_list)[1]]]-0.5,
                                         ymax = yrep_list[[names(scatter_examples_list)[1]]]+0.5),
                       ggplot2::aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
                       color=col_rect[[scatter_examples_list[[1]][2]]], fill="white", alpha=0, size=1) +
    ggplot2::geom_rect(inherit.aes = F, 
                       data = data.frame(protein = names(scatter_examples_list)[2],
                                         pca_type =  scatter_examples_list[[2]][1],
                                         xmin = as.numeric(gsub("rep", "", scatter_examples_list[[2]][3]))-0.5,
                                         xmax = as.numeric(gsub("rep", "", scatter_examples_list[[2]][3]))+0.5,
                                         ymin = yrep_list[[names(scatter_examples_list)[2]]]-0.5,
                                         ymax = yrep_list[[names(scatter_examples_list)[2]]]+0.5),
                       ggplot2::aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
                       color=col_rect[[scatter_examples_list[[2]][2]]], fill="white", alpha=0, size=1) +
    ggplot2::facet_grid(pca_type ~ protein) 
  ggplot2::ggsave(file.path(outpath,"fitness_replicates_cor_heatmap.pdf"), p1, width = 4, height = 4.3)
  
  
  p2 <- do.call("grid.arrange", c(lapply(1:length(scatter_examples_list), function(x){
    reps_vector <- c(paste0(gsub("rep", "fitness", scatter_examples_list[[x]][3]), "_uncorr"), paste0(gsub("rep", "fitness", scatter_examples_list[[x]][4]), "_uncorr"))
    plot_dt <- setnames(fitness_dt[protein == names(scatter_examples_list)[x] & pca_type == scatter_examples_list[[x]][1], ..reps_vector], c("x","y"))
    ggplot2::ggplot(plot_dt, aes(x=x, y=y)) +
      ggplot2::stat_binhex(bins = 50, size = 0.2, color = "grey") +
      ggplot2::scale_fill_gradientn(colours = c("white", "black"), trans = "log10") +
      ggplot2::xlab(scatter_examples_list[[x]][3]) +
      ggplot2::ylab(scatter_examples_list[[x]][4]) +
      ggplot2::geom_abline(color = "black", linetype = 2, size = 1) +
      ggplot2::geom_text(data = plot_dt[,.(label = paste(" r = ", round(cor(x, y, use = "pairwise.complete"), 2), sep=""))], ggplot2::aes(label=label, x=-Inf, y=Inf, hjust = 0, vjust = 1)) +
      ggplot2::theme_classic() +
      ggplot2::labs(title = paste0(names(scatter_examples_list)[x], "\n", scatter_examples_list[[x]][1], ", ", scatter_examples_list[[x]][2]), color=col_rect[[x]]) +
      ggplot2::theme(plot.title = element_text(color = col_rect[[x]])) +
      ggplot2::scale_x_continuous(limits = c(-2, 0.5), breaks = c(-2,-1,0)) +
      ggplot2::scale_y_continuous(limits = c(-2, 0.5), breaks = c(-2,-1,0))
  }), ncol=1))
  ggplot2::ggsave(file.path(outpath,"fitness_replicates_scatters.pdf"), p2, width = 2.5, height = 4.3)
}






