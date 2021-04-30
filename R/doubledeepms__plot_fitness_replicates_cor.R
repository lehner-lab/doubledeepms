
#' doubledeepms__plot_fitness_replicates_cor
#'
#' Create results folder for analysis script plots and saved objects.
#'
#' @param fitness_list list of folder paths to fitness data (required)
#' @param outpath output path for plots and saved objects (required)
#' @param execute whether or not the system command will be executed (required)
#'
#' @return Nothing
#' @export
doubledeepms__plot_fitness_replicates_cor <- function(
  fitness_list, 
  outpath,
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
  for(protein in names(fitness_list)){
    fit_list[[protein]] <- list()
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
    corm <- matrix(NA, nrow=length(fitness_reps_names), ncol = length(fitnes_reps_names))
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
  p1 <- ggplot2::ggplot(cor_df, aes(x=x, y=y, fill=r)) + geom_tile() +
    ggplot2::facet_grid(pca_type ~ protein) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.title = element_blank()) +
    ggplot2::scale_fill_gradient(low = "#F7941E", high = "#9061A7", breaks = seq(0.85, 1, 0.05), labels = format(seq(0.85, 1, 0.05)))
  ggplot2::ggsave(file.path(outpath,"fitness_replicates_cor_heatmap.pdf"), p1, width = 5, height = 3.8)
}








