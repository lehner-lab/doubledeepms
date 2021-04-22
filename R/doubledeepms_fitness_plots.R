
#' doubledeepms_fitness_plots
#'
#' Evaluate thermo model results.
#'
#' @param fitness_list list of folder paths to fitness data (required)
#' @param structure_metrics_list list of paths to structure metrics (required)
#' @param outpath output path for plots and saved objects (required)
#' @param colour_scheme colour scheme file (required)
#' @param execute whether or not to execute the analysis (default: TRUE)
#'
#' @return Nothing
#' @export
#' @import data.table
doubledeepms_fitness_plots <- function(
  fitness_list,
  structure_metrics_list,
  outpath,
  colour_scheme,
  execute = TRUE
  ){

  #Return if analysis not executed
  if(!execute){
    return()
  }

  #Display status
  message(paste("\n\n*******", "running stage: doubledeepms_fitness_plots", "*******\n\n"))

  #Create output directory
  doubledeepms__create_dir(doubledeepms_dir = outpath)

  #Load fitness and metrics data
  fit_list <- list()
  for(protein in names(fitness_list)){
    fit_list[[protein]] <- list()
    #Structure metrics
    metrics_dt <- unique(fread(structure_metrics_list[[protein]])[,.(Pos, Pos_class)])
    for(pca_type in c("Abundance", "Binding")){
      rdata_file <- list.files(file.path(fitness_list[[protein]], pca_type))
      load(file.path(fitness_list[[protein]], pca_type, rdata_file))
      all_variants[, protein := protein]
      all_variants[, pca_type := pca_type]
      #Single mutant position
      wt_seq <- unlist(strsplit(all_variants[WT==T,aa_seq], ""))
      #Single mutant position class
      all_variants[Nham_aa==1, Pos := which(unlist(strsplit(aa_seq, "")) != wt_seq), aa_seq]
      all_variants <- merge(all_variants, metrics_dt, by = "Pos", all = T)
      fit_list[[protein]][[pca_type]] <- all_variants
    }
  }
  fitness_dt <- rbindlist(unlist(fit_list, recursive = FALSE), fill = T)

  #Fitness distributions by protein and PCA type
  plot_dt <- copy(fitness_dt)
  #Detrimental STOPs
  plot_dt[STOP==T, STOP_pos := sapply(lapply(strsplit(aa_seq, "\\*"), "[", 1), nchar)/nchar(gsub("*$", "", aa_seq))]
  plot_dt[STOP==T, STOP_detrimental := STOP_pos>1/4 & STOP_pos<3/4]
  #Plot
  d <- ggplot2::ggplot(plot_dt[is.na(WT) & STOP==F & STOP_readthrough==F & Nham_aa==Nmut_codons,],ggplot2::aes(fitness, fill = as.factor(Nham_aa))) +
    ggplot2::geom_density(alpha = 0.5) +
    ggplot2::geom_vline(data = plot_dt[WT==T,], ggplot2::aes(xintercept = fitness)) +
    ggplot2::geom_vline(data = plot_dt[STOP_detrimental==T,.(fitness = median(fitness)),.(pca_type, protein)], ggplot2::aes(xintercept = fitness), linetype = 2) +
    ggplot2::facet_wrap(pca_type~protein, scales = "free") +
    ggplot2::geom_text(data = plot_dt[Nham_aa==1 & STOP==F,.(label = paste0("#Singles = ", .N), mut_group = NA),.(pca_type, protein, Nham_aa)], ggplot2::aes(label=label, x=-Inf, y=Inf, hjust = 0, vjust = 1)) +
    ggplot2::geom_text(data = plot_dt[Nham_aa==2 & STOP==F,.(label = paste0("#Doubles = ", .N), mut_group = NA),.(pca_type, protein, Nham_aa)], ggplot2::aes(label=label, x=-Inf, y=-Inf, hjust = 0, vjust = 0)) +
    ggplot2::xlab("Fitness") +
    ggplot2::ylab("Density") +
    ggplot2::labs(fill = "Num. AA\nsubstitutions") +
    ggplot2::theme_classic()
  if(!is.null(colour_scheme)){
    d <- d + ggplot2::scale_fill_manual(values = unlist(colour_scheme[["shade 0"]][c(1, 3)]))
  }
  ggplot2::ggsave(file.path(outpath, "fitness_densities.pdf"), d, width = 7, height = 5, useDingbats=FALSE)

  #Fitness scatter by protein and position class
  plot_dt_abundance <- fitness_dt[Nham_aa==1 & STOP==F & STOP_readthrough==F & pca_type=="Abundance",.(aa_seq, fitness_abundance = fitness, Pos_class, protein)]
  plot_dt_binding <- fitness_dt[Nham_aa==1 & STOP==F & STOP_readthrough==F & pca_type=="Binding",.(aa_seq, fitness_binding = fitness, Pos_class, protein)]
  plot_dt <- merge(plot_dt_abundance, plot_dt_binding, by = c("aa_seq", "Pos_class", "protein"))
  plot_dt[, Pos_class_plot := stringr::str_to_title(Pos_class)]
  plot_dt[Pos_class=="binding_interface", Pos_class_plot := "Binding\ninterface"]
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(fitness_abundance, fitness_binding, colour = Pos_class_plot)) +
    ggplot2::geom_point(alpha = 0.5) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::geom_vline(xintercept = 0) +
    ggplot2::geom_abline(linetype = 2) +
    ggplot2::xlab("Fitness (Abundance)") +
    ggplot2::ylab("Fitness (Binding)") +
    ggplot2::labs(color = "Residue\nposition") +
    ggplot2::facet_wrap(~protein, scales = "free", nrow = 1) +
    ggplot2::theme_classic()
    # ggplot2::labs(color = "Residue\nposition")
  if(!is.null(colour_scheme)){
    d <- d + ggplot2::scale_colour_manual(values = unlist(colour_scheme[["shade 0"]][c(1, 3, 4)]))
  }
  ggplot2::ggsave(file.path(outpath, "fitness_scatter_singles_overlay.pdf"), d, width = 7, height = 3, useDingbats=FALSE)


}