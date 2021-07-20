
#' doubledeepms_fitness_plots
#'
#' Plot fitness distributions and scatterplots.
#'
#' @param fitness_list list of folder paths to fitness data (required)
#' @param val_inpath  path to experimental fitness validations (required)
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
  val_inpath,
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
      if(length(rdata_file)!=0){
        load(file.path(fitness_list[[protein]], pca_type, rdata_file))
        all_variants[, protein := protein]
        all_variants[, pca_type := pca_type]
        #Single mutant position
        wt_seq <- unlist(strsplit(all_variants[WT==T,aa_seq], ""))
        #Single mutant position class
        all_variants[Nham_aa==1, Pos := which(unlist(strsplit(aa_seq, "")) != wt_seq), aa_seq]
        all_variants[Nham_aa==1, WT_AA := wt_seq[Pos], aa_seq]
        all_variants[Nham_aa==1, Mut := substr(aa_seq, Pos, Pos), aa_seq]
        all_variants <- merge(all_variants, metrics_dt, by = "Pos", all = T)
        fit_list[[protein]][[pca_type]] <- all_variants
      }
    }
  }
  fitness_dt <- rbindlist(unlist(fit_list, recursive = FALSE), fill = T)

  ### Plot fitness replicate correlations
  ###########################

  doubledeepms__plot_fitness_replicates_cor(
    input_dt = fitness_dt,
    outpath = outpath,
    colour_scheme = colour_scheme)
  
  ### Plot experimental growth validations vs. sequencing results
  ###########################
  doubledeepms__plot_growth_validations(
    input_dt = fitness_dt, 
    outpath = outpath,
    val_inpath = val_inpath)

  ### Plot fitness distributions by protein and PCA type
  ###########################

  #Fitness distributions by protein and PCA type
  plot_dt <- copy(fitness_dt[protein!="GB1"])
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
  ggplot2::ggsave(file.path(outpath, "fitness_densities.pdf"), d, width = 7, height = 3, useDingbats=FALSE)

  #Proportion of singles in bimodal peaks - PSD95-PDZ3
  strong_singles <- fitness_dt[pca_type=="Binding" & protein=="PSD95-PDZ3" & Nham_aa==1 & STOP==F & STOP_readthrough==F,sum(fitness<(-0.75))/.N*100]
  print(paste0("%Single AA substitutions strongly affecting binding in PSD95-PDZ3 (fitness < -0.75): ", round(strong_singles, 0)))
  weak_singles <- fitness_dt[pca_type=="Binding" & protein=="PSD95-PDZ3" & Nham_aa==1 & STOP==F & STOP_readthrough==F,sum(fitness>(-0.25))/.N*100]
  print(paste0("%Single AA substitutions strongly affecting binding in PSD95-PDZ3 (fitness > -0.25): ", round(weak_singles, 0)))
  #Proportion of singles in bimodal peaks - GRB2-SH3
  strong_singles <- fitness_dt[pca_type=="Binding" & protein=="GRB2-SH3" & Nham_aa==1 & STOP==F & STOP_readthrough==F,sum(fitness<(-0.75))/.N*100]
  print(paste0("%Single AA substitutions strongly affecting binding in GRB2-SH3 (fitness < -0.75): ", round(strong_singles, 0)))
  weak_singles <- fitness_dt[pca_type=="Binding" & protein=="GRB2-SH3" & Nham_aa==1 & STOP==F & STOP_readthrough==F,sum(fitness>(-0.25))/.N*100]
  print(paste0("%Single AA substitutions strongly affecting binding in GRB2-SH3 (fitness > -0.25): ", round(weak_singles, 0)))

  ### Plot fitness scatter by protein and position class
  ###########################

  #Fitness scatter by protein and position class
  plot_dt_abundance <- fitness_dt[protein!="GB1" & Nham_aa==1 & STOP==F & STOP_readthrough==F & pca_type=="Abundance",.(aa_seq, fitness_abundance = fitness, Pos_class, protein)]
  plot_dt_binding <- fitness_dt[protein!="GB1" & Nham_aa==1 & STOP==F & STOP_readthrough==F & pca_type=="Binding",.(aa_seq, fitness_binding = fitness, Pos_class, protein)]
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
  
  
  ###########################
  ### Where are strongest binding effects?
  ###########################

  fitness_dt_merge <- merge(
    fitness_dt[Nham_aa==1 & STOP==F & STOP_readthrough==F & pca_type=="Binding"], 
    fitness_dt[Nham_aa==1 & STOP==F & STOP_readthrough==F & pca_type=="Abundance", .(aa_seq, fitness_abundance = fitness, sigma_abundance = sigma)], by = "aa_seq", all.x = T)

  #Thresholded data
  plot_list <- list()
  for(i in c(1:40)/20){
    plot_list[[as.character(i)]] <- fitness_dt_merge[Nham_aa==1 & STOP==F & STOP_readthrough==F & (fitness-1.96*sigma)>i/10,]
    plot_list[[as.character(i)]][, fitness_class := i]
    plot_list[[as.character(-i)]] <- fitness_dt_merge[Nham_aa==1 & STOP==F & STOP_readthrough==F & (fitness-1.96*sigma)<(-i),]
    plot_list[[as.character(-i)]][, fitness_class := (-i)]
  }
  plot_dt <- rbindlist(plot_list)
  
  #Plot - all
  plot_dt_all <- plot_dt[,.(num_mutations = .N),.(fitness_class, Pos_class, protein)]
  d <- ggplot2::ggplot(plot_dt_all[order(Pos_class)],ggplot2::aes(fitness_class, num_mutations, color = Pos_class)) +
    ggplot2::geom_line(size = 1) +
    # ggplot2::geom_point(size = 0.5) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2) +
    ggplot2::facet_grid(protein~., scales = "free") +
    ggplot2::xlab("Binding fitness threshold") +
    ggplot2::ylab("#Mutations") +
    ggplot2::labs(color = "Residue\nposition") +
    ggplot2::theme_bw()
  if(!is.null(colour_scheme)){
    d <- d + ggplot2::scale_color_manual(values = unlist(colour_scheme[["shade 0"]][c(1, 3, 4)]))
  }
  ggplot2::ggsave(file.path(outpath, "binding_fitness_mutations_all.pdf"), d, width = 5, height = 7, useDingbats=FALSE)

  #Plot
  plot_dt_all <- plot_dt[protein!="GB1",.(num_mutations = .N),.(fitness_class, Pos_class, protein)]
  d <- ggplot2::ggplot(plot_dt_all[order(Pos_class)],ggplot2::aes(fitness_class, num_mutations, color = Pos_class)) +
    ggplot2::geom_line(size = 1) +
    # ggplot2::geom_point(size = 0.5) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2) +
    ggplot2::facet_wrap(~protein, nrow = 1, scales = "free") +
    ggplot2::xlab("Binding fitness threshold") +
    ggplot2::ylab("#Mutations") +
    ggplot2::labs(color = "Residue\nposition") +
    ggplot2::theme_bw()
  if(!is.null(colour_scheme)){
    d <- d + ggplot2::scale_color_manual(values = unlist(colour_scheme[["shade 0"]][c(1, 3, 4)]))
  }
  ggplot2::ggsave(file.path(outpath, "binding_fitness_mutations.pdf"), d, width = 7, height = 2, useDingbats=FALSE)

  #Plot - abundance unchanged
  plot_dt_all <- plot_dt[(abs(fitness_abundance)-1.96*sigma_abundance)<0,.(num_mutations = .N),.(fitness_class, Pos_class, protein)]
  d <- ggplot2::ggplot(plot_dt_all[order(Pos_class)],ggplot2::aes(fitness_class, num_mutations, color = Pos_class)) +
    ggplot2::geom_line(size = 1) +
    # ggplot2::geom_point(size = 0.5) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2) +
    ggplot2::facet_wrap(~protein, nrow = 1, scales = "free") +
    ggplot2::xlab("Binding fitness threshold") +
    ggplot2::ylab("#Mutations") +
    ggplot2::labs(color = "Residue\nposition") +
    ggplot2::theme_bw()
  if(!is.null(colour_scheme)){
    d <- d + ggplot2::scale_color_manual(values = unlist(colour_scheme[["shade 0"]][c(1, 3, 4)]))
  }
  ggplot2::ggsave(file.path(outpath, "binding_fitness_mutations_abundanceunchanged.pdf"), d, width = 7, height = 2, useDingbats=FALSE)
  
  
  ###########################
  ### Binding ROC plot
  ###########################
  
  #Per residue metrics
  fitness_dt_pos_list <- list()
  for(protein in names(fitness_list)){
    temp_dt <- fitness_dt_merge[protein == protein,]
    for(i in c("", "_abundance")){
      temp_dt[,paste0("fitness", i, "_posmeanabs") := mean(abs(.SD[[1]]), na.rm = T),Pos,.SDcols = paste0(c("fitness"), i)]
      temp_dt[,paste0("fitness", i, "_wposmeanabs") := sum(abs(.SD[[1]])/.SD[[2]]^2, na.rm = T)/sum(1/.SD[[2]]^2, na.rm = T),Pos,.SDcols = c(paste0("fitness", i), paste0("sigma", i))]
    }
    fitness_dt_pos_list[[protein]] <- temp_dt
  }
  fitness_dt_pos <- rbindlist(fitness_dt_pos_list)
  fitness_dt_pos[, Pos_ref := Pos]

  #Plot ROC curves
  for(i in unique(fitness_dt_pos[,protein])){
    doubledeepms__plot_binding_site_ROC(
      input_dt = copy(fitness_dt_pos)[protein==i & Pos!=0],
      outpath = file.path(outpath, paste0(i, "_fitness_binding_site_ROC.pdf")),
      colour_scheme = colour_scheme,
      metric_names <- c(
        "fitness_posmeanabs",
        "fitness_wposmeanabs",
        "fitness_abundance_posmeanabs",
        "fitness_abundance_wposmeanabs"),
      metric_names_plot <- c(
        "Mean |bindingPCA fitness|",
        "Weighted mean |bindingPCA fitness|",
        "Mean |abundancePCA fitness|",
        "Weighted mean |abundancePCA fitness|"))
  }
  
}