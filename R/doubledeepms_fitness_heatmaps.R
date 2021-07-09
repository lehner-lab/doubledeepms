
#' doubledeepms_fitness_heatmaps
#'
#' Plot fitness heatmaps.
#'
#' @param input_file path to MoCHI thermo model fit results (required)
#' @param input_file_fitness path to fitness data (required)
#' @param input_file_MSA path to MSA frequencies data (optional)
#' @param domain_name domain name (required)
#' @param outpath output path for plots and saved objects (required)
#' @param colour_scheme colour scheme file (required)
#' @param plot_width heatmap plot width in inches (default:10)
#' @param plot_height heatmap plot height in inches (default:4)
#' @param plot_traits traits to plot (default:Abundance, Binding)
#' @param execute whether or not to execute the analysis (default: TRUE)
#'
#' @return Nothing
#' @export
#' @import data.table
doubledeepms_fitness_heatmaps <- function(
  input_file,
  input_file_fitness,
  input_file_MSA = NULL,
  domain_name,
  outpath,
  colour_scheme,
  plot_width = 10,
  plot_height = 4,
  plot_traits = c("Abundance", "Binding"),
  execute = TRUE
  ){

  #Return if analysis not executed
  if(!execute){
    return()
  }

  #Display status
  message(paste("\n\n*******", paste("running stage: doubledeepms_fitness_heatmaps for", domain_name), "*******\n\n"))

  #Create output directory
  doubledeepms__create_dir(doubledeepms_dir = outpath)

  #Load fitness and metrics data
  fit_list <- list()
  #Load structure metrics
  metrics_dt <- unique(fread(input_file)[id!="-0-",.(Pos, Pos_class, scHAmin_ligand, Pos_ref)])
  for(pca_type in plot_traits){
    rdata_file <- list.files(file.path(input_file_fitness, pca_type))
    load(file.path(input_file_fitness, pca_type, rdata_file))
    all_variants[Nham_aa<=1 & STOP==F & STOP_readthrough==F, pca_type := pca_type]
    #Single mutant position
    wt_seq <- unlist(strsplit(all_variants[WT==T,aa_seq], ""))
    #Single mutant position class
    all_variants[Nham_aa==1, Pos := which(unlist(strsplit(aa_seq, "")) != wt_seq), aa_seq]
    all_variants[Nham_aa==1, WT_AA := wt_seq[Pos], aa_seq]
    all_variants[Nham_aa==1, Mut := substr(aa_seq, Pos, Pos), aa_seq]
    all_variants <- merge(all_variants, metrics_dt, by = "Pos", all = T)
    fit_list[[pca_type]] <- all_variants
  }
  fitness_dt <- rbindlist(fit_list, fill = T)

  ###########################
  ### Folding heatmap
  ###########################

  if("Abundance" %in% plot_traits){
    #Heatmap data
    heatmap_dt <- fitness_dt[pca_type=="Abundance" & Nham_aa==1,.(fitness = -fitness, fitness_conf = T, Pos_ref, scHAmin_ligand, WT_AA, Mut, Pos_class)]
    #Add dummy mutation to ensure all positions included
    dummy_dt <- fitness_dt[Nham_aa==1 & !is.na(Pos_ref)][!duplicated(Pos_ref) & Nham_aa==1,.(fitness = -fitness, fitness_conf = T, Pos_ref, scHAmin_ligand, WT_AA, Mut, Pos_class)]
    dummy_dt[, Mut := WT_AA]
    dummy_dt[, fitness := NA]
    heatmap_dt <- rbind(heatmap_dt, dummy_dt)

    doubledeepms__plot_heatmap(
      input_dt = heatmap_dt,
      variable_name = "fitness",
      output_file = file.path(outpath, "abundance_heatmap.pdf"),
      width = plot_width,
      height = plot_height,
      plot_title = paste0(domain_name, " amino acid position"),
      colour_clip = 1.5,
      colour_low = colour_scheme[["shade 0"]][[3]],
      colour_high = colour_scheme[["shade 0"]][[1]])
  }

  ###########################
  ### Binding heatmap
  ###########################

  if("Binding" %in% plot_traits){
    #Heatmap data
    heatmap_dt <- fitness_dt[pca_type=="Binding" & Nham_aa==1,.(fitness = -fitness, fitness_conf = T, Pos_ref, scHAmin_ligand, WT_AA, Mut, Pos_class)]
    #Add dummy mutation to ensure all positions included
    dummy_dt <- fitness_dt[Nham_aa==1 & !is.na(Pos_ref)][!duplicated(Pos_ref) & Nham_aa==1,.(fitness = -fitness, fitness_conf = T, Pos_ref, scHAmin_ligand, WT_AA, Mut, Pos_class)]
    dummy_dt[, Mut := WT_AA]
    dummy_dt[, fitness := NA]
    heatmap_dt <- rbind(heatmap_dt, dummy_dt)

    doubledeepms__plot_heatmap(
      input_dt = heatmap_dt,
      variable_name = "fitness",
      output_file = file.path(outpath, "binding_heatmap.pdf"),
      width = plot_width,
      height = plot_height,
      plot_title = paste0(domain_name, " amino acid position"),
      colour_clip = 1.5,
      colour_low = colour_scheme[["shade 0"]][[3]],
      colour_high = colour_scheme[["shade 0"]][[1]])
  }

  ###########################
  ### Barplot with proportion of mutations significantly affecting binding fitness (due to change in stability/binding free energy)
  ###########################

  if(!"Binding" %in% plot_traits){
    return()
  }

  #Load structure metrics
  dg_dt <- fread(input_file)[id!="-0-",]
  pca_type <- "Binding"
  #Load fitness and metrics data
  rdata_file <- list.files(file.path(input_file_fitness, pca_type))
  load(file.path(input_file_fitness, pca_type, rdata_file))
  #Only up to single mutants
  all_variants <- all_variants[Nham_aa<=1 & STOP==F & STOP_readthrough==F]
  #Single mutant position
  wt_seq <- unlist(strsplit(all_variants[WT==T,aa_seq], ""))
  #Single mutant position class
  all_variants[, Pos := which(unlist(strsplit(aa_seq, "")) != wt_seq), aa_seq]
  all_variants[, WT_AA := wt_seq[Pos], aa_seq]
  all_variants[, Mut := substr(aa_seq, Pos, Pos), aa_seq]
  all_variants[, id := paste0(WT_AA, Pos, Mut), aa_seq]
  singles_dt <- merge(all_variants[Nham_aa==1], dg_dt, by = "id", all.x = T)

  #Significant change in binding fitness
  singles_dt[, fitness_pvalue := doubledeepms__pvalue(fitness, sigma)]
  singles_dt[, fitness_FDR := p.adjust(fitness_pvalue, method = "BH")]

  #Significant change in ddG binding
  singles_dt[b_ddg_pred_conf==T, b_ddg_pred_pvalue := doubledeepms__pvalue(b_ddg_pred, b_ddg_pred_sd)]
  singles_dt[b_ddg_pred_conf==T, b_ddg_pred_FDR := p.adjust(b_ddg_pred_pvalue, method = "BH")]

  #Significant change in ddG folding
  singles_dt[f_ddg_pred_conf==T, f_ddg_pred_pvalue := doubledeepms__pvalue(f_ddg_pred, f_ddg_pred_sd)]
  singles_dt[f_ddg_pred_conf==T, f_ddg_pred_FDR := p.adjust(f_ddg_pred_pvalue, method = "BH")]

  #Subset to significant binding fitness effects and confident ddGs
  singles_dt <- singles_dt[fitness_FDR<0.05 & f_ddg_pred_conf==T & b_ddg_pred_conf==T]

  #Fitness effect
  singles_dt[, fitness_decrease := fitness<0]
  #Fitness decrease
  singles_dt[fitness_decrease==T, ddg_class := "Remainder"]
  singles_dt[fitness_decrease==T & (b_ddg_pred>0 & b_ddg_pred_FDR<0.05) & (f_ddg_pred>0 & f_ddg_pred_FDR<0.05) & abs(b_ddg_pred) > abs(f_ddg_pred), ddg_class := "Binding synergistic"]
  singles_dt[fitness_decrease==T & (b_ddg_pred>0 & b_ddg_pred_FDR<0.05) & (f_ddg_pred>0 & f_ddg_pred_FDR<0.05) & abs(b_ddg_pred) < abs(f_ddg_pred), ddg_class := "Folding synergistic"]
  singles_dt[fitness_decrease==T & (b_ddg_pred<0 & b_ddg_pred_FDR<0.05) & (f_ddg_pred>0 & f_ddg_pred_FDR<0.05), ddg_class := "Folding antagonistic"]
  singles_dt[fitness_decrease==T & (b_ddg_pred>0 & b_ddg_pred_FDR<0.05) & (f_ddg_pred<0 & f_ddg_pred_FDR<0.05), ddg_class := "Binding antagonistic"]
  singles_dt[fitness_decrease==T & (b_ddg_pred_FDR>=0.05) & (f_ddg_pred>0 & f_ddg_pred_FDR<0.05), ddg_class := "Folding only"]
  singles_dt[fitness_decrease==T & (b_ddg_pred>0 & b_ddg_pred_FDR<0.05) & (f_ddg_pred_FDR>=0.05), ddg_class := "Binding only"]
  #Fitness increase
  singles_dt[fitness_decrease==F, ddg_class := "Remainder"]
  singles_dt[fitness_decrease==F & (b_ddg_pred<0 & b_ddg_pred_FDR<0.05) & (f_ddg_pred<0 & f_ddg_pred_FDR<0.05) & abs(b_ddg_pred) > abs(f_ddg_pred), ddg_class := "Binding synergistic"]
  singles_dt[fitness_decrease==F & (b_ddg_pred<0 & b_ddg_pred_FDR<0.05) & (f_ddg_pred<0 & f_ddg_pred_FDR<0.05) & abs(b_ddg_pred) < abs(f_ddg_pred), ddg_class := "Folding synergistic"]
  singles_dt[fitness_decrease==F & (b_ddg_pred>0 & b_ddg_pred_FDR<0.05) & (f_ddg_pred<0 & f_ddg_pred_FDR<0.05), ddg_class := "Folding antagonistic"]
  singles_dt[fitness_decrease==F & (b_ddg_pred<0 & b_ddg_pred_FDR<0.05) & (f_ddg_pred>0 & f_ddg_pred_FDR<0.05), ddg_class := "Binding antagonistic"]
  singles_dt[fitness_decrease==F & (b_ddg_pred_FDR>=0.05) & (f_ddg_pred<0 & f_ddg_pred_FDR<0.05), ddg_class := "Folding only"]
  singles_dt[fitness_decrease==F & (b_ddg_pred<0 & b_ddg_pred_FDR<0.05) & (f_ddg_pred_FDR>=0.05), ddg_class := "Binding only"]

  #Plot - All protein positions together
  plot_dt <- singles_dt[,.(count = as.numeric(.N)), .(fitness_decrease, ddg_class)][order(fitness_decrease)]
  #Calculate percentage
  plot_dt[fitness_decrease==T, percentage := count/sum(count)*100]
  plot_dt[fitness_decrease==F, percentage := count/sum(count)*100]
  plot_dt[, ddg_class := factor(ddg_class, levels = rev(c("Binding only", "Binding antagonistic", "Binding synergistic", "Folding synergistic", "Folding antagonistic", "Folding only", "Remainder")))]
  #Rename binding phenotype
  plot_dt[fitness_decrease==F, binding_phenotype := "Increase"]
  plot_dt[fitness_decrease==T, binding_phenotype := "Decrease"]
  plot_dt[, binding_phenotype := factor(binding_phenotype, levels = c("Increase", "Decrease"))]
  
  plot_cols = c(
    colour_scheme[["shade 3"]][[1]],
    colour_scheme[["shade 0"]][[1]],
    colour_scheme[["shade 1"]][[1]],
    colour_scheme[["shade 1"]][[3]],
    colour_scheme[["shade 0"]][[3]],
    colour_scheme[["shade 3"]][[3]],
     "grey")
  names(plot_cols) <- c("Binding only", "Binding antagonistic", "Binding synergistic", "Folding synergistic", "Folding antagonistic", "Folding only", "Remainder")
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(y = binding_phenotype, percentage, fill = as.factor(ddg_class))) +
    ggplot2::geom_bar(stat="identity") +
    ggplot2::ylab("Binding fitness") +
    ggplot2::xlab("% Mutations") +
    ggplot2::labs(fill = "Biophysical\nmechanism") +
    ggplot2::theme_classic()
  if(!is.null(colour_scheme)){
    d <- d + ggplot2::scale_fill_manual(values = plot_cols)
  }
  ggplot2::ggsave(file.path(outpath, "fitness_binding_mechanism_barplot.pdf"), d, width = 7, height = 3, useDingbats=FALSE)
  #Pleitropic % for binding increase and decrease
  print(paste0("%mutations pleiotropic (increased binding): ", round(plot_dt[fitness_decrease==F & !grepl("only|Remainder", as.character(ddg_class)),sum(count)]/plot_dt[fitness_decrease==F,sum(count)]*100, 2)))
  print(paste0("%mutations pleiotropic (decreased binding): ", round(plot_dt[fitness_decrease==T & !grepl("only|Remainder", as.character(ddg_class)),sum(count)]/plot_dt[fitness_decrease==T,sum(count)]*100, 2)))

  #Plot - split per protein position class
  plot_dt_list <- lapply(unique(singles_dt$Pos_class), function(posclass) {
    temp_plot_dt <- singles_dt[Pos_class == posclass,.(count = as.numeric(.N)), .(fitness_decrease, ddg_class)][order(fitness_decrease)]
    #Calculate percentage
    temp_plot_dt[fitness_decrease==T, percentage := count/sum(count)*100]
    temp_plot_dt[fitness_decrease==F, percentage := count/sum(count)*100]
    temp_plot_dt[, ddg_class := factor(ddg_class, levels = rev(c("Binding only", "Binding antagonistic", "Binding synergistic", "Folding synergistic", "Folding antagonistic", "Folding only", "Remainder")))]
    temp_plot_dt[, Pos_class := posclass]
    temp_plot_dt
  })  
  plot_dt <- do.call("rbind", plot_dt_list)
  
  #Rename binding phenotype
  plot_dt[fitness_decrease==F, binding_phenotype := "Increase"]
  plot_dt[fitness_decrease==T, binding_phenotype := "Decrease"]
  plot_dt[, binding_phenotype := factor(binding_phenotype, levels = c("Increase", "Decrease"))]
  #Rename position class
  plot_dt[, Pos_class_plot := stringr::str_to_title(Pos_class)]
  plot_dt[Pos_class=="binding_interface", Pos_class_plot := "Binding interface"]
  plot_cols = c(
    colour_scheme[["shade 3"]][[1]],
    colour_scheme[["shade 0"]][[1]],
    colour_scheme[["shade 1"]][[1]],
    colour_scheme[["shade 1"]][[3]],
    colour_scheme[["shade 0"]][[3]],
    colour_scheme[["shade 3"]][[3]],
    "grey")
  names(plot_cols) <- c("Binding only", "Binding antagonistic", "Binding synergistic", "Folding synergistic", "Folding antagonistic", "Folding only", "Remainder")
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(y = binding_phenotype, percentage, fill = as.factor(ddg_class))) +
    ggplot2::geom_bar(stat="identity") +
    ggplot2::ylab("Binding fitness") +
    ggplot2::xlab("% Mutations") +
    ggplot2::labs(fill = "Biophysical\nmechanism") +
    ggplot2::theme_classic() +
    ggplot2::facet_grid(Pos_class_plot ~ ., scales = "free")
  if(!is.null(colour_scheme)){
    d <- d + ggplot2::scale_fill_manual(values = plot_cols)
  }
  ggplot2::ggsave(file.path(outpath, "fitness_binding_mechanism_barplot_per_pos_class.pdf"), d, width = 6, height = 5, useDingbats=FALSE)
  
  
  ###########################
  ### Frequency of mutations across multiple sequence alignments
  ###########################

  if (!is.null(input_file_MSA)){
    #Get MSA conservation and frequency
    MSA_dt <- as.data.table(reshape2::melt(fread(file = input_file_MSA), 
                                     id.vars = c("i", "A_i", "conservation"), 
                                     value.name = "Freq_MSA",
                                     variable.name = "AA_sub"))
    #Merge with singles data
    plot_dt <- data.table::merge.data.table(singles_dt, MSA_dt[, .(i, conservation, AA_sub, Freq_MSA)], by.x = c("Pos_ref", "Mut"), by.y = c("i", "AA_sub"))
    #Rename position class
    plot_dt[, Pos_class_plot := stringr::str_to_title(Pos_class)]
    plot_dt[Pos_class=="binding_interface", Pos_class_plot := "Binding interface"]
    plot_dt[, ddg_class := factor(ddg_class, levels = rev(c("Binding only", "Binding antagonistic", "Binding synergistic", "Folding synergistic", "Folding antagonistic", "Folding only", "Remainder")))]
    

    d <- ggplot2::ggplot(plot_dt,ggplot2::aes(y = Freq_MSA, as.factor(ddg_class))) +
      ggplot2::geom_violin(ggplot2::aes(fill = as.factor(ddg_class)), draw_quantiles = c(0.25, 0.5, 0.75), show.legend = F, scale = "width") +
      ggplot2::ylab("Biophisical mechanism") +
      ggplot2::xlab("Frequency AA substitution in MSA") +
      ggplot2::theme_classic() +
      ggplot2::theme(axis.text.x  = ggplot2::element_text(angle = 45, hjust=1)) +
      ggplot2::facet_grid(Pos_class_plot ~ ., scales = "free")
    if(!is.null(colour_scheme)){
      d <- d + ggplot2::scale_fill_manual(values = plot_cols)
    }
    suppressWarnings(ggplot2::ggsave(file.path(outpath, "fitness_binding_mechanism_FreqMSA_violins.pdf"), d, width = 4, height = 5, useDingbats=FALSE))
  }
  
}

