
#' doubledeepms_protein_stability_plots
#'
#' Plot free energy heatmaps.
#'
#' @param input_list path to MoCHI thermo model fit results (required)
#' @param aaprop_file path to amino acid properties file (required)
#' @param aaprop_file_selected path to file with selected subset of identifiers
#' @param outpath output path for plots and saved objects (required)
#' @param colour_scheme colour scheme file (required)
#' @param execute whether or not to execute the analysis (default: TRUE)
#'
#' @return Nothing
#' @export
#' @import data.table
doubledeepms_protein_stability_plots <- function(
  input_list,
  aaprop_file,
  aaprop_file_selected,
  outpath,
  colour_scheme,
  execute = TRUE
  ){

  #Return if analysis not executed
  if(!execute){
    return()
  }

  #Display status
  message(paste("\n\n*******", "running stage: doubledeepms_protein_stability_plots", "*******\n\n"))

  #Create output directory
  doubledeepms__create_dir(doubledeepms_dir = outpath)

  ### AA properties PCA
  ###########################

  #AA properties PCA
  exp_pca_list <- doubledeepms__aa_properties_pca(
    aa_properties_file = aaprop_file, 
    selected_identifiers = unlist(fread(aaprop_file_selected, header = F)),
    return_evidences = T)
  exp_pca <- exp_pca_list[["PCA"]]
  aa_evidences <- exp_pca_list[["evidences"]]

  #% variance explained by top 5 PCs
  top5pc_var <- sum((exp_pca$sdev^2/sum(exp_pca$sdev^2))[1:5])

  #Screeplot
  plot_df <- data.frame(perc_var = exp_pca$sdev^2/sum(exp_pca$sdev^2)*100, pc = 1:20)
  plot_df[,"pc"] <- factor(plot_df[,"pc"], levels = 1:20)
  d <- ggplot2::ggplot(plot_df, ggplot2::aes(pc, perc_var)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::geom_vline(xintercept = 5.5, linetype = 2) +
    ggplot2::theme_bw() +
    ggplot2::annotate("text", x = 10, y = 30, label = paste0("Var. explained by PC1-5 = ", round(top5pc_var*100, 0), "%"))
  ggplot2::ggsave(file=file.path(outpath, 'PCA_screeplot.pdf'), width=4, height=4, useDingbats=FALSE)

  #Top features on top 5 PCs
  aa_evidences_name <- as.list(paste(names(aa_evidences), unlist(aa_evidences), sep = ": "))
  names(aa_evidences_name) <- names(aa_evidences)
  for(i in 1:20){
    output_file = file.path(outpath, paste0("PCA_loadings_PC", i, "_highlow.txt"))
    write.table(data.frame(), file=output_file, col.names=FALSE)
    lapply(aa_evidences_name[rownames(exp_pca$rotation)[order(exp_pca$rotation[,i], decreasing = T)[1:20]]], write, output_file, append=TRUE, ncolumns=1000)
    write("...", file = output_file, append=TRUE)
    lapply(rev(aa_evidences_name[rownames(exp_pca$rotation)[order(exp_pca$rotation[,i], decreasing = F)[1:20]]]), write, output_file, append=TRUE, ncolumns=1000)
  }

  #Feature type
  temp_cols <- c("black", unlist(colour_scheme[["shade 0"]][1:4]), "grey")
  feature_type <- rep("6_Remainder", dim(exp_pca$rotation)[1])
  feature_type[grep("hydrophobic|Hydrophobic", unlist(aa_evidences))] <- "1_hydrophobic/Hydrophobic"
  feature_type[grep("helix|helical", unlist(aa_evidences))] <- "2_helix/helical"
  feature_type[grep("composition|Composition", unlist(aa_evidences))] <- "3_composition/Composition"
  feature_type[grep("linker|Linker", unlist(aa_evidences))] <- "4_linker/Linker"
  feature_type[grep("beta-sheet|beta-strand|Beta-sheet|Beta-strand", unlist(aa_evidences))] <- "5_beta-sheet/beta-strand/Beta-sheet/Beta-strand"
  names(temp_cols) <- unique(feature_type[order(feature_type)])

  #Plot PCA variable loadings
  doubledeepms__plot_loadings(
    pca_obj=exp_pca, 
    output_file=file.path(outpath, paste0('PCA_biplots_5_symbols.pdf')),
    plot_categories=feature_type,
    plot_colours=temp_cols, 
    comps=1:5, 
    plot_symbol=19,
    width=15, height=15)

  #Plot PCA variable loadings legend
  doubledeepms__plot_loadings(
    pca_obj=exp_pca, 
    output_file=file.path(outpath, paste0('PCA_biplots_5_symbols_PC1_PC2.pdf')),
    plot_categories=feature_type,
    plot_colours=temp_cols, 
    comps=1:2, 
    plot_symbol=19,
    width=5, height=5)

  ### Single mutant loadings
  ###########################

  #Load dg data
  dg_list <- list()
  for(protein in names(input_list)){
    temp_dt <- fread(input_list[[protein]])
    temp_dt[, protein := protein]
    dg_list[[protein]] <- temp_dt
  }
  dg_dt <- rbindlist(dg_list)

  #Add WT and mutant AAs
  dg_dt[, WT_AA := substr(id, 1, 1)]
  dg_dt[, Mut := substr(id, nchar(id), nchar(id))]

  #Amino acid properties PCA scores for each single mutant
  aa_pca_dt <- doubledeepms__aa_properties_pca_singles_loadings(
    input_dt = dg_dt[id!="-0-"],
    aa_properties_file = aaprop_file,
    selected_identifiers = unlist(fread(aaprop_file_selected, header = F)))

  #Rename PCs
  names(aa_pca_dt) <- gsub("_score", "", names(aa_pca_dt))
  top_PCs <- 5
  top_PC_signs <- c(-1, -1, 1, -1, -1)
  top_PC_names <- c("Hydrophobicity", "Helix propensity", "Commonness", "Linker propensity", "Beta-sheet propensity")
  for(i in 1:top_PCs){aa_pca_dt[, (paste0('PC', i)) := scale(.SD, scale=top_PC_signs[i], center=F),,.SDcols = paste0('PC', i)]}
  names(aa_pca_dt)[grep("^PC", names(aa_pca_dt))][1:5] <- paste0(names(aa_pca_dt)[grep("^PC", names(aa_pca_dt))][1:5], " (", top_PC_names, ")")

  #Correlate PCs with folding ddGs
  cor_list <- list()
  for(i in aa_pca_dt[,unique(protein)]){
    cor_core <- aa_pca_dt[Pos_class=="core" & protein==i & f_ddg_pred_conf==T,cor(.SD)[-1,1],,.SDcols = grepl("^PC|f_ddg_pred$", names(aa_pca_dt))]
    cor_surf <- aa_pca_dt[Pos_class=="surface" & protein==i & f_ddg_pred_conf==T,cor(.SD)[-1,1],,.SDcols = grepl("^PC|f_ddg_pred$", names(aa_pca_dt))]
    cor_bind <- aa_pca_dt[Pos_class=="binding_interface" & protein==i & f_ddg_pred_conf==T,cor(.SD)[-1,1],,.SDcols = grepl("^PC|f_ddg_pred$", names(aa_pca_dt))]
    cor_list[[i]] <- copy(rbind(
      data.table(cor = cor_core, Pos_class = "core", PC = names(aa_pca_dt)[grep("^PC", names(aa_pca_dt))]),
      data.table(cor = cor_surf, Pos_class = "surface", PC = names(aa_pca_dt)[grep("^PC", names(aa_pca_dt))]),
      data.table(cor = cor_bind, Pos_class = "binding_interface", PC = names(aa_pca_dt)[grep("^PC", names(aa_pca_dt))])))
    cor_list[[i]][, protein := i]
  }

  #Plot correlations for all proteins and position classes
  plot_dt <- rbindlist(cor_list)
  plot_dt[, PC_plot := factor(PC, levels = plot_dt[order(abs(cor), decreasing = T),][,rev(unique(PC))])]
  plot_dt[, Pos_class_plot := stringr::str_to_title(Pos_class)]
  plot_dt[Pos_class=="binding_interface", Pos_class_plot := "Binding\ninterface"]
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(cor, PC_plot, fill = Pos_class_plot)) +
    ggplot2::geom_bar(stat = "identity", position = 'dodge') +
    ggplot2::xlab(expression("Correlation with "*Delta*Delta*"G of Folding")) +
    ggplot2::ylab("Amino acid property PC") +
    ggplot2::facet_grid(~protein, scales = "free") + 
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x=ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    ggplot2::labs(fill = "Residue\nposition")
  if(!is.null(colour_scheme)){
    d <- d + ggplot2::scale_fill_manual(values = unlist(colour_scheme[["shade 0"]][c(1, 3, 4)]))
  }
  ggplot2::ggsave(file.path(outpath, "folding_aa_prop_barplot.pdf"), d, width = 10, height = 7, useDingbats=FALSE)

}

