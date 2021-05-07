
#' doubledeepms_protein_stability_plots
#'
#' Plot free energy heatmaps.
#'
#' @param input_list path to MoCHI thermo model fit results (required)
#' @param pdb_file_list path to PDB file (required)
#' @param pdb_chain_query_list query chain id (required)
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
  pdb_file_list,
  pdb_chain_query_list,
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

  ### Load single mutant free energies
  ###########################

  #Load dg data
  dg_list <- list()
  for(protein in names(input_list)){
    temp_dt <- fread(input_list[[protein]])
    temp_dt[, protein := protein]
    dg_list[[protein]] <- temp_dt
  }
  dg_dt <- rbindlist(dg_list)

  #Set to loop if not helix nor sheet
  dg_dt[is.na(SS), SS := "loop"]

  #Add WT and mutant AAs
  dg_dt[, WT_AA := substr(id, 1, 1)]
  dg_dt[, Mut := substr(id, nchar(id), nchar(id))]

  ###########################
  ### Position class violin plots
  ###########################

  plot_dt <- dg_dt[f_ddg_pred_conf==T & id!="-0-",.(protein, f_ddg_pred, Pos_class)]
  plot_dt[, Pos_class_plot := stringr::str_to_title(Pos_class)]
  plot_dt[Pos_class=="binding_interface", Pos_class_plot := "Binding\ninterface"]
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(Pos_class_plot, f_ddg_pred, fill = Pos_class_plot)) +
    ggplot2::geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
    ggplot2::xlab("Residue position") +
    ggplot2::ylab(expression(Delta*Delta*"G of Folding")) +
    ggplot2::facet_grid(~protein, scales = "free") + 
    ggplot2::theme_bw() +
    ggplot2::coord_cartesian(ylim = c(-2, 5)) +
    ggplot2::labs(fill = "Residue\nposition")
  if(!is.null(colour_scheme)){
    d <- d + ggplot2::scale_fill_manual(values = unlist(colour_scheme[["shade 0"]][c(1, 3, 4)]))
  }
  suppressWarnings(ggplot2::ggsave(file.path(outpath, "position_violins_all.pdf"), d, width = 7, height = 3, useDingbats=FALSE))

  plot_dt <- dg_dt[f_ddg_pred_conf==T & id!="-0-" & protein!="GB1",.(protein, f_ddg_pred, Pos_class)]
  plot_dt[, Pos_class_plot := stringr::str_to_title(Pos_class)]
  plot_dt[Pos_class=="binding_interface", Pos_class_plot := "Binding\ninterface"]
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(Pos_class_plot, f_ddg_pred, fill = Pos_class_plot)) +
    ggplot2::geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
    ggplot2::xlab("") +
    ggplot2::ylab(expression(Delta*Delta*"G of Folding")) +
    ggplot2::facet_grid(~protein, scales = "free") + 
    ggplot2::theme_bw() +
    ggplot2::coord_cartesian(ylim = c(-1.5, 4)) +
    # ggplot2::theme(
    #   axis.text.x = ggplot2::element_blank(),
    #   axis.ticks.x = ggplot2::element_blank()) +
    ggplot2::labs(fill = "Residue\nposition")
  if(!is.null(colour_scheme)){
    d <- d + ggplot2::scale_fill_manual(values = unlist(colour_scheme[["shade 0"]][c(1, 3, 4)]))
  }
  suppressWarnings(ggplot2::ggsave(file.path(outpath, "position_violins.pdf"), d, width = 5, height = 3, useDingbats=FALSE))

  ###########################
  ### Save median folding ddG to file
  ###########################

  for(i in names(pdb_file_list)){
    #load PDB structure
    sink(file = "/dev/null")
    pdb <- bio3d::read.pdb(pdb_file_list[[i]], rm.alt = TRUE)
    sink()

    #Replace B factor with median folding ddG 
    pdb_atom_dt <- as.data.table(pdb$atom)
    f_ddg_pred_med_dt <- dg_dt[protein==i & f_ddg_pred_conf==T & id!="-0-",.(b_new = median(f_ddg_pred), resno = Pos_ref),Pos_ref][,.(resno, b_new)]
    old_colnames <- names(pdb_atom_dt)
    pdb_atom_dt <- merge(pdb_atom_dt, f_ddg_pred_med_dt, by = "resno", all.x = T)
    pdb_atom_dt[, b := b_new]
    pdb_atom_dt[is.na(b) | chain!=pdb_chain_query_list[[i]], b := 0]
    pdb$atom <- as.data.frame(pdb_atom_dt[order(eleno),.SD,,.SDcols = old_colnames])

    bio3d::write.pdb(pdb, file = file.path(outpath, gsub(".pdb", "_f_ddg_pred_median.pdb", basename(pdb_file_list[[i]]))))
  }

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

  ###########################
  ### Stabilising mutations
  ###########################

  for(i in aa_pca_dt[,unique(protein)]){
    aa_pca_dt[protein==i & f_ddg_pred_conf & id!="-0-", f_ddg_pred_fdr := p.adjust(doubledeepms__pvalue(f_ddg_pred, f_ddg_pred_sd), method = "BH")]
  }

  #Stabilising mutations
  aa_pca_dt[, f_ddg_pred_sig := 0]
  aa_pca_dt[f_ddg_pred<0 & f_ddg_pred_fdr<0.05, f_ddg_pred_sig := 1]

  #Add charge to PCs
  aa_pca_dt[, WT_AA_charge := 0]
  aa_pca_dt[grepl("^[RHK]", id), WT_AA_charge := 1]
  aa_pca_dt[grepl("^[DE]", id), WT_AA_charge := -1]
  aa_pca_dt[, Mut_charge := 0]
  aa_pca_dt[grepl("[RHK]$", id), Mut_charge := 1]
  aa_pca_dt[grepl("[DE]$", id), Mut_charge := -1]
  aa_pca_dt[, delta_charge := Mut_charge-WT_AA_charge]
  aa_pca_dt[, Mut_charge_abs := abs(Mut_charge)]
  aa_pca_dt[, WT_AA_charge_abs := abs(WT_AA_charge)]
  aa_pca_dt[, Mut_charge_bin := Mut_charge!=0]
  aa_pca_dt[, WT_AA_charge_bin := WT_AA_charge!=0]

  #Model data
  model_data <- aa_pca_dt[f_ddg_pred_conf==T & protein=="GRB2-SH3",.SD,,.SDcols = grepl("^PC|f_ddg_pred_sig$|RSASA|SS", names(aa_pca_dt))]
  model_data <- model_data[,.SD,,.SDcols = !grepl("PC20", names(model_data))]

  #Model list
  model_list <- list()
  #RSASA model
  model_list[["RSASA"]] <- glm(f_ddg_pred_sig~RSASA, family="binomial", data=model_data)
  #Hydrophobicity model
  model_list[["Hydrophobicity"]] <- glm(f_ddg_pred_sig~`PC1 (Hydrophobicity)`, family="binomial", data=model_data)
  #Full model
  model_list[["Full model"]] <- glm(f_ddg_pred_sig~., family="binomial", data=model_data)

  #Performance
  perf_list <- list()
  for(i in names(model_list)){
    model_predict <- predict(model_list[[i]], data = model_data)
    pred <- ROCR::prediction(model_predict, model_data[,f_ddg_pred_sig])
    perf <- ROCR::performance(pred,"tpr","fpr")
    auc <- round(ROCR::performance(pred, measure = "auc")@'y.values'[[1]], 2)
    #Save
    perf_list[[i]] <- data.table(
      FPR = perf@'x.values'[[1]],
      TPR = perf@'y.values'[[1]],
      measure = i,
      auc = auc)
  }
  plot_dt <- rbindlist(perf_list)
  plot_dt[, measure := factor(measure, levels = c("Hydrophobicity", "RSASA", "Full model"))]
  plot_cols <- c(colour_scheme[["shade 0"]][[2]], colour_scheme[["shade 0"]][[1]], colour_scheme[["shade 0"]][[3]])
  names(plot_cols) <- c("Hydrophobicity", "RSASA", "Full model")

  #Plot
  auc_dt <- plot_dt[!duplicated(measure)][order(measure, decreasing = T)]
  auc_dt[, FPR := 0.5]
  auc_dt[, TPR := seq(0, 1, 1/(.N+1))[2:(.N+1)]]
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(FPR, TPR, color = measure)) +
    ggplot2::geom_line() +
    ggplot2::geom_abline(linetype = 2) +
    ggplot2::xlab("FPR") +
    ggplot2::ylab("TPR") +
    ggplot2::geom_text(data = auc_dt, ggplot2::aes(label=paste("AUC = ", auc, sep=""))) +
    ggplot2::theme_bw() +
    ggplot2::scale_colour_manual(values=plot_cols) +
    ggplot2::labs(color = "Model")   
  ggplot2::ggsave(file.path(outpath, "stabilising_ROC_GRB2.pdf"), d, width = 4.5, height = 3, useDingbats=FALSE)

  # ### Correlate PCs with folding ddGs - stratify by position class
  # ###########################

  # #Correlate PCs with folding ddGs
  # plot_dt <- aa_pca_dt[f_ddg_pred_conf==T,lapply(as.list(.SD), function(x){cor(x, f_ddg_pred)}),.(protein, Pos_class),.SDcols = grepl("^PC", names(aa_pca_dt))]
  # plot_dt <- reshape2::melt(plot_dt, id = c("protein", "Pos_class"))
  # plot_dt[, PC_plot := factor(variable, levels = plot_dt[order(abs(value), decreasing = T),][,rev(unique(variable))])]
  # plot_dt[, Pos_class_plot := stringr::str_to_title(Pos_class)]
  # plot_dt[Pos_class=="binding_interface", Pos_class_plot := "Binding\ninterface"]
  # d <- ggplot2::ggplot(plot_dt,ggplot2::aes(value, PC_plot, fill = Pos_class_plot)) +
  #   ggplot2::geom_bar(stat = "identity", position = 'dodge') +
  #   ggplot2::xlab(expression("Correlation with "*Delta*Delta*"G of Folding")) +
  #   ggplot2::ylab("Amino acid property PC") +
  #   ggplot2::facet_grid(~protein, scales = "free") + 
  #   ggplot2::theme_classic() +
  #   ggplot2::theme(axis.text.x=ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  #   ggplot2::labs(fill = "Residue\nposition")
  # if(!is.null(colour_scheme)){
  #   d <- d + ggplot2::scale_fill_manual(values = unlist(colour_scheme[["shade 0"]][c(1, 3, 4)]))
  # }
  # ggplot2::ggsave(file.path(outpath, "folding_aa_prop_barplot.pdf"), d, width = 10, height = 7, useDingbats=FALSE)

  # ### Correlate PCs with folding ddGs - stratify by secondary structure element
  # ###########################

  # #Correlate PCs with folding ddGs
  # plot_dt <- aa_pca_dt[f_ddg_pred_conf==T,lapply(as.list(.SD), function(x){cor(x, f_ddg_pred)}),.(protein, SS),.SDcols = grepl("^PC", names(aa_pca_dt))]
  # plot_dt <- reshape2::melt(plot_dt, id = c("protein", "SS"))
  # plot_dt[, PC_plot := factor(variable, levels = plot_dt[order(abs(value), decreasing = T),][,rev(unique(variable))])]
  # plot_dt[, SS_plot := stringr::str_to_title(SS)]
  # d <- ggplot2::ggplot(plot_dt,ggplot2::aes(value, PC_plot, fill = SS_plot)) +
  #   ggplot2::geom_bar(stat = "identity", position = 'dodge') +
  #   ggplot2::xlab(expression("Correlation with "*Delta*Delta*"G of Folding")) +
  #   ggplot2::ylab("Amino acid property PC") +
  #   ggplot2::facet_grid(~protein, scales = "free") + 
  #   ggplot2::theme_classic() +
  #   ggplot2::theme(axis.text.x=ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  #   ggplot2::labs(fill = "Residue\nposition")
  # if(!is.null(colour_scheme)){
  #   d <- d + ggplot2::scale_fill_manual(values = unlist(colour_scheme[["shade 0"]][c(1, 3, 4)]))
  # }
  # ggplot2::ggsave(file.path(outpath, "folding_aa_prop_barplot_SS.pdf"), d, width = 10, height = 7, useDingbats=FALSE)

  # ### Correlate PCs with folding ddGs - stratify by secondary structure element and position class and pool proteins
  # ###########################

  # #Correlate PCs with folding ddGs
  # plot_dt <- aa_pca_dt[f_ddg_pred_conf==T,lapply(as.list(.SD), function(x){cor(x, f_ddg_pred)}),.(Pos_class, SS),.SDcols = grepl("^PC", names(aa_pca_dt))]
  # plot_dt <- reshape2::melt(plot_dt, id = c("Pos_class", "SS"))
  # plot_dt[, PC_plot := factor(variable, levels = plot_dt[order(abs(value), decreasing = T),][,rev(unique(variable))])]
  # plot_dt[, SS_plot := stringr::str_to_title(SS)]
  # d <- ggplot2::ggplot(plot_dt,ggplot2::aes(value, PC_plot, fill = SS_plot)) +
  #   ggplot2::geom_bar(stat = "identity", position = 'dodge') +
  #   ggplot2::xlab(expression("Correlation with "*Delta*Delta*"G of Folding")) +
  #   ggplot2::ylab("Amino acid property PC") +
  #   ggplot2::facet_grid(~Pos_class, scales = "free") + 
  #   ggplot2::theme_classic() +
  #   ggplot2::theme(axis.text.x=ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  #   ggplot2::labs(fill = "Residue\nposition")
  # if(!is.null(colour_scheme)){
  #   d <- d + ggplot2::scale_fill_manual(values = unlist(colour_scheme[["shade 0"]][c(1, 3, 4)]))
  # }
  # ggplot2::ggsave(file.path(outpath, "folding_aa_prop_barplot_SS_Pos_class.pdf"), d, width = 10, height = 7, useDingbats=FALSE)

  # ###########################
  # ### Hydrophobicity scatterplot
  # ###########################

  # plot_dt <- aa_pca_dt[f_ddg_pred_conf==T,.(protein, f_ddg_pred, Hydrophobicity = .SD[[1]], Pos_class),,.SDcols = c("PC1 (Hydrophobicity)")]
  # d <- ggplot2::ggplot(plot_dt,ggplot2::aes(Hydrophobicity, f_ddg_pred)) +
  #   ggplot2::stat_binhex(bins = 20, size = 0.2, color = "grey") +
  #   ggplot2::scale_fill_gradientn(colours = c("white", "black")) +
  #   ggplot2::geom_smooth(method = "lm", color = colour_scheme[["shade 0"]][1]) +
  #   ggplot2::xlab(expression(Delta*" Hydrophobicity")) +
  #   ggplot2::ylab(expression(Delta*Delta*"G of Folding")) +
  #   ggplot2::facet_grid(protein~Pos_class, scales = "free") + 
  #   ggplot2::geom_text(data = plot_dt[,.(label = paste("Pearson's r = ", round(cor(Hydrophobicity, f_ddg_pred, use = "pairwise.complete"), 2), sep="")),.(protein, Pos_class)], ggplot2::aes(label=label, x=-Inf, y=Inf, hjust = 0, vjust = 1)) +
  #   ggplot2::theme_bw() +
  #   ggplot2::labs(color = "Model")   
  # ggplot2::ggsave(file.path(outpath, "hydrophobicity_scatter_all.pdf"), d, width = 7, height = 7, useDingbats=FALSE)

  # plot_dt <- aa_pca_dt[f_ddg_pred_conf==T & Pos_class!="binding_interface" & protein!="GB1",.(protein, f_ddg_pred, Hydrophobicity = .SD[[1]], Pos_class),,.SDcols = c("PC1 (Hydrophobicity)")]
  # d <- ggplot2::ggplot(plot_dt,ggplot2::aes(Hydrophobicity, f_ddg_pred)) +
  #   ggplot2::stat_binhex(bins = 20, size = 0.2, color = "grey") +
  #   ggplot2::scale_fill_gradientn(colours = c("white", "black")) +
  #   ggplot2::geom_smooth(method = "lm", color = colour_scheme[["shade 0"]][1]) +
  #   ggplot2::xlab(expression(Delta*" Hydrophobicity")) +
  #   ggplot2::ylab(expression(Delta*Delta*"G of Folding")) +
  #   ggplot2::facet_grid(protein~Pos_class, scales = "free") + 
  #   ggplot2::geom_text(data = plot_dt[,.(label = paste("r = ", round(cor(Hydrophobicity, f_ddg_pred, use = "pairwise.complete"), 2), sep="")),.(protein, Pos_class)], ggplot2::aes(label=label, x=-Inf, y=Inf, hjust = 0, vjust = 1)) +
  #   ggplot2::theme_bw() +
  #   ggplot2::labs(color = "Model")   
  # ggplot2::ggsave(file.path(outpath, "hydrophobicity_scatter.pdf"), d, width = 4, height = 3, useDingbats=FALSE)

  # ###########################
  # ### RSASA scatterplot
  # ###########################

  # plot_dt <- aa_pca_dt[f_ddg_pred_conf==T,.(protein, f_ddg_pred, RSASA, Pos_class)]
  # d <- ggplot2::ggplot(plot_dt,ggplot2::aes(RSASA, f_ddg_pred)) +
  #   ggplot2::stat_binhex(bins = 30, size = 0.2, color = "grey") +
  #   ggplot2::scale_fill_gradientn(colours = c("white", "black")) +
  #   ggplot2::geom_smooth(method = "lm", color = colour_scheme[["shade 0"]][1]) +
  #   ggplot2::xlab("RSASA") +
  #   ggplot2::ylab(expression(Delta*Delta*"G of Folding")) +
  #   ggplot2::facet_grid(~protein, scales = "free") + 
  #   ggplot2::geom_text(data = plot_dt[,.(label = paste("r = ", round(cor(RSASA, f_ddg_pred, use = "pairwise.complete"), 2), sep="")),.(protein)], ggplot2::aes(label=label, x=-Inf, y=Inf, hjust = 0, vjust = 1)) +
  #   ggplot2::theme_bw() +
  #   ggplot2::labs(color = "Model")   
  # ggplot2::ggsave(file.path(outpath, "RSASA_scatter.pdf"), d, width = 7, height = 3, useDingbats=FALSE)

  # ###########################
  # ### Folding ROC plot
  # ###########################

  # #Destabilising mutations
  # aa_pca_dt[f_ddg_pred_conf==T, f_ddg_pred_sig := f_ddg_pred>2]

  # #Model data
  # model_data <- aa_pca_dt[f_ddg_pred_conf==T & protein!="GB1",.SD,,.SDcols = grepl("^PC|f_ddg_pred_sig$|RSASA", names(aa_pca_dt))]

  # #Model list
  # model_list <- list()
  # #RSASA model
  # model_list[["RSASA"]] <- glm(f_ddg_pred_sig~RSASA, family="binomial", data=model_data)
  # #Hydrophobicity model
  # model_list[["Hydrophobicity"]] <- glm(f_ddg_pred_sig~`PC1 (Hydrophobicity)`, family="binomial", data=model_data)
  # #Full model
  # model_list[["Full model"]] <- glm(f_ddg_pred_sig~., family="binomial", data=model_data)

  # #Performance
  # perf_list <- list()
  # for(i in names(model_list)){
  #   model_predict <- predict(model_list[[i]], data = model_data)
  #   pred <- ROCR::prediction(model_predict, model_data[,f_ddg_pred_sig])
  #   perf <- ROCR::performance(pred,"tpr","fpr")
  #   auc <- round(ROCR::performance(pred, measure = "auc")@'y.values'[[1]], 2)
  #   #Save
  #   perf_list[[i]] <- data.table(
  #     FPR = perf@'x.values'[[1]],
  #     TPR = perf@'y.values'[[1]],
  #     measure = i,
  #     auc = auc)
  # }
  # plot_dt <- rbindlist(perf_list)
  # plot_dt[, measure := factor(measure, levels = c("Hydrophobicity", "RSASA", "Full model"))]
  # plot_cols <- c(colour_scheme[["shade 0"]][[2]], colour_scheme[["shade 0"]][[1]], colour_scheme[["shade 0"]][[3]])
  # names(plot_cols) <- c("Hydrophobicity", "RSASA", "Full model")

  # #Plot
  # auc_dt <- plot_dt[!duplicated(measure)][order(measure, decreasing = T)]
  # auc_dt[, FPR := 0.5]
  # auc_dt[, TPR := seq(0, 1, 1/(.N+1))[2:(.N+1)]]
  # d <- ggplot2::ggplot(plot_dt,ggplot2::aes(FPR, TPR, color = measure)) +
  #   ggplot2::geom_line() +
  #   ggplot2::geom_abline(linetype = 2) +
  #   ggplot2::xlab("FPR") +
  #   ggplot2::ylab("TPR") +
  #   ggplot2::geom_text(data = auc_dt, ggplot2::aes(label=paste("AUC = ", auc, sep=""))) +
  #   ggplot2::theme_bw() +
  #   ggplot2::scale_colour_manual(values=plot_cols) +
  #   ggplot2::labs(color = "Model")   
  # ggplot2::ggsave(file.path(outpath, "folding_ROC.pdf"), d, width = 4.5, height = 3, useDingbats=FALSE)

  # ###########################
  # ### Folding ROC plot - GB1
  # ###########################

  # #Destabilising mutations
  # aa_pca_dt[f_ddg_pred_conf==T, f_ddg_pred_sig := f_ddg_pred>2]

  # #Model data
  # model_data <- aa_pca_dt[f_ddg_pred_conf==T & protein=="GB1",.SD,,.SDcols = grepl("^PC|f_ddg_pred_sig$|RSASA", names(aa_pca_dt))]

  # #Model list
  # model_list <- list()
  # #RSASA model
  # model_list[["RSASA"]] <- glm(f_ddg_pred_sig~RSASA, family="binomial", data=model_data)
  # #Hydrophobicity model
  # model_list[["Hydrophobicity"]] <- glm(f_ddg_pred_sig~`PC1 (Hydrophobicity)`, family="binomial", data=model_data)
  # #Full model
  # model_list[["Full model"]] <- glm(f_ddg_pred_sig~., family="binomial", data=model_data)

  # #Performance
  # perf_list <- list()
  # for(i in names(model_list)){
  #   model_predict <- predict(model_list[[i]], data = model_data)
  #   pred <- ROCR::prediction(model_predict, model_data[,f_ddg_pred_sig])
  #   perf <- ROCR::performance(pred,"tpr","fpr")
  #   auc <- round(ROCR::performance(pred, measure = "auc")@'y.values'[[1]], 2)
  #   #Save
  #   perf_list[[i]] <- data.table(
  #     FPR = perf@'x.values'[[1]],
  #     TPR = perf@'y.values'[[1]],
  #     measure = i,
  #     auc = auc)
  # }
  # plot_dt <- rbindlist(perf_list)
  # plot_dt[, measure := factor(measure, levels = c("Hydrophobicity", "RSASA", "Full model"))]
  # plot_cols <- c(colour_scheme[["shade 0"]][[2]], colour_scheme[["shade 0"]][[1]], colour_scheme[["shade 0"]][[3]])
  # names(plot_cols) <- c("Hydrophobicity", "RSASA", "Full model")

  # #Plot
  # auc_dt <- plot_dt[!duplicated(measure)][order(measure, decreasing = T)]
  # auc_dt[, FPR := 0.5]
  # auc_dt[, TPR := seq(0, 1, 1/(.N+1))[2:(.N+1)]]
  # d <- ggplot2::ggplot(plot_dt,ggplot2::aes(FPR, TPR, color = measure)) +
  #   ggplot2::geom_line() +
  #   ggplot2::geom_abline(linetype = 2) +
  #   ggplot2::xlab("FPR") +
  #   ggplot2::ylab("TPR") +
  #   ggplot2::geom_text(data = auc_dt, ggplot2::aes(label=paste("AUC = ", auc, sep=""))) +
  #   ggplot2::theme_bw() +
  #   ggplot2::scale_colour_manual(values=plot_cols) +
  #   ggplot2::labs(color = "Model")   
  # ggplot2::ggsave(file.path(outpath, "folding_ROC_GB1.pdf"), d, width = 4.5, height = 3, useDingbats=FALSE)

}

