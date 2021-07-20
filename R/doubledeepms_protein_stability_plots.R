
#' doubledeepms_protein_stability_plots
#'
#' Plot free energy heatmaps.
#'
#' @param input_list path to MoCHI thermo model fit results (required)
#' @param pdb_file_list path to PDB file (required)
#' @param pdb_chain_query_list query chain id (required)
#' @param aaprop_file path to amino acid properties file (required)
#' @param aaprop_file_selected path to file with selected subset of identifiers
#' @param input_MSA_list path to MSA frequencies data (required)
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
  input_MSA_list,
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
    temp_dt <- fread(input_list[[protein]])[id!="-0-"]
    temp_dt[, protein := protein]

    #Add WT and mutant AAs
    temp_dt[, WT_AA := substr(id, 1, 1)]
    temp_dt[, Mut := substr(id, nchar(id), nchar(id))]

    #Per residue metrics
    for(i in c("f_ddg", "b_ddg")){
      temp_dt[,paste0(i, "_wposmean") := sum(.SD[[1]]/.SD[[2]]^2, na.rm = T)/sum(1/.SD[[2]]^2, na.rm = T),Pos_ref,.SDcols = paste0(i, c("_pred", "_pred_sd"))]
      temp_dt[,paste0(i, "_wposse") := sqrt(1/sum(1/.SD[[2]]^2, na.rm = T)),Pos_ref,.SDcols = paste0(i, c("_pred", "_pred_sd"))]
    }

    #Per residue metrics - confident only
    for(i in c("f_ddg", "b_ddg")){
      temp_dt[get(paste0(i, "_pred_conf"))==TRUE,paste0(i, "_pred_filtered") := .SD[[1]],,.SDcols = paste0(i, "_pred")]
      temp_dt[get(paste0(i, "_pred_conf"))==TRUE,paste0(i, "_pred_sd_filtered") := .SD[[1]],,.SDcols = paste0(i, "_pred_sd")]
      temp_dt[,paste0(i, "_wposmean_conf") := sum(.SD[[1]]/.SD[[2]]^2, na.rm = T)/sum(1/.SD[[2]]^2, na.rm = T),Pos_ref,.SDcols = paste0(i, c("_pred_filtered", "_pred_sd_filtered"))]
      temp_dt[,paste0(i, "_wposse_conf") := sqrt(1/sum(1/.SD[[2]]^2, na.rm = T)),Pos_ref,.SDcols = paste0(i, c("_pred_filtered", "_pred_sd_filtered"))]
    }

    dg_list[[protein]] <- temp_dt
  }
  dg_dt <- rbindlist(dg_list)

  ###########################
  ### Save mean folding ddG to file
  ###########################

  for(i in names(pdb_file_list)){
    #load PDB structure
    sink(file = "/dev/null")
    pdb <- bio3d::read.pdb(pdb_file_list[[i]], rm.alt = TRUE)
    sink()

    #Replace B factor with mean folding ddG
    pdb_atom_dt <- as.data.table(pdb$atom)
    f_ddg_pred_mean_dt <- dg_dt[protein==i & id!="-0-"][!duplicated(Pos_ref),.(b_new = f_ddg_wposmean, resno = Pos_ref)]
    old_colnames <- names(pdb_atom_dt)
    pdb_atom_dt <- merge(pdb_atom_dt, f_ddg_pred_mean_dt, by = "resno", all.x = T)
    pdb_atom_dt[, b := b_new]
    pdb_atom_dt[is.na(b) | chain!=pdb_chain_query_list[[i]], b := 0]
    pdb$atom <- as.data.frame(pdb_atom_dt[order(eleno),.SD,,.SDcols = old_colnames])

    bio3d::write.pdb(pdb, file = file.path(outpath, gsub(".pdb", "_f_ddg_pred_mean.pdb", basename(pdb_file_list[[i]]))))
  }

  ###########################
  ### Position class violin plots
  ###########################

  plot_dt <- dg_dt[f_ddg_pred_conf==T & id!="-0-",.(protein, f_ddg_pred, Pos_class, RSASA, scHAmin_ligand, Pos_ref)]
  plot_dt[, Pos_class_plot := stringr::str_to_title(Pos_class)]
  plot_dt[Pos_class=="binding_interface", Pos_class_plot := "Binding\ninterface"]
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(Pos_class_plot, f_ddg_pred, fill = Pos_class_plot)) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
    ggplot2::xlab("Residue position") +
    ggplot2::ylab(expression(Delta*Delta*"G of Folding")) +
    ggplot2::facet_grid(~protein, scales = "free") + 
    ggplot2::theme_classic() +
    ggplot2::coord_cartesian(ylim = c(-2, 5)) +
    ggplot2::labs(fill = "Residue\nposition")
  if(!is.null(colour_scheme)){
    d <- d + ggplot2::scale_fill_manual(values = unlist(colour_scheme[["shade 0"]][c(1, 3, 4)]))
  }
  suppressWarnings(ggplot2::ggsave(file.path(outpath, "position_violins_all.pdf"), d, width = 7, height = 3, useDingbats=FALSE))

  plot_dt <- dg_dt[f_ddg_pred_conf==T & id!="-0-" & protein!="GB1",.(protein, f_ddg_pred, Pos_class, RSASA, scHAmin_ligand)]
  plot_dt[, Pos_class_plot := stringr::str_to_title(Pos_class)]
  plot_dt[Pos_class=="binding_interface", Pos_class_plot := "Binding\ninterface"]
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(Pos_class_plot, f_ddg_pred, fill = Pos_class_plot)) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
    ggplot2::xlab("") +
    ggplot2::ylab(expression(Delta*Delta*"G of Folding")) +
    ggplot2::facet_grid(~protein, scales = "free") + 
    ggplot2::theme_classic() +
    ggplot2::coord_cartesian(ylim = c(-1.5, 4)) +
    ggplot2::labs(fill = "Residue\nposition")
  if(!is.null(colour_scheme)){
    d <- d + ggplot2::scale_fill_manual(values = unlist(colour_scheme[["shade 0"]][c(1, 3, 4)]))
  }
  suppressWarnings(ggplot2::ggsave(file.path(outpath, "position_violins.pdf"), d, width = 5, height = 3, useDingbats=FALSE))

  #P-value for each protein separately
  plot_dt <- dg_dt[f_ddg_pred_conf==T & id!="-0-",.(protein, f_ddg_pred, Pos_class, RSASA, scHAmin_ligand, Pos_ref)]
  for(i in plot_dt[,unique(protein)]){
    print(paste0("Folding free energy change of mutations in core residues vs. the remainder, Mann–Whitney U test p-value (", i, "): ", 
    doubledeepms__mann_whitney_U_wrapper(
      plot_dt[protein==i & Pos_class=="core",f_ddg_pred],
      plot_dt[protein==i & Pos_class!="core",f_ddg_pred])['p_value']))
  }

  ###########################
  ### Correlation of mean folding ddG with RSASA
  ###########################

  #Plot correlation with RSASA
  plot_dt <- dg_dt[id!="-0-"][!duplicated(paste(Pos_ref, protein, sep = ":"))]
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(RSASA, f_ddg_wposmean)) +
    ggplot2::geom_point() +
    ggplot2::geom_linerange(ggplot2::aes(ymin = f_ddg_wposmean-f_ddg_wposse*1.96, ymax = f_ddg_wposmean+f_ddg_wposse*1.96)) +
    ggplot2::geom_smooth(method = "lm", formula = 'y~x', color = colour_scheme[["shade 0"]][c(1)], se = F) +
    ggplot2::xlab("Relative solvent accessible surface area (SASA)") +
    ggplot2::ylab(expression("Mean "*Delta*Delta*"G of Folding")) +
    ggplot2::geom_text(data = plot_dt[,.(label = paste("Pearson's r = ", round(cor(RSASA, f_ddg_wposmean, use = "pairwise.complete"), 2), sep="")),.(protein)], ggplot2::aes(label=label, x=Inf, y=Inf, hjust = 1, vjust = 1)) +
    ggplot2::facet_grid(~protein, scales = "free") + 
    ggplot2::theme_bw()
  if(!is.null(colour_scheme)){
    d <- d + ggplot2::scale_color_manual(values = unlist(colour_scheme[["shade 0"]][c(1, 3, 4)]))
  }
  ggplot2::ggsave(file.path(outpath, "mean_ddg_pred_RSASA_scatter_all.pdf"), d, width = 7, height = 3, useDingbats=FALSE)

  #Plot correlation with RSASA
  plot_dt <- dg_dt[f_ddg_pred_conf==T & id!="-0-" & protein!="GB1"][!duplicated(paste(Pos_ref, protein, sep = ":"))]
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(RSASA, f_ddg_wposmean)) +
    ggplot2::geom_point() +
    ggplot2::geom_linerange(ggplot2::aes(ymin = f_ddg_wposmean-f_ddg_wposse*1.96, ymax = f_ddg_wposmean+f_ddg_wposse*1.96)) +
    ggplot2::geom_smooth(method = "lm", formula = 'y~x', color = colour_scheme[["shade 0"]][c(1)], se = F) +
    ggplot2::xlab("Relative solvent accessible surface area (SASA)") +
    ggplot2::ylab(expression("Mean "*Delta*Delta*"G of Folding")) +
    ggplot2::geom_text(data = plot_dt[,.(label = paste("Pearson's r = ", round(cor(RSASA, f_ddg_wposmean, use = "pairwise.complete"), 2), sep="")),.(protein)], ggplot2::aes(label=label, x=Inf, y=Inf, hjust = 1, vjust = 1)) +
    ggplot2::facet_grid(~protein, scales = "free") + 
    ggplot2::theme_bw()
  if(!is.null(colour_scheme)){
    d <- d + ggplot2::scale_color_manual(values = unlist(colour_scheme[["shade 0"]][c(1, 3, 4)]))
  }
  ggplot2::ggsave(file.path(outpath, "mean_ddg_pred_RSASA_scatter.pdf"), d, width = 4, height = 3, useDingbats=FALSE)

  ###########################
  ### Stabilising mutations
  ###########################

  #Significant stabilising mutations
  dg_dt[f_ddg_pred_conf==T & id!="-0-", f_ddg_pred_fdr := p.adjust(doubledeepms__pvalue(f_ddg_pred, f_ddg_pred_sd), method = "BH"),.(protein)]

  #Stabilising mutations
  dg_dt[, f_ddg_pred_stab := 0]
  dg_dt[f_ddg_pred<0 & f_ddg_pred_fdr<0.05, f_ddg_pred_stab := 1]

  #Stabilising residues (at least 5)
  dg_dt[, f_ddg_pred_stab_res5 := F]
  for(i in names(pdb_file_list)){
    temp <- dg_dt[protein==i & f_ddg_pred_stab==1,table(Pos_ref)]
    stab_res <- as.integer(names(temp)[temp>=5])
    if(length(stab_res)==0){stab_res <- as.integer(names(temp))}
    dg_dt[protein==i & Pos_ref %in% stab_res, f_ddg_pred_stab_res5 := T]
    print(paste0("Destabilising residues for ", i, ": ", paste(dg_dt[protein==i & f_ddg_pred_stab_res5][!duplicated(Pos_ref),Pos_ref], collapse = ",")))
    print(paste0("Surface destabilising residues for ", i, ": ", paste(dg_dt[protein==i & f_ddg_pred_stab_res5 & Pos_class=="surface"][!duplicated(Pos_ref),Pos_ref], collapse = ",")))
  }
  ###########################
  ### Position of de-stabilising residues
  ###########################

  #All proteins
  plot_dt <- dg_dt[!duplicated(paste(Pos_ref, protein, sep = ":"))][,.(f_ddg_pred_stab_res5, Pos_class, protein)]
  plot_dt <- plot_dt[,.(count = .N),.(protein, Pos_class, f_ddg_pred_stab_res5)]
  #Calculate percentage
  plot_dt[f_ddg_pred_stab_res5==T & protein=="GB1", percentage := count/sum(count)*100]
  plot_dt[f_ddg_pred_stab_res5==F & protein=="GB1", percentage := count/sum(count)*100]
  plot_dt[f_ddg_pred_stab_res5==T & protein=="GRB2-SH3", percentage := count/sum(count)*100]
  plot_dt[f_ddg_pred_stab_res5==F & protein=="GRB2-SH3", percentage := count/sum(count)*100]
  plot_dt[f_ddg_pred_stab_res5==T & protein=="PSD95-PDZ3", percentage := count/sum(count)*100]
  plot_dt[f_ddg_pred_stab_res5==F & protein=="PSD95-PDZ3", percentage := count/sum(count)*100]
  plot_dt[f_ddg_pred_stab_res5==T, class := "De-stabilising"]
  plot_dt[f_ddg_pred_stab_res5==F, class := "Remainder"]
  plot_dt[, class := factor(class, levels = c("Remainder", "De-stabilising"))]
  
  #All domains
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(y = class, percentage, fill = Pos_class)) +
    ggplot2::geom_bar(stat="identity") +
    ggplot2::geom_text(data = plot_dt[,.(percentage = 50, text = paste0(unlist(count), collapse = ","), Pos_class),.(class, protein)], ggplot2::aes(label = text)) +
    ggplot2::xlab("% Residues") +
    ggplot2::ylab("Residue class") +
    ggplot2::facet_grid(~protein, scales = "free") + 
    ggplot2::theme_classic() +
    ggplot2::labs(fill = "Residue\nposition")
  if(!is.null(colour_scheme)){
    d <- d + ggplot2::scale_fill_manual(values = unlist(colour_scheme[["shade 0"]][c(1, 3, 4)]))
  }
  ggplot2::ggsave(file.path(outpath, "destabilising_position_barplot_all.pdf"), d, width = 7, height = 2, useDingbats=FALSE)

  #GRB2-SH3 and PSD95-PDZ3 only
  d <- ggplot2::ggplot(plot_dt[protein!="GB1"],ggplot2::aes(y = class, percentage, fill = Pos_class)) +
    ggplot2::geom_bar(stat="identity") +
    ggplot2::geom_text(data = plot_dt[protein!="GB1",.(percentage = 50, text = paste0(unlist(count), collapse = ","), Pos_class),.(class, protein)], ggplot2::aes(label = text)) +
    ggplot2::xlab("% Residues") +
    ggplot2::ylab("Residue class") +
    ggplot2::facet_grid(~protein, scales = "free") + 
    ggplot2::theme_classic() +
    ggplot2::labs(fill = "Residue\nposition")
  if(!is.null(colour_scheme)){
    d <- d + ggplot2::scale_fill_manual(values = unlist(colour_scheme[["shade 0"]][c(1, 3, 4)]))
  }
  ggplot2::ggsave(file.path(outpath, "destabilising_position_barplot.pdf"), d, width = 7, height = 2, useDingbats=FALSE)

  ###########################
  ### Amino acid properties of WT residues
  ###########################

  #Amino acid properties PCA scores for each WT residue
  aa_pca_dt <- doubledeepms__aa_properties_pca_WT_loadings(
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
  ### Hydrophobicity of de-stabilising residues
  ###########################

  #All proteins
  plot_dt <- aa_pca_dt[!duplicated(paste(Pos_ref, protein, sep = ":")) & id!="-0-"]
  plot_dt[, Hydrophobicity := .SD[[1]],,.SDcols = "PC1 (Hydrophobicity)"]
  plot_dt[, Hydrophobic_AA := WT_AA %in% unlist(strsplit("AVILMFYW", ""))]
  set.seed(2)
  plot_dt[, Pos_class_plot := Pos_class]
  plot_dt[Pos_class=="surface" & f_ddg_pred_stab_res5, Pos_class_plot := "surface_s"]
  d <- ggplot2::ggplot(plot_dt[Pos_class_plot!="binding_interface"],ggplot2::aes(y = Pos_class_plot, Hydrophobicity, fill = Pos_class_plot)) +
    ggplot2::geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
    ggplot2::geom_jitter(width = 0, height = 0.2, pch = 1) +
    ggplot2::geom_vline(xintercept = 0) +
    ggplot2::facet_grid(~protein, scales = "free") + 
    ggplot2::theme_classic() +
    ggrepel::geom_text_repel(data = plot_dt[Pos_class_plot!="binding_interface" & f_ddg_pred_stab_res5,], ggplot2::aes(label=paste0(WT_AA,Pos_ref))) + 
    ggplot2::labs(fill = "Residue\ntype")
  if(!is.null(colour_scheme)){
    d <- d + ggplot2::scale_fill_manual(values = c(unlist(colour_scheme[["shade 0"]][c(3)]), "grey", unlist(colour_scheme[["shade 0"]][c(4)])))
  }
  ggplot2::ggsave(file.path(outpath, "destabilising_hydrophobicity_violin_all.pdf"), d, width = 8, height = 3, useDingbats=FALSE)

  #GRB2-SH3 and PSD95-PDZ3 only
  set.seed(1)
  d <- ggplot2::ggplot(plot_dt[Pos_class_plot!="binding_interface" & protein!="GB1"],ggplot2::aes(y = Pos_class_plot, Hydrophobicity, fill = Pos_class_plot)) +
    ggplot2::geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
    ggplot2::geom_jitter(width = 0, height = 0.2, pch = 1) +
    ggplot2::geom_vline(xintercept = 0) +
    ggplot2::facet_grid(~protein, scales = "free") + 
    ggplot2::theme_classic() +
    ggplot2::labs(fill = "Residue\ntype")
  if(!is.null(colour_scheme)){
    d <- d + ggplot2::scale_fill_manual(values = c(unlist(colour_scheme[["shade 0"]][c(3)]), "grey", unlist(colour_scheme[["shade 0"]][c(4)])))
  }
  ggplot2::ggsave(file.path(outpath, "destabilising_hydrophobicity_violin.pdf"), d, width = 5, height = 2, useDingbats=FALSE)

  #P-value for each protein separately
  for(i in plot_dt[,unique(protein)]){
    print(paste0("Hydrophobicity of surface destabilising residues vs. remainder, Mann–Whitney U test p-value (", i, "): ", 
    doubledeepms__mann_whitney_U_wrapper(
      plot_dt[protein==i & Pos_class=="surface" & f_ddg_pred_stab_res5,Hydrophobicity],
      plot_dt[protein==i & Pos_class=="surface" & !f_ddg_pred_stab_res5,Hydrophobicity])['p_value']))
  }
  #P-value for PSD95-PDZ3 and GRB2-SH3 pooled
  print(paste0("Hydrophobicity of surface destabilising residues vs. remainder, Mann–Whitney U test p-value (GRB2-SH3 & PSD95-PDZ3): ", 
  doubledeepms__mann_whitney_U_wrapper(
    plot_dt[protein %in% c("GRB2-SH3", "PSD95-PDZ3") & Pos_class=="surface" & f_ddg_pred_stab_res5,Hydrophobicity],
    plot_dt[protein %in% c("GRB2-SH3", "PSD95-PDZ3") & Pos_class=="surface" & !f_ddg_pred_stab_res5,Hydrophobicity])['p_value']))
  
  
  ###########################
  ### Hydrophobicity of destabilizing residues vs. conservation
  ###########################
  
  MSA_list <- list()
  for (domain in names(input_MSA_list)) {
    if(is.null(input_MSA_list[[domain]])){
      msa_dt = data.table::data.table(Pos_ref = plot_dt[protein == domain, Pos_ref])
      msa_dt[, conservation := NA]
      temp_dt <- data.table::merge.data.table(plot_dt[protein == domain, ], 
                                           msa_dt, 
                                           by = "Pos_ref", all.x = T)
    } else {
      msa_dt <- fread(input_MSA_list[[domain]])
      temp_dt <- data.table::merge.data.table(plot_dt[protein == domain, ], 
                                              msa_dt[, .(i, conservation)], 
                                              by.x = "Pos_ref", by.y="i", all.x = T)
    }
    MSA_list[[domain]] <- temp_dt
  }
  MSA_dt<- rbindlist(MSA_list)[!is.na(conservation)]
  
  #Plot
  d <- ggplot2::ggplot(MSA_dt[Pos_class_plot!="binding_interface",], ggplot2::aes(x=Hydrophobicity, y=conservation, color=Pos_class_plot)) +
    ggplot2::geom_point() + 
    ggplot2::geom_smooth(method = "lm", se=FALSE, formula = "y ~ x", linetype=2) +
    ggplot2::facet_grid(Pos_class_plot~protein, scales = "free") +
    ggplot2::theme_classic() +
    ggrepel::geom_label_repel(data = MSA_dt[Pos_class_plot!="binding_interface" & f_ddg_pred_stab_res5,], ggplot2::aes(label=paste0(WT_AA,Pos_ref)))
  if(!is.null(colour_scheme)){
    d <- d + ggplot2::scale_color_manual(values = c(unlist(colour_scheme[["shade 0"]][c(3)]), "grey", unlist(colour_scheme[["shade 0"]][c(4)])))
  }
  ggplot2::ggsave(file.path(outpath, "destabilising_hydrophobicity_vs_conservation.pdf"), d, width = 7, height = 5, useDingbats=FALSE)
  
  
  # ###########################
  # ### Conservation of destabilizing residues
  # ###########################
  
  set.seed(1)
  d <- ggplot2::ggplot(MSA_dt[Pos_class_plot!="binding_interface"],ggplot2::aes(y = Pos_class_plot, conservation, fill = Pos_class_plot)) +
    ggplot2::geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
    ggplot2::geom_jitter(width = 0, height = 0.2, pch = 1) +
    ggrepel::geom_text_repel(data = MSA_dt[Pos_class_plot!="binding_interface" & f_ddg_pred_stab_res5,], ggplot2::aes(label=paste0(WT_AA,Pos_ref))) + 
    ggplot2::facet_grid(~protein, scales = "free") + 
    ggplot2::xlab("Conservation") +
    ggplot2::theme_classic() +
    ggplot2::labs(fill = "Residue\ntype")
  if(!is.null(colour_scheme)){
    d <- d + ggplot2::scale_fill_manual(values = c(unlist(colour_scheme[["shade 0"]][c(3)]), "grey", unlist(colour_scheme[["shade 0"]][c(4)])))
  }
  ggplot2::ggsave(file.path(outpath, "destabilising_conservation_violin.pdf"), d, width = 8, height = 3, useDingbats=FALSE)
  
  # ###########################
  # ### Plot examples quaternary structure destabilizing residues interactions
  # ###########################
  
  # GRB2-SH3
  
  pymol_script = "reinitialize"
  pymol_script[length(pymol_script)+1] = "fetch 1GRI"
  pymol_script[length(pymol_script)+1] = "bg_color white"
  pymol_script[length(pymol_script)+1] = "color paleyellow, chain A"
  pymol_script[length(pymol_script)+1] = "color grey80, chain B"
  pymol_script[length(pymol_script)+1] = "select V53, resi 211 and chain A"
  pymol_script[length(pymol_script)+1] = "show sticks, V53"
  pymol_script[length(pymol_script)+1] = "color atomic,  V53 & (not elem C)"
  pymol_script[length(pymol_script)+1] = "select S31, resi 189 and chain B"
  pymol_script[length(pymol_script)+1] = "show sticks, S31"
  pymol_script[length(pymol_script)+1] = "color atomic,  S31 & (not elem C)"
  pymol_script[length(pymol_script)+1] = "set_view (-0.957919061,-0.034897089,-0.284906834,-0.187024251,0.828841031,0.527293384,0.217741489,0.558390021,-0.800491571,0.000000000,0.000000000, -195.653869629,29.452850342,72.807014465,21.840946198,154.255004883,  237.052734375,  -20.000000000 )"
  pymol_script[length(pymol_script)+1] = "select Y2, resi 160 and chain A"
  pymol_script[length(pymol_script)+1] = "show sticks, Y2"
  pymol_script[length(pymol_script)+1] = "color atomic, Y2 & (not elem C)"
  pymol_script[length(pymol_script)+1] = "select SH2_Glu, resi 87 & chain B"
  pymol_script[length(pymol_script)+1] = "select SH2_Glu, resi 87 & chain B"
  pymol_script[length(pymol_script)+1] = "show sticks, SH2_Glu"
  pymol_script[length(pymol_script)+1] = "color atomic, SH2_Glu & (not elem C)"
  pymol_script[length(pymol_script)+1] = "zoom center, 25"
  pymol_script[length(pymol_script)+1] = "zoom center, 25"
  pymol_script[length(pymol_script)+1] = "ray 2400,2400"
  pymol_script[length(pymol_script)+1] = "png pymol_SH3_example_destabilising_residues_quaternary_struct.png, dpi=600"
  
  write(x = pymol_script, file = file.path(outpath, "GRB2-SH3_example_destabilising_residues_quaternary_struct.txt"))
  
  # PSD95-PDZ3
  pymol_script = "reinitialize"
  pymol_script[length(pymol_script)+1] = "fetch 1be9"
  pymol_script[length(pymol_script)+1] = "bg_color white"
  pymol_script[length(pymol_script)+1] = "hide everything"
  pymol_script[length(pymol_script)+1] = "show cartoon, chain A"
  pymol_script[length(pymol_script)+1] = "color grey80, chain A"
  pymol_script[length(pymol_script)+1] = "select PDZ, resi  311-394"
  pymol_script[length(pymol_script)+1] = "color paleyellow, PDZ"
  pymol_script[length(pymol_script)+1] = "show sticks, chain B"
  pymol_script[length(pymol_script)+1] = "color black, chain B"
  pymol_script[length(pymol_script)+1] = "select Y392, resi 392 and chain A"
  pymol_script[length(pymol_script)+1] = "show sticks, Y392"
  pymol_script[length(pymol_script)+1] = "color atomic,  Y392 & (not elem C)"
  pymol_script[length(pymol_script)+1] = "select P394, resi 394 and chain A"
  pymol_script[length(pymol_script)+1] = "show sticks, P394"
  pymol_script[length(pymol_script)+1] = "color atomic,  Y394 & (not elem C)"
  pymol_script[length(pymol_script)+1] = "set_view (     0.823759854,   -0.005460032,   -0.566911936,    -0.272597969,    0.872961581,   -0.404510409,     0.497101754,    0.487758458,    0.717621565,     0.000000000,    0.000000000, -125.236885071,    40.318359375,   59.448665619,   32.835613251,    98.737716675,  151.736053467,  -20.000000000 )"
  pymol_script[length(pymol_script)+1] = "ray 2400,2400"
  pymol_script[length(pymol_script)+1] = "png pymol_SH3_example_destabilising_residues_quaternary_struct.png, dpi=600"
  
  write(x = pymol_script, file = file.path(outpath, "PSD95-PDZ3_example_destabilising_residues_quaternary_struct.txt"))
  
}

