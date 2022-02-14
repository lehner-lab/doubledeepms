
#' doubledeepms_allostery_plots
#'
#' Plot free energy heatmaps.
#'
#' @param input_file path to input file (required)
#' @param pdb_file_list path to PDB file (required)
#' @param pdb_chain_query_list query chain id (required)
#' @param annotation_list annotations of allosteric sites (required)
#' @param ohm_file_list ohm output file list (required)
#' @param outpath output path for plots and saved objects (required)
#' @param colour_scheme colour scheme file (required)
#' @param execute whether or not to execute the analysis (default: TRUE)
#'
#' @return Nothing
#' @export
#' @import data.table
doubledeepms_allostery_plots <- function(
  input_file,
  pdb_file_list,
  pdb_chain_query_list,
  annotation_list,
  ohm_file_list,
  outpath,
  colour_scheme,
  execute = TRUE
  ){

  #Return if analysis not executed
  if(!execute){
    return()
  }

  #Display status
  message(paste("\n\n*******", "running stage: doubledeepms_allostery_plots", "*******\n\n"))

  #Create output directory
  doubledeepms__create_dir(doubledeepms_dir = outpath)

  ### Load single mutant free energies
  ###########################

  #Load dg data
  dg_dt <- fread(input_file)

  ###########################
  ### Literature sites
  ###########################

  literature_list <- list()
  for(protein in names(pdb_file_list)){
    literature_list[[protein]] <- list(class_switching = c(), sector = c())
    if(protein %in% names(annotation_list)){
      anno_dt <- fread(annotation_list[[protein]])
      literature_list[[protein]][["class_switching"]] <- anno_dt[get("CS;Mclaughlin2012")==1,Pos_ref]
      literature_list[[protein]][["sector"]] <- anno_dt[get("SCA;Mclaughlin2012")==1,Pos_ref]
    }
  }

  ###########################
  ### Save mean binding ddG to file
  ###########################

  for(i in names(pdb_file_list)){
    #load PDB structure
    sink(file = "/dev/null")
    pdb <- bio3d::read.pdb(pdb_file_list[[i]], rm.alt = TRUE)
    sink()

    #Replace B factor with mean folding ddG 
    pdb_atom_dt <- as.data.table(pdb$atom)
    b_ddg_pred_mean_dt <- dg_dt[protein==i & id!="-0-"][!duplicated(Pos_ref),.(b_new = b_ddg_wposmeanabs, resno = Pos_ref)]
    old_colnames <- names(pdb_atom_dt)
    pdb_atom_dt <- merge(pdb_atom_dt, b_ddg_pred_mean_dt, by = "resno", all.x = T)
    pdb_atom_dt[, b := b_new]
    pdb_atom_dt[is.na(b) | chain!=pdb_chain_query_list[[i]], b := 0]
    pdb$atom <- as.data.frame(pdb_atom_dt[order(eleno),.SD,,.SDcols = old_colnames])

    bio3d::write.pdb(pdb, file = file.path(outpath, gsub(".pdb", "_b_ddg_pred_abs_mean.pdb", basename(pdb_file_list[[i]]))))
  }

  ###########################
  ### Folding energy distance correlation plots
  ###########################

  for(i in dg_dt[,unique(protein)]){
    doubledeepms__persite_energy_vs_distance_plot(
      input_dt = copy(dg_dt)[protein==i & id!="-0-"],
      outpath = file.path(outpath, paste0(i, "_persite_folding_energy_vs_distance_scatter.pdf")),
      colour_scheme = colour_scheme,
      trait_name = "folding")
  }

  ###########################
  ### Binding energy distance correlation plots (no absolute value)
  ###########################

  for(i in dg_dt[,unique(protein)]){
    allostery_pos <- doubledeepms__persite_energy_vs_distance_plot(
      input_dt = copy(dg_dt)[protein==i & id!="-0-"],
      literature_sites = literature_list[[i]][["class_switching"]],
      outpath = file.path(outpath, paste0(i, "_persite_binding_energy_vs_distance_scatter_noabs.pdf")),
      colour_scheme = colour_scheme,
      trait_name = "binding",
      absolute_value = F)
    #Limit y axis to sites with weighted mean ddG<3
    if(dg_dt[,max(b_ddg_wposmeanabs)]>3){
      doubledeepms__persite_energy_vs_distance_plot(
        input_dt = copy(dg_dt)[protein==i & id!="-0-" & b_ddg_wposmeanabs<3],
        literature_sites = literature_list[[i]][["class_switching"]],
        outpath = file.path(outpath, paste0(i, "_persite_binding_energy_vs_distance_scatter_noabs_l3.pdf")),
        colour_scheme = colour_scheme,
        trait_name = "binding",
        absolute_value = F)
    }
    dg_dt[protein==i, allosteric := NA]
    dg_dt[protein==i, orthosteric := NA]
    dg_dt[protein==i & Pos_ref %in% allostery_pos & Pos_class!="binding_interface", allosteric := T]
    dg_dt[protein==i & Pos_ref %in% allostery_pos & Pos_class=="binding_interface", orthosteric := T]
    print(paste0("Allosteric residues for ", i, " (no absolute value): ", paste(dg_dt[protein==i & allosteric][!duplicated(Pos_ref),Pos_ref], collapse = ",")))
    print(paste0("Orthosteric residues for ", i, " (no absolute value): ", paste(dg_dt[protein==i & orthosteric][!duplicated(Pos_ref),Pos_ref], collapse = ",")))
  }

  ###########################
  ### Binding energy distance correlation plots (distal threshold residues)
  ###########################

  for(i in dg_dt[,unique(protein)]){
    allostery_pos <- doubledeepms__persite_energy_vs_distance_plot(
      input_dt = copy(dg_dt)[protein==i & id!="-0-"],
      literature_sites = literature_list[[i]][["class_switching"]],
      outpath = file.path(outpath, paste0(i, "_persite_binding_energy_vs_distance_scatter_distaltreshold.pdf")),
      colour_scheme = colour_scheme,
      trait_name = "binding",
      threshold_residues = "distal")
    #Limit y axis to sites with weighted mean ddG<3
    if(dg_dt[,max(b_ddg_wposmeanabs)]>3){
      doubledeepms__persite_energy_vs_distance_plot(
        input_dt = copy(dg_dt)[protein==i & id!="-0-" & b_ddg_wposmeanabs<3],
        literature_sites = literature_list[[i]][["class_switching"]],
        outpath = file.path(outpath, paste0(i, "_persite_binding_energy_vs_distance_scatter_distaltreshold_l3.pdf")),
        colour_scheme = colour_scheme,
        trait_name = "binding",
        threshold_residues = "distal")
    }
    dg_dt[protein==i, allosteric := NA]
    dg_dt[protein==i, orthosteric := NA]
    dg_dt[protein==i & Pos_ref %in% allostery_pos & Pos_class!="binding_interface", allosteric := T]
    dg_dt[protein==i & Pos_ref %in% allostery_pos & Pos_class=="binding_interface", orthosteric := T]
    print(paste0("Allosteric residues for ", i, " (distal threshold): ", paste(dg_dt[protein==i & allosteric][!duplicated(Pos_ref),Pos_ref], collapse = ",")))
    print(paste0("Orthosteric residues for ", i, " (distal threshold): ", paste(dg_dt[protein==i & orthosteric][!duplicated(Pos_ref),Pos_ref], collapse = ",")))
  }

  ###########################
  ### Binding energy distance correlation plots
  ###########################

  for(i in dg_dt[,unique(protein)]){
    allostery_pos <- doubledeepms__persite_energy_vs_distance_plot(
      input_dt = copy(dg_dt)[protein==i & id!="-0-"],
      literature_sites = literature_list[[i]][["class_switching"]],
      outpath = file.path(outpath, paste0(i, "_persite_binding_energy_vs_distance_scatter.pdf")),
      colour_scheme = colour_scheme,
      trait_name = "binding")
    #Limit y axis to sites with weighted mean ddG<3
    if(dg_dt[,max(b_ddg_wposmeanabs)]>3){
      doubledeepms__persite_energy_vs_distance_plot(
        input_dt = copy(dg_dt)[protein==i & id!="-0-" & b_ddg_wposmeanabs<3],
        literature_sites = literature_list[[i]][["class_switching"]],
        outpath = file.path(outpath, paste0(i, "_persite_binding_energy_vs_distance_scatter_l3.pdf")),
        colour_scheme = colour_scheme,
        trait_name = "binding")
    }
    dg_dt[protein==i, allosteric := NA]
    dg_dt[protein==i, orthosteric := NA]
    dg_dt[protein==i & Pos_ref %in% allostery_pos & Pos_class!="binding_interface", allosteric := T]
    dg_dt[protein==i & Pos_ref %in% allostery_pos & Pos_class=="binding_interface", orthosteric := T]
    print(paste0("Allosteric residues for ", i, ": ", paste(dg_dt[protein==i & allosteric][!duplicated(Pos_ref),Pos_ref], collapse = ",")))
    print(paste0("Orthosteric residues for ", i, ": ", paste(dg_dt[protein==i & orthosteric][!duplicated(Pos_ref),Pos_ref], collapse = ",")))
  }

  ###########################
  ### Binding ROC plot
  ###########################

  for(i in dg_dt[,unique(protein)]){
    doubledeepms__plot_binding_site_ROC(
      input_dt = copy(dg_dt)[protein==i & id!="-0-"],
      outpath = file.path(outpath, paste0(i, "_binding_site_ROC_noabs.pdf")),
      colour_scheme = colour_scheme,
      metric_names <- c(
        "b_ddg_posmean", 
        "b_ddg_posmean_conf",
        "b_ddg_wposmean", 
        "b_ddg_wposmean_conf"),
      metric_names_plot <- c(
        "Mean Binding ddG",
        "Mean Binding ddG (conf.)",
        "Weighted mean Binding ddG",
        "Weighted mean Binding ddG (conf.)"))
    doubledeepms__plot_binding_site_ROC(
      input_dt = copy(dg_dt)[protein==i & id!="-0-"],
      outpath = file.path(outpath, paste0(i, "_binding_site_ROC.pdf")),
      colour_scheme = colour_scheme,
      metric_names <- c(
        "b_ddg_posmeanabs", 
        "b_ddg_posmeanabs_conf",
        "b_ddg_wposmeanabs", 
        "b_ddg_wposmeanabs_conf"),
      metric_names_plot <- c(
        "Mean |Binding ddG|",
        "Mean |Binding ddG| (conf.)",
        "Weighted mean |Binding ddG|",
        "Weighted mean |Binding ddG| (conf.)"))
  }

  ###########################
  ### Allosteric mutations
  ###########################

  dg_list <- list()
  for(i in dg_dt[,unique(protein)]){
    yaxis_limits <- c(-3,3)
    if(i=="PSD95-PDZ3"){yaxis_limits <- c(-2,2)}
    dg_list[[i]] <- doubledeepms__allosteric_mutations_scatterplot(
      input_dt = copy(dg_dt)[protein==i],
      outpath = file.path(outpath, paste0(i, "_allosteric_mutations_scatter.pdf")),
      colour_scheme = colour_scheme,
      yaxis_limits = yaxis_limits)
    print(paste0("Surface allosteric sites for ", i, ": ", paste(dg_list[[i]][Pos_class=="surface" & allosteric][!duplicated(Pos_ref),Pos_ref], collapse = ",")))
    print(paste0("Surface sites with allosteric mutations (not in allosteric sites) for ", i, ": ", paste(dg_list[[i]][allosteric_mutation & Pos_class=="surface" & is.na(allosteric)][!duplicated(Pos_ref),Pos_ref], collapse = ",")))
  }
  dg_dt <- rbindlist(dg_list)

  ###########################
  ### Numbers of allosteric mutations
  ###########################

  print("Numbers of allosteric mutations: ")
  print(dg_dt[allosteric_mutation==T,.(total_mutations = .N),.(protein, Pos_class)][order(protein, Pos_class)])

  print("Numbers of unique residues with allosteric mutations: ")
  print(dg_dt[allosteric_mutation==T,.(total_residues = length(unique(Pos_ref))),.(protein, Pos_class)])

  plot_dt <- dg_dt[!is.na(allosteric_mutation),.(prop_mut = sum(allosteric_mutation)/length(allosteric_mutation)*100, scHAmin_ligand = unique(scHAmin_ligand), Pos_class=unique(Pos_class), allosteric=unique(allosteric)),.(Pos_ref, protein)]
  cor_dt <- plot_dt[,.(cor = round(cor(scHAmin_ligand, prop_mut, use = "pairwise.complete", method = "spearman"), 2), scHAmin_ligand = 12, prop_mut = 100),protein]
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(scHAmin_ligand, prop_mut)) +
    ggplot2::geom_smooth(method = "lm", se = F, color = "grey", linetype = 2, formula = 'y ~ x') + 
    ggplot2::geom_point(ggplot2::aes(shape = !is.na(allosteric), color = Pos_class), size = 2) +
    ggplot2::xlab(expression("Distance to ligand ("*ring(A)*")")) +
    ggplot2::ylab("%Allosteric mutations per residue") +
    ggplot2::facet_wrap(protein~., scales = "free", ncol = 1) +
    ggplot2::labs(color = "Residue\nposition", shape = "Major allosteric\nsite") +
    ggplot2::geom_text(data = cor_dt, ggplot2::aes(label=paste("Spearman's rho = ", cor, sep=""))) +
    ggplot2::theme_bw()
  if(!is.null(colour_scheme)){
    d <- d + ggplot2::scale_colour_manual(values = unlist(colour_scheme[["shade 0"]][c(3, 4)]))
  }
  ggplot2::ggsave(file.path(outpath, "allosteric_mutations_vs_distance.pdf"), d, width = 4, height = 6, useDingbats=FALSE)

  ###########################
  ### Proportion of mutations in allosteric sites that are allosteric mutations
  ###########################

  for(i in dg_dt[,unique(protein)]){
    temp_prop <- dg_dt[protein==i & allosteric & !is.na(allosteric_mutation),sum(allosteric_mutation)/.N]
    print(paste0("Proportion of allosteric mutations at allosteric sites (", i, "): ", format(temp_prop, digits=2, scientific=T)))
  }

  ###########################
  ### Correlation of mean absolute binding free energy changes with Ohm results
  ###########################

  for(i in dg_dt[,unique(protein)]){
    ohm_dt <- fread(ohm_file_list[[i]])
    names(ohm_dt) <- c("Pos_ref", "ACI")
    ohm_dt[, protein := i]
    plot_dt <- merge(dg_dt[,.(Pos_ref, b_ddg_wposmeanabs, Pos_class, protein, allosteric)], ohm_dt, by = c("Pos_ref", "protein"), all.x = T)
    plot_dt <- plot_dt[protein==i & Pos_class!="binding_interface"][!duplicated(Pos_ref)]
    cor_dt <- plot_dt[,.(cor = round(cor(ACI, b_ddg_wposmeanabs, use = "pairwise.complete", method = "spearman"), 2), ACI = 0.6, b_ddg_wposmeanabs = 0.9),protein]
    #Plot
    d <- ggplot2::ggplot(plot_dt,ggplot2::aes(ACI, b_ddg_wposmeanabs)) +
      ggplot2::geom_smooth(method = "lm", se = F, color = "black", linetype = 2, formula = 'y ~ x') + 
      ggplot2::geom_point(ggplot2::aes(color = !is.na(allosteric)), size = 2) +
      ggplot2::xlab("Allosteric coupling intensity") +
      ggplot2::ylab(expression("Weighted mean |Binding "*Delta*Delta*"G|")) +
      ggplot2::labs(color = "Major\nallosteric\nsite") +
      ggplot2::geom_text(data = cor_dt, ggplot2::aes(label=paste("Spearman's rho = ", cor, sep=""))) +
      ggplot2::theme_bw()
    if(!is.null(colour_scheme)){
      d <- d + ggplot2::scale_colour_manual(values = c("darkgrey", unlist(colour_scheme[["shade 0"]][[2]])))
    }
    ggplot2::ggsave(file.path(outpath, paste0("ohm_scatter_", i, ".pdf")), d, width = 4, height = 3, useDingbats=FALSE)
    #Remove binding interface adjacent residues
    plot_dt <- merge(dg_dt[,.(Pos_ref, b_ddg_wposmeanabs, Pos_class, protein, allosteric)], ohm_dt, by = c("Pos_ref", "protein"), all.x = T)
    bi_residues <- plot_dt[protein==i & Pos_class=="binding_interface",unique(Pos_ref)]
    plot_dt <- plot_dt[protein==i & !(Pos_ref-1) %in% bi_residues & !(Pos_ref+1) %in% bi_residues & Pos_class!="binding_interface"][!duplicated(Pos_ref)]
    cor_dt <- plot_dt[,.(cor = round(cor(ACI, b_ddg_wposmeanabs, use = "pairwise.complete", method = "spearman"), 2), ACI = 0.6, b_ddg_wposmeanabs = 0.9),protein]
    #Plot
    d <- ggplot2::ggplot(plot_dt,ggplot2::aes(ACI, b_ddg_wposmeanabs)) +
      ggplot2::geom_smooth(method = "lm", se = F, color = "black", linetype = 2, formula = 'y ~ x') +
      ggplot2::geom_point(ggplot2::aes(color = !is.na(allosteric)), size = 2) +
      ggplot2::xlab("Allosteric coupling intensity") +
      ggplot2::ylab(expression("Weighted mean |Binding "*Delta*Delta*"G|")) +
      ggplot2::labs(color = "Major\nallosteric\nsite") +
      ggplot2::geom_text(data = cor_dt, ggplot2::aes(label=paste("Spearman's rho = ", cor, sep=""))) +
      ggplot2::theme_bw()
    if(!is.null(colour_scheme)){
      d <- d + ggplot2::scale_colour_manual(values = c("darkgrey", unlist(colour_scheme[["shade 0"]][[2]])))
    }
    ggplot2::ggsave(file.path(outpath, paste0("ohm_scatter_", i, "_noadjacent.pdf")), d, width = 4, height = 3, useDingbats=FALSE)
    #Set binding interface adjacent residue scores to 1
    plot_dt <- merge(dg_dt[,.(Pos_ref, b_ddg_wposmeanabs, Pos_class, protein, allosteric)], ohm_dt, by = c("Pos_ref", "protein"), all.x = T)
    bi_residues <- plot_dt[protein==i & Pos_class=="binding_interface",unique(Pos_ref)]
    plot_dt <- plot_dt[protein==i & Pos_class!="binding_interface"][!duplicated(Pos_ref)]
    plot_dt[(Pos_ref-1) %in% bi_residues | (Pos_ref+1) %in% bi_residues, ACI := 1]
    cor_dt <- plot_dt[,.(cor = round(cor(ACI, b_ddg_wposmeanabs, use = "pairwise.complete", method = "spearman"), 2), ACI = 0.6, b_ddg_wposmeanabs = 0.9),protein]
    #Plot
    d <- ggplot2::ggplot(plot_dt,ggplot2::aes(ACI, b_ddg_wposmeanabs)) +
      ggplot2::geom_smooth(method = "lm", se = F, color = "black", linetype = 2, formula = 'y ~ x') + 
      ggplot2::geom_point(ggplot2::aes(color = !is.na(allosteric)), size = 2) +
      ggplot2::xlab("Allosteric coupling intensity") +
      ggplot2::ylab(expression("Weighted mean |Binding "*Delta*Delta*"G|")) +
      ggplot2::labs(color = "Major\nallosteric\nsite") +
      ggplot2::geom_text(data = cor_dt, ggplot2::aes(label=paste("Spearman's rho = ", cor, sep=""))) +
      ggplot2::theme_bw()
    if(!is.null(colour_scheme)){
      d <- d + ggplot2::scale_colour_manual(values = c("darkgrey", unlist(colour_scheme[["shade 0"]][[2]])))
    }
    ggplot2::ggsave(file.path(outpath, paste0("ohm_scatter_", i, "_adjacent1.pdf")), d, width = 4, height = 3, useDingbats=FALSE)
  }

  ###########################
  ### Correlation of proportion of mutations in allosteric sites with Ohm results
  ###########################

  for(i in dg_dt[,unique(protein)]){
    ohm_dt <- fread(ohm_file_list[[i]])
    names(ohm_dt) <- c("Pos_ref", "ACI")
    ohm_dt[, protein := i]
    plot_dt <- dg_dt[,.(prop_mut = sum(allosteric_mutation, na.rm = T)/sum(!is.na(allosteric_mutation))*100, Pos_class=unique(Pos_class), allosteric=sum(allosteric, na.rm = T)!=0),.(Pos_ref, protein)]
    plot_dt <- merge(plot_dt[,.(Pos_ref, Pos_class, protein, prop_mut, allosteric)], ohm_dt, by = c("Pos_ref", "protein"), all.x = T)
    plot_dt <- plot_dt[protein==i & Pos_class!="binding_interface"][!duplicated(Pos_ref)]
    cor_dt <- plot_dt[,.(cor = round(cor(ACI, prop_mut, use = "pairwise.complete", method = "spearman"), 2), ACI = 0.6, prop_mut = 90),protein]
    #Plot
    d <- ggplot2::ggplot(plot_dt,ggplot2::aes(ACI, prop_mut)) +
      ggplot2::geom_smooth(method = "lm", se = F, color = "black", linetype = 2, formula = 'y ~ x') + 
      ggplot2::geom_point(data = plot_dt[allosteric==F], ggplot2::aes(color = allosteric), size = 2) +
      ggplot2::geom_point(data = plot_dt[allosteric==T],ggplot2::aes(color = allosteric), size = 2) +
      ggplot2::xlab("Allosteric coupling intensity") +
      ggplot2::ylab("%Allosteric mutations per residue") +
      ggplot2::labs(color = "Major\nallosteric\nsite") +
      ggplot2::geom_text(data = cor_dt, ggplot2::aes(label=paste("Spearman's rho = ", cor, sep=""))) +
      ggplot2::theme_bw()
    if(!is.null(colour_scheme)){
      d <- d + ggplot2::scale_colour_manual(values = c("darkgrey", unlist(colour_scheme[["shade 0"]][[2]])))
    }
    ggplot2::ggsave(file.path(outpath, paste0("ohm_scatter_", i, "_allosteric_mutations.pdf")), d, width = 4, height = 3, useDingbats=FALSE)
    #Remove binding interface adjacent residues
    plot_dt <- dg_dt[,.(prop_mut = sum(allosteric_mutation, na.rm = T)/sum(!is.na(allosteric_mutation))*100, Pos_class=unique(Pos_class), allosteric=sum(allosteric, na.rm = T)!=0),.(Pos_ref, protein)]
    plot_dt <- merge(plot_dt[,.(Pos_ref, Pos_class, protein, prop_mut, allosteric)], ohm_dt, by = c("Pos_ref", "protein"), all.x = T)
    bi_residues <- plot_dt[protein==i & Pos_class=="binding_interface",unique(Pos_ref)]
    plot_dt <- plot_dt[protein==i & !(Pos_ref-1) %in% bi_residues & !(Pos_ref+1) %in% bi_residues & Pos_class!="binding_interface"][!duplicated(Pos_ref)]
    cor_dt <- plot_dt[,.(cor = round(cor(ACI, prop_mut, use = "pairwise.complete", method = "spearman"), 2), ACI = 0.6, prop_mut = 90),protein]
    #Plot
    d <- ggplot2::ggplot(plot_dt,ggplot2::aes(ACI, prop_mut)) +
      ggplot2::geom_smooth(method = "lm", se = F, color = "black", linetype = 2, formula = 'y ~ x') + 
      ggplot2::geom_point(data = plot_dt[allosteric==F], ggplot2::aes(color = allosteric), size = 2) +
      ggplot2::geom_point(data = plot_dt[allosteric==T],ggplot2::aes(color = allosteric), size = 2) +
      ggplot2::xlab("Allosteric coupling intensity") +
      ggplot2::ylab("%Allosteric mutations per residue") +
      ggplot2::labs(color = "Major\nallosteric\nsite") +
      ggplot2::geom_text(data = cor_dt, ggplot2::aes(label=paste("Spearman's rho = ", cor, sep=""))) +
      ggplot2::theme_bw()
    if(!is.null(colour_scheme)){
      d <- d + ggplot2::scale_colour_manual(values = c("darkgrey", unlist(colour_scheme[["shade 0"]][[2]])))
    }
    ggplot2::ggsave(file.path(outpath, paste0("ohm_scatter_", i, "_allosteric_mutations_noadjacent.pdf")), d, width = 4, height = 3, useDingbats=FALSE)
    #Set binding interface adjacent residue scores to 1
    plot_dt <- dg_dt[,.(prop_mut = sum(allosteric_mutation, na.rm = T)/sum(!is.na(allosteric_mutation))*100, Pos_class=unique(Pos_class), allosteric=sum(allosteric, na.rm = T)!=0),.(Pos_ref, protein)]
    plot_dt <- merge(plot_dt[,.(Pos_ref, Pos_class, protein, prop_mut, allosteric)], ohm_dt, by = c("Pos_ref", "protein"), all.x = T)
    bi_residues <- plot_dt[protein==i & Pos_class=="binding_interface",unique(Pos_ref)]
    plot_dt <- plot_dt[protein==i & Pos_class!="binding_interface"][!duplicated(Pos_ref)]
    plot_dt[(Pos_ref-1) %in% bi_residues | (Pos_ref+1) %in% bi_residues, ACI := 1]
    cor_dt <- plot_dt[,.(cor = round(cor(ACI, prop_mut, use = "pairwise.complete", method = "spearman"), 2), ACI = 0.6, prop_mut = 90),protein]
    #Plot
    d <- ggplot2::ggplot(plot_dt,ggplot2::aes(ACI, prop_mut)) +
      ggplot2::geom_smooth(method = "lm", se = F, color = "black", linetype = 2, formula = 'y ~ x') + 
      ggplot2::geom_point(data = plot_dt[allosteric==F], ggplot2::aes(color = allosteric), size = 2) +
      ggplot2::geom_point(data = plot_dt[allosteric==T],ggplot2::aes(color = allosteric), size = 2) +
      ggplot2::xlab("Allosteric coupling intensity") +
      ggplot2::ylab("%Allosteric mutations per residue") +
      ggplot2::labs(color = "Major\nallosteric\nsite") +
      ggplot2::geom_text(data = cor_dt, ggplot2::aes(label=paste("Spearman's rho = ", cor, sep=""))) +
      ggplot2::theme_bw()
    if(!is.null(colour_scheme)){
      d <- d + ggplot2::scale_colour_manual(values = c("darkgrey", unlist(colour_scheme[["shade 0"]][[2]])))
    }
    ggplot2::ggsave(file.path(outpath, paste0("ohm_scatter_", i, "_allosteric_mutations_adjacent1.pdf")), d, width = 4, height = 3, useDingbats=FALSE)
  }

  ###########################
  ### Predicting literature sites using mean absolute binding free energy changes
  ###########################

  #Sector residues
  for(i in dg_dt[,unique(protein)]){
    if(length(literature_list[[i]][["sector"]]!=0)){
      in_sector_means <- dg_dt[protein==i & !is.na(b_ddg_wposmeanabs) & id!="-0-"][!duplicated(Pos_ref)][Pos_ref %in% literature_list[[i]][["sector"]],b_ddg_wposmeanabs]
      out_sector_means <- dg_dt[protein==i & !is.na(b_ddg_wposmeanabs) & id!="-0-"][!duplicated(Pos_ref)][!Pos_ref %in% literature_list[[i]][["sector"]],b_ddg_wposmeanabs]
      temp_test <- doubledeepms__mann_whitney_U_wrapper(in_sector_means, out_sector_means)
      print(paste0("Predicting sector residues using position-wise weighted mean of binding free energies changes (", i, "): p-value=", 
        format(temp_test[["p_value"]], digits=2, scientific=T), " AUC=", 
        round(temp_test[["effect_size"]], 2), " n=",
        length(c(in_sector_means, out_sector_means))))
   }
  }

  #Class switching residues
  for(i in dg_dt[,unique(protein)]){
    if(length(literature_list[[i]][["class_switching"]]!=0)){
      in_sector_means <- dg_dt[protein==i & !is.na(b_ddg_wposmeanabs) & id!="-0-"][!duplicated(Pos_ref)][Pos_ref %in% literature_list[[i]][["class_switching"]],b_ddg_wposmeanabs]
      out_sector_means <- dg_dt[protein==i & !is.na(b_ddg_wposmeanabs) & id!="-0-"][!duplicated(Pos_ref)][!Pos_ref %in% literature_list[[i]][["class_switching"]],b_ddg_wposmeanabs]
      temp_test <- doubledeepms__mann_whitney_U_wrapper(in_sector_means, out_sector_means)
      print(paste0("Predicting class-switching residues using position-wise weighted mean of binding free energies changes (", i, "): p-value=", 
        format(temp_test[["p_value"]], digits=2, scientific=T), " AUC=", 
        round(temp_test[["effect_size"]], 2), " n=",
        length(c(in_sector_means, out_sector_means))))
      #Predicting class-switching residues not classified as allosteric
      in_sector_means <- dg_dt[protein==i & !is.na(b_ddg_pred) & id!="-0-"][Pos_ref %in% literature_list[[i]][["class_switching"]] & is.na(allosteric) & is.na(orthosteric),abs(b_ddg_pred)]
      out_sector_means <- dg_dt[protein==i & !is.na(b_ddg_pred) & id!="-0-"][!Pos_ref %in% literature_list[[i]][["class_switching"]] & is.na(allosteric) & is.na(orthosteric),abs(b_ddg_pred)]
      temp_test <- doubledeepms__mann_whitney_U_wrapper(in_sector_means, out_sector_means)
      print(paste0("Binding free energies of mutations at non-allosteric/orthosteric residues that are class-switching versus remainder (", i, "): p-value=", 
        format(temp_test[["p_value"]], digits=2, scientific=T), " AUC=", 
        round(temp_test[["effect_size"]], 2), " n=",
        length(c(in_sector_means, out_sector_means))))
    }
  }

  ###########################
  ### Enrichment of allosteric mutations in literature sites
  ###########################

  #Literature sites
  result_list <- list()
  for(i in dg_dt[,unique(protein)]){
    if(i %in% names(annotation_list)){
      anno_dt <- fread(annotation_list[[i]])
      for(lset_name in names(anno_dt[,-1])){
        lset <- anno_dt[get(lset_name)==1,Pos_ref]
        in_lset_allo <- dg_dt[protein==i & Pos_class!="binding_interface"][Pos_ref %in% lset & allosteric_mutation==T,.N]
        out_lset_allo <- dg_dt[protein==i & Pos_class!="binding_interface"][!Pos_ref %in% lset & allosteric_mutation==T,.N]
        in_lset_nallo <- dg_dt[protein==i & Pos_class!="binding_interface"][Pos_ref %in% lset & allosteric_mutation==F,.N]
        out_lset_nallo <- dg_dt[protein==i & Pos_class!="binding_interface"][!Pos_ref %in% lset & allosteric_mutation==F,.N]
        temp_test <- fisher.test(matrix(c(in_lset_allo, out_lset_allo, in_lset_nallo, out_lset_nallo), nrow = 2))
        print(paste0("Enrichment of allosteric mutations in ", lset_name, " (", i, "): p-value=", 
          format(temp_test$p.value, digits=2, scientific=T), " odds ratio=", 
          round(temp_test$estimate, 2), " n=",
          sum(c(in_lset_allo, out_lset_allo, in_lset_nallo, out_lset_nallo))))
        result_list <- c(result_list, list(data.table(protein = i, set_name = paste0(lset_name, " (n = ", length(lset), ")"), mutation = "in", odds_ratio = temp_test$estimate, p_value = temp_test$p.value)))
      }
    }
  }

  #Plot
  plot_dt <- rbindlist(result_list)
  plot_dt[, set_name := factor(set_name, levels = plot_dt[order(odds_ratio, decreasing = T),set_name])]
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(set_name, log2(odds_ratio), alpha = p_value<2.2e-16)) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::facet_wrap(protein~., scales = "free", ncol = 1) +
    ggplot2::xlab("Dataset") +
    ggplot2::ylab("log2(odds ratio)") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
  if(!is.null(colour_scheme)){
    d <- d + ggplot2::scale_colour_manual(values = unlist(colour_scheme[["shade 0"]][c(1, 2)]))
    d <- d + ggplot2::scale_alpha_manual(values = c(0.5, 1))
  }
  ggplot2::ggsave(file.path(outpath, "allosteric_mutations_literature_associations.pdf"), d, width = 5, height = 4, useDingbats=FALSE)

  ###########################
  ### Enrichment of allosteric mutations in certain mutant residues
  ###########################

  #Mutant residues
  result_list <- list()
  bset_list <- list(
    "Charged" = c("R", "H", "D", "E", "K"), 
    "Hydrophobic" = c("A", "V", "I", "L", "M", "F", "Y", "W"))
  aa_list <- as.list(unlist(strsplit("GAVLMIFYWKRHDESTCNQP", "")))
  names(aa_list) <- unlist(aa_list)
  bset_list <- c(bset_list, aa_list)
  for(lset_name in names(bset_list)){
    for(i in dg_dt[,unique(protein)]){
      lset <- bset_list[[lset_name]]
      in_lset_allo <- dg_dt[protein==i & Pos_class!="binding_interface"][WT_AA %in% lset & allosteric_mutation==T,.N]
      out_lset_allo <- dg_dt[protein==i & Pos_class!="binding_interface"][!WT_AA %in% lset & allosteric_mutation==T,.N]
      in_lset_nallo <- dg_dt[protein==i & Pos_class!="binding_interface"][WT_AA %in% lset & allosteric_mutation==F,.N]
      out_lset_nallo <- dg_dt[protein==i & Pos_class!="binding_interface"][!WT_AA %in% lset & allosteric_mutation==F,.N]
      temp_test <- fisher.test(matrix(c(in_lset_allo, out_lset_allo, in_lset_nallo, out_lset_nallo), nrow = 2))
      if(lset_name %in% c("G", "P")){
        print(paste0("Enrichment of allosteric mutations from ", lset_name, " (", i, "): p-value=", 
          format(temp_test$p.value, digits=2, scientific=T), " odds ratio=", 
          round(temp_test$estimate, 2), " n=",
          sum(c(in_lset_allo, out_lset_allo, in_lset_nallo, out_lset_nallo))))
      }
      result_list <- c(result_list, list(data.table(protein = i, set_name = lset_name, mutation = "WT", odds_ratio = temp_test$estimate, p_value = temp_test$p.value)))
      in_lset_allo <- dg_dt[protein==i & Pos_class!="binding_interface"][Mut %in% lset & allosteric_mutation==T,.N]
      out_lset_allo <- dg_dt[protein==i & Pos_class!="binding_interface"][!Mut %in% lset & allosteric_mutation==T,.N]
      in_lset_nallo <- dg_dt[protein==i & Pos_class!="binding_interface"][Mut %in% lset & allosteric_mutation==F,.N]
      out_lset_nallo <- dg_dt[protein==i & Pos_class!="binding_interface"][!Mut %in% lset & allosteric_mutation==F,.N]
      temp_test <- fisher.test(matrix(c(in_lset_allo, out_lset_allo, in_lset_nallo, out_lset_nallo), nrow = 2))
      if(lset_name %in% c("G", "P")){
        print(paste0("Enrichment of allosteric mutations to ", lset_name, " (", i, "): p-value=", 
          format(temp_test$p.value, digits=2, scientific=T), " odds ratio=", 
          round(temp_test$estimate, 2), " n=",
          sum(c(in_lset_allo, out_lset_allo, in_lset_nallo, out_lset_nallo))))
      }
      result_list <- c(result_list, list(data.table(protein = i, set_name = lset_name, mutation = "Mutant", odds_ratio = temp_test$estimate, p_value = temp_test$p.value)))
    }
  }

  #Plot
  plot_dt <- rbindlist(result_list)
  plot_dt[, set_name := factor(set_name, levels = plot_dt[,.(moddsratio = mean(odds_ratio)),set_name][order(moddsratio, decreasing = T),set_name])]
  plot_dt[, mutation := factor(mutation, levels = c("WT", "Mutant"))]
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(set_name, log2(odds_ratio), fill = mutation, alpha = p_value<0.05)) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::facet_wrap(protein~., scales = "free", ncol = 1) +
    ggplot2::xlab("Residue type") +
    ggplot2::ylab("log2(odds ratio)") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
  if(!is.null(colour_scheme)){
    d <- d + ggplot2::scale_fill_manual(values = unlist(colour_scheme[["shade 0"]][c(1, 3)]))
    d <- d + ggplot2::scale_alpha_manual(values = c(0.5, 1))
  }
  ggplot2::ggsave(file.path(outpath, "allosteric_mutations_mutated_residues.pdf"), d, width = 6, height = 7, useDingbats=FALSE)
  #Save
  plot_dt_allres <- plot_dt

  ###########################
  ### Enrichment of allosteric mutations in certain mutant residues - no loops
  ###########################

  #Mutant residues
  result_list <- list()
  bset_list <- list(
    "Charged" = c("R", "H", "D", "E", "K"), 
    "Hydrophobic" = c("A", "V", "I", "L", "M", "F", "Y", "W"))
  aa_list <- as.list(unlist(strsplit("GAVLMIFYWKRHDESTCNQP", "")))
  names(aa_list) <- unlist(aa_list)
  bset_list <- c(bset_list, aa_list)
  for(lset_name in names(bset_list)){
    for(i in dg_dt[,unique(protein)]){
      lset <- bset_list[[lset_name]]
      in_lset_allo <- dg_dt[protein==i & Pos_class!="binding_interface" & !is.na(SS)][WT_AA %in% lset & allosteric_mutation==T,.N]
      out_lset_allo <- dg_dt[protein==i & Pos_class!="binding_interface" & !is.na(SS)][!WT_AA %in% lset & allosteric_mutation==T,.N]
      in_lset_nallo <- dg_dt[protein==i & Pos_class!="binding_interface" & !is.na(SS)][WT_AA %in% lset & allosteric_mutation==F,.N]
      out_lset_nallo <- dg_dt[protein==i & Pos_class!="binding_interface" & !is.na(SS)][!WT_AA %in% lset & allosteric_mutation==F,.N]
      temp_matrix <- matrix(c(in_lset_allo, out_lset_allo, in_lset_nallo, out_lset_nallo), nrow = 2)
      rownames(temp_matrix) <- c("In set", "Out set")
      colnames(temp_matrix) <- c("Allosteric", "Not allosteric")
      temp_test <- fisher.test(temp_matrix)
      if(lset_name %in% c("G", "P")){
        print(paste0("Enrichment of allosteric mutations (not in loops) from ", lset_name, " (", i, "): p-value=", 
          format(temp_test$p.value, digits=2, scientific=T), " odds ratio=", 
          round(temp_test$estimate, 2), " n=",
          sum(c(in_lset_allo, out_lset_allo, in_lset_nallo, out_lset_nallo))))
        print(temp_matrix)
      }
      result_list <- c(result_list, list(data.table(protein = i, set_name = lset_name, mutation = "WT", odds_ratio = temp_test$estimate, p_value = temp_test$p.value)))
      in_lset_allo <- dg_dt[protein==i & Pos_class!="binding_interface" & !is.na(SS)][Mut %in% lset & allosteric_mutation==T,.N]
      out_lset_allo <- dg_dt[protein==i & Pos_class!="binding_interface" & !is.na(SS)][!Mut %in% lset & allosteric_mutation==T,.N]
      in_lset_nallo <- dg_dt[protein==i & Pos_class!="binding_interface" & !is.na(SS)][Mut %in% lset & allosteric_mutation==F,.N]
      out_lset_nallo <- dg_dt[protein==i & Pos_class!="binding_interface" & !is.na(SS)][!Mut %in% lset & allosteric_mutation==F,.N]
      temp_test <- fisher.test(matrix(c(in_lset_allo, out_lset_allo, in_lset_nallo, out_lset_nallo), nrow = 2))
      if(lset_name %in% c("G", "P")){
        print(paste0("Enrichment of allosteric mutations (not in loops) to ", lset_name, " (", i, "): p-value=", 
          format(temp_test$p.value, digits=2, scientific=T), " odds ratio=", 
          round(temp_test$estimate, 2), " n=",
          sum(c(in_lset_allo, out_lset_allo, in_lset_nallo, out_lset_nallo))))
      }
      result_list <- c(result_list, list(data.table(protein = i, set_name = lset_name, mutation = "Mutant", odds_ratio = temp_test$estimate, p_value = temp_test$p.value)))
    }
  }

  #Plot
  plot_dt <- rbindlist(result_list)
  plot_dt[, set_name := factor(set_name, levels = plot_dt[,.(moddsratio = mean(odds_ratio)),set_name][order(moddsratio, decreasing = T),set_name])]
  plot_dt[, mutation := factor(mutation, levels = c("WT", "Mutant"))]
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(set_name, log2(odds_ratio), fill = mutation, alpha = p_value<0.05)) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::facet_wrap(protein~., scales = "free", ncol = 1) +
    ggplot2::xlab("Residue type") +
    ggplot2::ylab("log2(odds ratio)") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
  if(!is.null(colour_scheme)){
    d <- d + ggplot2::scale_fill_manual(values = unlist(colour_scheme[["shade 0"]][c(1, 3)]))
    d <- d + ggplot2::scale_alpha_manual(values = c(0.5, 1))
  }
  ggplot2::ggsave(file.path(outpath, "allosteric_mutations_mutated_residues_noloops.pdf"), d, width = 6, height = 7, useDingbats=FALSE)

  #Plot restricting to G and P
  plot_dt <- plot_dt[set_name %in% c("G", "P")]
  plot_dt[, subset := "noloop"]
  plot_dt_allres <- plot_dt_allres[set_name %in% c("G", "P")]
  plot_dt_allres[, subset := "all"]
  plot_dt_combined <- rbind(plot_dt, plot_dt_allres)
  d <- ggplot2::ggplot(plot_dt_combined,ggplot2::aes(mutation, log2(odds_ratio), fill = subset)) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::facet_grid(protein~set_name, scales = "free") +
    ggplot2::xlab("Residue type") +
    ggplot2::ylab("log2(odds ratio)") +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
  if(!is.null(colour_scheme)){
    d <- d + ggplot2::scale_fill_manual(values = c("grey", unlist(colour_scheme[["shade 0"]][c(1)])))
  }
  ggplot2::ggsave(file.path(outpath, "allosteric_mutations_mutated_residues_noloops_GP.pdf"), d, width = 3, height = 4, useDingbats=FALSE)

  # Enrichment of major allosteric sites in Glycines
  g_allo <- dg_dt[!duplicated(paste0(protein, ":", Pos_ref))][Pos_class != "binding_interface"][allosteric==T & WT_AA=="G",.N]
  ng_allo <- dg_dt[!duplicated(paste0(protein, ":", Pos_ref))][Pos_class != "binding_interface"][allosteric==T & WT_AA!="G",.N]
  g_nallo <- dg_dt[!duplicated(paste0(protein, ":", Pos_ref))][Pos_class != "binding_interface"][is.na(allosteric) & WT_AA=="G",.N]
  ng_nallo <- dg_dt[!duplicated(paste0(protein, ":", Pos_ref))][Pos_class != "binding_interface"][is.na(allosteric) & WT_AA!="G",.N]
  temp_matrix <- matrix(c(g_allo, ng_allo, g_nallo, ng_nallo), nrow = 2)
  rownames(temp_matrix) <- c("In set", "Out set")
  colnames(temp_matrix) <- c("Allosteric", "Not allosteric")
  temp_test <- fisher.test(temp_matrix)
  print(paste0("Amongst residues not in binding interface, enrichment of Glycine residues for major allosteric sites: p-value=", 
    format(temp_test$p.value, digits=2, scientific=T), " odds ratio=", 
    round(temp_test$estimate, 2), " n=",
    sum(c(g_allo, ng_allo, g_nallo, ng_nallo))))
  print(temp_matrix)

  # Amongst Glycines, enrichment of major allosteric sites in non-loop residues
  nloop_allo <- dg_dt[!duplicated(paste0(protein, ":", Pos_ref))][WT_AA=="G" & Pos_class != "binding_interface"][allosteric==T & !is.na(SS),.N]
  loop_allo <- dg_dt[!duplicated(paste0(protein, ":", Pos_ref))][WT_AA=="G" & Pos_class != "binding_interface"][allosteric==T & is.na(SS),.N]
  nloop_nallo <- dg_dt[!duplicated(paste0(protein, ":", Pos_ref))][WT_AA=="G" & Pos_class != "binding_interface"][is.na(allosteric) & !is.na(SS),.N]
  loop_nallo <- dg_dt[!duplicated(paste0(protein, ":", Pos_ref))][WT_AA=="G" & Pos_class != "binding_interface"][is.na(allosteric) & is.na(SS),.N]
  temp_matrix <- matrix(c(nloop_allo, loop_allo, nloop_nallo, loop_nallo), nrow = 2)
  rownames(temp_matrix) <- c("In set", "Out set")
  colnames(temp_matrix) <- c("Allosteric", "Not allosteric")
  temp_test <- fisher.test(temp_matrix)
  print(paste0("Amongst Glycines not in binding interface, enrichment of non-loop residues for major allosteric sites: p-value=", 
    format(temp_test$p.value, digits=2, scientific=T), " odds ratio=", 
    round(temp_test$estimate, 2), " n=",
    sum(c(nloop_allo, loop_allo, nloop_nallo, loop_nallo))))
  print(temp_matrix)

  # Enrichment of major allosteric sites in non-loop Glycines
  nloop_allo <- dg_dt[!duplicated(paste0(protein, ":", Pos_ref))][Pos_class != "binding_interface"][allosteric==T & (!is.na(SS) & WT_AA=="G"),.N]
  loop_allo <- dg_dt[!duplicated(paste0(protein, ":", Pos_ref))][Pos_class != "binding_interface"][allosteric==T & !(!is.na(SS) & WT_AA=="G"),.N]
  nloop_nallo <- dg_dt[!duplicated(paste0(protein, ":", Pos_ref))][Pos_class != "binding_interface"][is.na(allosteric) & (!is.na(SS) & WT_AA=="G"),.N]
  loop_nallo <- dg_dt[!duplicated(paste0(protein, ":", Pos_ref))][Pos_class != "binding_interface"][is.na(allosteric) & !(!is.na(SS) & WT_AA=="G"),.N]
  temp_matrix <- matrix(c(nloop_allo, loop_allo, nloop_nallo, loop_nallo), nrow = 2)
  rownames(temp_matrix) <- c("In set", "Out set")
  colnames(temp_matrix) <- c("Allosteric", "Not allosteric")
  temp_test <- fisher.test(temp_matrix)
  print(paste0("Amongst residues not in binding interface, enrichment of non-loop Glycine residues for major allosteric sites: p-value=", 
    format(temp_test$p.value, digits=2, scientific=T), " odds ratio=", 
    round(temp_test$estimate, 2), " n=",
    sum(c(nloop_allo, loop_allo, nloop_nallo, loop_nallo))))
  print(temp_matrix)

  ### Enrichment of allosteric mutations in negatively charged surface sites
  ###########################

  for(i in dg_dt[,unique(protein)]){
    temp_in <- dg_dt[protein==i & allosteric_mutation==T & Pos_class=="surface",WT_AA %in% c("D", "E")]
    temp_out <- dg_dt[protein==i & allosteric_mutation==F & Pos_class=="surface",WT_AA %in% c("D", "E")]
    temp_test <- fisher.test(matrix(c(sum(temp_in), sum(!temp_in), sum(temp_out), sum(!temp_out)), nrow = 2))
    print(paste0("Enrichment of allosteric mutations in negatively charged surface sites (", i, "): p-value=", format(temp_test$p.value, digits=2, scientific=T), " odds ratio=", round(temp_test$estimate, 2)))
  }

  ### Enrichment of allosteric mutations to negatively charged surface sites
  ###########################

  for(i in dg_dt[,unique(protein)]){
    temp_in <- dg_dt[protein==i & allosteric_mutation==T & Pos_class=="surface",Mut %in% c("D", "E")]
    temp_out <- dg_dt[protein==i & allosteric_mutation==F & Pos_class=="surface",Mut %in% c("D", "E")]
    temp_test <- fisher.test(matrix(c(sum(temp_in), sum(!temp_in), sum(temp_out), sum(!temp_out)), nrow = 2))
    print(paste0("Enrichment of allosteric mutations to negatively charged surface sites (", i, "): p-value=", format(temp_test$p.value, digits=2, scientific=T), " odds ratio=", round(temp_test$estimate, 2)))
  }

  ### Enrichment of allosteric mutations in net -1 charged mutations at surface sites
  ###########################

  #Delta charge
  dg_dt[id!="-0-", WT_AA_charge := 0]
  dg_dt[WT_AA %in% c("R", "H", "K"), WT_AA_charge := 1]
  dg_dt[WT_AA %in% c("D", "E"), WT_AA_charge := -1]
  dg_dt[id!="-0-", Mut_charge := 0]
  dg_dt[Mut %in% c("R", "H", "K"), Mut_charge := 1]
  dg_dt[Mut %in% c("D", "E"), Mut_charge := -1]
  dg_dt[id!="-0-", delta_charge := Mut_charge-WT_AA_charge]

  for(i in dg_dt[,unique(protein)]){
    temp_in <- dg_dt[protein==i & allosteric_mutation==T & Pos_class=="surface",delta_charge==(-1)]
    temp_out <- dg_dt[protein==i & allosteric_mutation==F & Pos_class=="surface",delta_charge==(-1)]
    temp_test <- fisher.test(matrix(c(sum(temp_in), sum(!temp_in), sum(temp_out), sum(!temp_out)), nrow = 2))
    print(paste0("Enrichment of allosteric mutations in net -1 charged mutations at surface sites (", i, "): p-value=", format(temp_test$p.value, digits=2, scientific=T), " odds ratio=", round(temp_test$estimate, 2)))
  }

  ###########################
  ### Where are strongest binding effects?
  ###########################

  #Thresholded data
  plot_list <- list()
  plot_list_conf <- list()
  for(i in c(1:40)/20){
    plot_list[[as.character(i)]] <- dg_dt[(b_ddg_pred-1.96*b_ddg_pred_sd)>i & id!="-0-",]
    plot_list[[as.character(i)]][, b_ddg_pred_class := i]
    plot_list[[as.character(-i)]] <- dg_dt[(b_ddg_pred+1.96*b_ddg_pred_sd)<(-i) & id!="-0-",]
    plot_list[[as.character(-i)]][, b_ddg_pred_class := (-i)]
    plot_list_conf[[as.character(i)]] <- dg_dt[b_ddg_pred_conf==T & (b_ddg_pred-1.96*b_ddg_pred_sd)>i & id!="-0-",]
    plot_list_conf[[as.character(i)]][, b_ddg_pred_class := i]
    plot_list_conf[[as.character(-i)]] <- dg_dt[b_ddg_pred_conf==T & (b_ddg_pred+1.96*b_ddg_pred_sd)<(-i) & id!="-0-",]
    plot_list_conf[[as.character(-i)]][, b_ddg_pred_class := (-i)]
  }
  plot_dt <- rbindlist(plot_list)
  plot_dt_conf <- rbindlist(plot_list_conf)

  #Plot - all
  plot_dt_all <- plot_dt[,.(num_mutations = .N),.(b_ddg_pred_class, Pos_class, protein)]
  d <- ggplot2::ggplot(plot_dt_all[order(Pos_class)],ggplot2::aes(b_ddg_pred_class, num_mutations, color = Pos_class)) +
    ggplot2::geom_line(size = 1) +
    # ggplot2::geom_point(size = 0.5) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2) +
    ggplot2::facet_wrap(protein~., scales = "free", ncol = 1) +
    ggplot2::xlab(expression("Binding "*Delta*Delta*"G threshold")) +
    ggplot2::ylab("#Mutations") +
    ggplot2::labs(color = "Residue\nposition") +
    ggplot2::theme_bw()
  if(!is.null(colour_scheme)){
    d <- d + ggplot2::scale_colour_manual(values = unlist(colour_scheme[["shade 0"]][c(1, 3, 4)]))
  }
  ggplot2::ggsave(file.path(outpath, "binding_affinity_mutations_all.pdf"), d, width = 5, height = 7, useDingbats=FALSE)

  #Plot - only PSD95-PDZ3 and GRB2-SH3
  plot_dt_all <- plot_dt[protein!="GB1",.(num_mutations = .N),.(b_ddg_pred_class, Pos_class, protein)]
  d <- ggplot2::ggplot(plot_dt_all[order(Pos_class)],ggplot2::aes(b_ddg_pred_class, num_mutations, color = Pos_class)) +
    ggplot2::geom_line(size = 1) +
    # ggplot2::geom_point(size = 0.5) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2) +
    ggplot2::facet_wrap(~protein, nrow = 1, scales = "free") +
    ggplot2::xlab(expression("Binding "*Delta*Delta*"G threshold")) +
    ggplot2::ylab("#Mutations") +
    ggplot2::labs(color = "Residue\nposition") +
    ggplot2::theme_bw()
  if(!is.null(colour_scheme)){
    d <- d + ggplot2::scale_colour_manual(values = unlist(colour_scheme[["shade 0"]][c(1, 3, 4)]))
  }
  ggplot2::ggsave(file.path(outpath, "binding_affinity_mutations.pdf"), d, width = 7, height = 2, useDingbats=FALSE)

  #Plot - conf - all
  plot_dt_all <- plot_dt_conf[,.(num_mutations = .N),.(b_ddg_pred_class, Pos_class, protein)]
  d <- ggplot2::ggplot(plot_dt_all[order(Pos_class)],ggplot2::aes(b_ddg_pred_class, num_mutations, color = Pos_class)) +
    ggplot2::geom_line(size = 1) +
    # ggplot2::geom_point(size = 0.5) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2) +
    ggplot2::facet_wrap(protein~., scales = "free", ncol = 1) +
    ggplot2::xlab(expression("Binding "*Delta*Delta*"G threshold")) +
    ggplot2::ylab("#Mutations") +
    ggplot2::labs(color = "Residue\nposition") +
    ggplot2::theme_bw()
  if(!is.null(colour_scheme)){
    d <- d + ggplot2::scale_colour_manual(values = unlist(colour_scheme[["shade 0"]][c(1, 3, 4)]))
  }
  ggplot2::ggsave(file.path(outpath, "binding_affinity_mutations_all_conf.pdf"), d, width = 5, height = 7, useDingbats=FALSE)

  #Plot - conf - only PSD95-PDZ3 and GRB2-SH3
  plot_dt_all <- plot_dt_conf[protein!="GB1",.(num_mutations = .N),.(b_ddg_pred_class, Pos_class, protein)]
  d <- ggplot2::ggplot(plot_dt_all[order(Pos_class)],ggplot2::aes(b_ddg_pred_class, num_mutations, color = Pos_class)) +
    ggplot2::geom_line(size = 1) +
    # ggplot2::geom_point(size = 0.5) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2) +
    ggplot2::facet_wrap(~protein, nrow = 1, scales = "free") +
    ggplot2::xlab(expression("Binding "*Delta*Delta*"G threshold")) +
    ggplot2::ylab("#Mutations") +
    ggplot2::labs(color = "Residue\nposition") +
    ggplot2::theme_bw()
  if(!is.null(colour_scheme)){
    d <- d + ggplot2::scale_colour_manual(values = unlist(colour_scheme[["shade 0"]][c(1, 3, 4)]))
  }
  ggplot2::ggsave(file.path(outpath, "binding_affinity_mutations_conf.pdf"), d, width = 7, height = 2, useDingbats=FALSE)

  #Enrichment of charged/polar residues in surface mutations that increase binding affinity
  for(i in dg_dt[,unique(protein)]){
    surf_charge <- dg_dt[b_ddg_pred_conf==T & (b_ddg_pred+1.96*b_ddg_pred_sd)<0 & id!="-0-" & protein==i & Pos_class=="surface", sum(WT_AA %in% unlist(strsplit("RHKDESTNQ", ""))),]
    surf_ncharge <-dg_dt[b_ddg_pred_conf==T & (b_ddg_pred-1.96*b_ddg_pred_sd)>0 & id!="-0-" & protein==i & Pos_class=="surface", sum(WT_AA %in% unlist(strsplit("RHKDESTNQ", ""))),]
    nsurf_charge <- dg_dt[b_ddg_pred_conf==T & (b_ddg_pred+1.96*b_ddg_pred_sd)<0 & id!="-0-" & protein==i & Pos_class=="surface", sum(!WT_AA %in% unlist(strsplit("RHKDESTNQ", ""))),]
    nsurf_ncharge <- dg_dt[b_ddg_pred_conf==T & (b_ddg_pred-1.96*b_ddg_pred_sd)>0 & id!="-0-" & protein==i & Pos_class=="surface", sum(!WT_AA %in% unlist(strsplit("RHKDESTNQ", ""))),]
    temp_test <- fisher.test(matrix(c(surf_charge, surf_ncharge, nsurf_charge, nsurf_ncharge), nrow = 2))
    print(paste0("Enrichment of charged/polar residues in surface mutations that increase binding affinity (", i, "): p-value=", format(temp_test$p.value, digits=2, scientific=T), " odds ratio=", round(temp_test$estimate, 2)))
  }

  #Number of surface and binding interface mutations that increase and decrease binding affinity respectively
  for(i in dg_dt[,unique(protein)]){
    surf_inc <- dg_dt[b_ddg_pred_conf==T & (b_ddg_pred+1.96*b_ddg_pred_sd)<0 & id!="-0-" & protein==i & Pos_class=="surface", .N,]
    bind_dec <- dg_dt[b_ddg_pred_conf==T & (b_ddg_pred-1.96*b_ddg_pred_sd)>0 & id!="-0-" & protein==i & Pos_class=="binding_interface", .N,]
    print(paste0("Number of surface mutations that increase binding affinity (", i, "): ", surf_inc))
    print(paste0("Number of binding interface mutations that decrease binding affinity (", i, "): ", bind_dec))
  }

  ###########################
  ### Save
  ###########################

  #Save
  write.table(dg_dt, 
    file = file.path(outpath, "dg_singles.txt"), 
    quote = F, sep = "\t", row.names = F)

}

