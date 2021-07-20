
#' doubledeepms_allostery_plots
#'
#' Plot free energy heatmaps.
#'
#' @param input_list path to MoCHI thermo model fit results (required)
#' @param pdb_file_list path to PDB file (required)
#' @param pdb_chain_query_list query chain id (required)
#' @param literature_list literature allosteric sites (required)
#' @param ohm_file_list ohm output file list (required)
#' @param aaprop_file path to amino acid properties file (required)
#' @param aaprop_file_selected path to file with selected subset of identifiers
#' @param outpath output path for plots and saved objects (required)
#' @param colour_scheme colour scheme file (required)
#' @param execute whether or not to execute the analysis (default: TRUE)
#'
#' @return Nothing
#' @export
#' @import data.table
doubledeepms_allostery_plots <- function(
  input_list,
  pdb_file_list,
  pdb_chain_query_list,
  literature_list,
  ohm_file_list,
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
  message(paste("\n\n*******", "running stage: doubledeepms_allostery_plots", "*******\n\n"))

  #Create output directory
  doubledeepms__create_dir(doubledeepms_dir = outpath)

  ### Load single mutant free energies
  ###########################

  #Load dg data
  dg_list <- list()
  for(protein in names(input_list)){
    temp_dt <- fread(input_list[[protein]])
    temp_dt[, protein := protein]

    #Add WT and mutant AAs
    temp_dt[id!="-0-", WT_AA := substr(id, 1, 1)]
    temp_dt[id!="-0-", Mut := substr(id, nchar(id), nchar(id))]

    #Per residue metrics
    for(i in c("f_ddg", "b_ddg")){
      temp_dt[,paste0(i, "_posmeanabs") := mean(abs(.SD[[1]]), na.rm = T),Pos_ref,.SDcols = paste0(i, c("_pred"))]
      temp_dt[,paste0(i, "_posse") := sd(abs(.SD[[1]]), na.rm = T)/sqrt(sum(!is.na(.SD[[1]]))),Pos_ref,.SDcols = paste0(i, c("_pred"))]
      temp_dt[,paste0(i, "_wposmeanabs") := sum(abs(.SD[[1]])/.SD[[2]]^2, na.rm = T)/sum(1/.SD[[2]]^2, na.rm = T),Pos_ref,.SDcols = paste0(i, c("_pred", "_pred_sd"))]
      temp_dt[,paste0(i, "_wposse") := sqrt(1/sum(1/.SD[[2]]^2, na.rm = T)),Pos_ref,.SDcols = paste0(i, c("_pred", "_pred_sd"))]
    }

    #Per residue metrics - confident only
    for(i in c("f_ddg", "b_ddg")){
      temp_dt[get(paste0(i, "_pred_conf"))==TRUE,paste0(i, "_pred_filtered") := .SD[[1]],,.SDcols = paste0(i, "_pred")]
      temp_dt[get(paste0(i, "_pred_conf"))==TRUE,paste0(i, "_pred_sd_filtered") := .SD[[1]],,.SDcols = paste0(i, "_pred_sd")]
      temp_dt[,paste0(i, "_posmeanabs_conf") := mean(abs(.SD[[1]]), na.rm = T),Pos_ref,.SDcols = paste0(i, c("_pred_filtered"))]
      temp_dt[,paste0(i, "_posse_conf") := sd(abs(.SD[[1]]), na.rm = T)/sqrt(sum(!is.na(.SD[[1]]))),Pos_ref,.SDcols = paste0(i, c("_pred_filtered"))]
      temp_dt[,paste0(i, "_wposmeanabs_conf") := sum(abs(.SD[[1]])/.SD[[2]]^2, na.rm = T)/sum(1/.SD[[2]]^2, na.rm = T),Pos_ref,.SDcols = paste0(i, c("_pred_filtered", "_pred_sd_filtered"))]
      temp_dt[,paste0(i, "_wposse_conf") := sqrt(1/sum(1/.SD[[2]]^2, na.rm = T)),Pos_ref,.SDcols = paste0(i, c("_pred_filtered", "_pred_sd_filtered"))]
    }

    # #Mann whitney U test + randomisation
    # result_list <- list()
    # for(i in temp_dt[order(Pos_ref), unique(Pos_ref)]){
    #   print(i)
    #   result_list[[as.character(i)]] <- doubledeepms__mann_whitney_U_wrapper_rand(
    #     temp_dt[Pos_ref==i & b_ddg_pred_conf,abs(b_ddg_pred)],
    #     temp_dt[Pos_ref!=i & b_ddg_pred_conf,abs(b_ddg_pred)],
    #     temp_dt[Pos_ref==i & b_ddg_pred_conf,b_ddg_pred_sd],
    #     temp_dt[Pos_ref!=i & b_ddg_pred_conf,b_ddg_pred_sd],
    #     100)
    # }
    # result_dt <- as.data.table(do.call('rbind', result_list))
    # result_dt[, Pos_ref := as.integer(names(result_list))]
    # temp_dt <- merge(temp_dt, result_dt, by = "Pos_ref")

    dg_list[[protein]] <- temp_dt
  }
  dg_dt <- rbindlist(dg_list)

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
  ### Binding energy distance correlation plots
  ###########################

  for(i in dg_dt[,unique(protein)]){
    allostery_pos <- doubledeepms__persite_energy_vs_distance_plot(
      input_dt = copy(dg_dt)[protein==i & id!="-0-"],
      literature_sites = literature_list[[i]][["class_switching"]],
      outpath = file.path(outpath, paste0(i, "_persite_binding_energy_vs_distance_scatter.pdf")),
      colour_scheme = colour_scheme,
      trait_name = "binding")
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
      outpath = file.path(outpath, paste0(i, "_binding_site_ROC.pdf")),
      colour_scheme = colour_scheme)
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
  cor_dt <- plot_dt[,.(cor = round(cor(scHAmin_ligand, prop_mut, use = "pairwise.complete"), 2), scHAmin_ligand = 12, prop_mut = 100),protein]
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(scHAmin_ligand, prop_mut)) +
    ggplot2::geom_smooth(method = "lm", se = F, color = "grey", linetype = 2, formula = 'y ~ x') + 
    ggplot2::geom_point(ggplot2::aes(shape = !is.na(allosteric), color = Pos_class), size = 2) +
    ggplot2::xlab(expression("Distance to ligand ("*ring(A)*")")) +
    ggplot2::ylab("%Allosteric mutations per residue") +
    ggplot2::facet_wrap(protein~., scales = "free", ncol = 1) +
    ggplot2::labs(color = "Residue\nposition", shape = "Major allosteric\nsite") +
    ggplot2::geom_text(data = cor_dt, ggplot2::aes(label=paste("Pearson's r = ", cor, sep=""))) +
    ggplot2::theme_bw()
  if(!is.null(colour_scheme)){
    d <- d + ggplot2::scale_colour_manual(values = unlist(colour_scheme[["shade 0"]][c(3, 4)]))
  }
  ggplot2::ggsave(file.path(outpath, "allosteric_mutations_vs_distance.pdf"), d, width = 4, height = 6, useDingbats=FALSE)

  ###########################
  ### Proportion of mutations in allosteric sites that are allosteric mutations
  ###########################

  #Sector residues
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
      in_sector_means <- dg_dt[protein==i & !is.na(b_ddg_wposmeanabs)][!duplicated(Pos_ref)][Pos_ref %in% literature_list[[i]][["sector"]],b_ddg_wposmeanabs]
      out_sector_means <- dg_dt[protein==i & !is.na(b_ddg_wposmeanabs)][!duplicated(Pos_ref)][!Pos_ref %in% literature_list[[i]][["sector"]],b_ddg_wposmeanabs]
      temp_test <- doubledeepms__mann_whitney_U_wrapper(in_sector_means, out_sector_means)
      print(paste0("Predicting sector residues using position-wise weighted mean of binding free energies changes (", i, "): p-value=", format(temp_test[["p_value"]], digits=2, scientific=T), " AUC=", round(temp_test[["effect_size"]], 2)))
   }
  }

  #Class switching residues
  for(i in dg_dt[,unique(protein)]){
    if(length(literature_list[[i]][["class_switching"]]!=0)){
      in_sector_means <- dg_dt[protein==i & !is.na(b_ddg_wposmeanabs)][!duplicated(Pos_ref)][Pos_ref %in% literature_list[[i]][["class_switching"]],b_ddg_wposmeanabs]
      out_sector_means <- dg_dt[protein==i & !is.na(b_ddg_wposmeanabs)][!duplicated(Pos_ref)][!Pos_ref %in% literature_list[[i]][["class_switching"]],b_ddg_wposmeanabs]
      temp_test <- doubledeepms__mann_whitney_U_wrapper(in_sector_means, out_sector_means)
      print(paste0("Predicting class-switching residues using position-wise weighted mean of binding free energies changes (", i, "): p-value=", format(temp_test[["p_value"]], digits=2, scientific=T), " AUC=", round(temp_test[["effect_size"]], 2)))
      #Predicting class-switching residues not classified as allosteric
      in_sector_means <- dg_dt[protein==i & !is.na(b_ddg_pred)][Pos_ref %in% literature_list[[i]][["class_switching"]] & is.na(allosteric) & is.na(orthosteric),abs(b_ddg_pred)]
      out_sector_means <- dg_dt[protein==i & !is.na(b_ddg_pred)][!Pos_ref %in% literature_list[[i]][["class_switching"]] & is.na(allosteric) & is.na(orthosteric),abs(b_ddg_pred)]
      temp_test <- doubledeepms__mann_whitney_U_wrapper(in_sector_means, out_sector_means)
      print(paste0("Predicting class-switching residues using binding free energies changes of mutations at non-allosteric/orthosteric sites (", i, "): p-value=", format(temp_test[["p_value"]], digits=2, scientific=T), " AUC=", round(temp_test[["effect_size"]], 2)))
    }
  }

  ###########################
  ### Enrichment of allosteric mutations in literature sites
  ###########################

  #Sector residues
  for(i in dg_dt[,unique(protein)]){
    if(length(literature_list[[i]][["sector"]]!=0)){
      in_sector_allo <- dg_dt[protein==i & Pos_class!="binding_interface"][Pos_ref %in% literature_list[[i]][["sector"]] & allosteric_mutation==T,.N]
      out_sector_allo <- dg_dt[protein==i & Pos_class!="binding_interface"][!Pos_ref %in% literature_list[[i]][["sector"]] & allosteric_mutation==T,.N]
      in_sector_nallo <- dg_dt[protein==i & Pos_class!="binding_interface"][Pos_ref %in% literature_list[[i]][["sector"]] & is.na(allosteric_mutation),.N]
      out_sector_nallo <- dg_dt[protein==i & Pos_class!="binding_interface"][!Pos_ref %in% literature_list[[i]][["sector"]] & is.na(allosteric_mutation),.N]
      temp_test <- fisher.test(matrix(c(in_sector_allo, out_sector_allo, in_sector_nallo, out_sector_nallo), nrow = 2))
      print(paste0("Enrichment of allosteric mutations in sector (", i, "): p-value=", format(temp_test$p.value, digits=2, scientific=T), " odds ratio=", round(temp_test$estimate, 2)))
    }
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

