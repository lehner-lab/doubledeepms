
#' doubledeepms_allostery_plots
#'
#' Plot free energy heatmaps.
#'
#' @param input_list path to MoCHI thermo model fit results (required)
#' @param pdb_file_list path to PDB file (required)
#' @param pdb_chain_query_list query chain id (required)
#' @param literature_list literature allosteric sites (required)
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
    temp_dt <- fread(input_list[[protein]])[id!="-0-"]
    temp_dt[, protein := protein]

    #Add WT and mutant AAs
    temp_dt[, WT_AA := substr(id, 1, 1)]
    temp_dt[, Mut := substr(id, nchar(id), nchar(id))]

    #Per residue metrics
    for(i in c("f_ddg", "b_ddg")){
      # temp_dt[,paste0(i, "_posmaxabs") := max(abs(.SD[[1]]), na.rm = T),Pos_ref,.SDcols = paste0(i, c("_pred"))]
      temp_dt[,paste0(i, "_posmedianabs") := median(abs(.SD[[1]]), na.rm = T),Pos_ref,.SDcols = paste0(i, c("_pred"))]
      temp_dt[,paste0(i, "_posmeanabs") := mean(abs(.SD[[1]]), na.rm = T),Pos_ref,.SDcols = paste0(i, c("_pred"))]
      temp_dt[,paste0(i, "_posse") := sd(abs(.SD[[1]]), na.rm = T)/sqrt(sum(!is.na(.SD[[1]]))),Pos_ref,.SDcols = paste0(i, c("_pred"))]
      temp_dt[,paste0(i, "_wposmeanabs") := sum(abs(.SD[[1]])/.SD[[2]]^2, na.rm = T)/sum(1/.SD[[2]]^2, na.rm = T),Pos_ref,.SDcols = paste0(i, c("_pred", "_pred_sd"))]
      temp_dt[,paste0(i, "_wposse") := sqrt(1/sum(1/.SD[[2]]^2, na.rm = T)),Pos_ref,.SDcols = paste0(i, c("_pred", "_pred_sd"))]
    }

    #Per residue metrics - confident only
    for(i in c("f_ddg", "b_ddg")){
      temp_dt[get(paste0(i, "_pred_conf"))==TRUE,paste0(i, "_pred_filtered") := .SD[[1]],,.SDcols = paste0(i, "_pred")]
      temp_dt[get(paste0(i, "_pred_conf"))==TRUE,paste0(i, "_pred_sd_filtered") := .SD[[1]],,.SDcols = paste0(i, "_pred_sd")]
      # temp_dt[,paste0(i, "_posmaxabs_conf") := max(abs(.SD[[1]]), na.rm = T),Pos_ref,.SDcols = paste0(i, c("_pred_filtered"))]
      temp_dt[,paste0(i, "_posmedianabs_conf") := median(abs(.SD[[1]]), na.rm = T),Pos_ref,.SDcols = paste0(i, c("_pred_filtered"))]
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
    b_ddg_pred_med_dt <- dg_dt[protein==i & b_ddg_pred_conf==T & id!="-0-",.(b_new = mean(abs(b_ddg_pred)), resno = Pos_ref),Pos_ref][,.(resno, b_new)]
    old_colnames <- names(pdb_atom_dt)
    pdb_atom_dt <- merge(pdb_atom_dt, b_ddg_pred_med_dt, by = "resno", all.x = T)
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
      input_dt = copy(dg_dt)[protein==i],
      outpath = file.path(outpath, paste0(i, "_persite_folding_energy_vs_distance_scatter.pdf")),
      colour_scheme = colour_scheme,
      trait_name = "folding")
  }

  ###########################
  ### Binding energy distance correlation plots
  ###########################

  for(i in dg_dt[,unique(protein)]){
    allostery_pos <- doubledeepms__persite_energy_vs_distance_plot(
      input_dt = copy(dg_dt)[protein==i],
      literature_sites = literature_list[[i]],
      outpath = file.path(outpath, paste0(i, "_persite_binding_energy_vs_distance_scatter.pdf")),
      colour_scheme = colour_scheme,
      trait_name = "binding")
    dg_dt[protein==i & Pos_ref %in% allostery_pos, allosteric := T]
  }

  ###########################
  ### Binding ROC plot
  ###########################

  for(i in dg_dt[,unique(protein)]){
    doubledeepms__plot_binding_site_ROC(
      input_dt = copy(dg_dt)[protein==i],
      outpath = file.path(outpath, paste0(i, "_binding_site_ROC.pdf")),
      colour_scheme = colour_scheme)
  }

  ###########################
  ### Allosteric mutations
  ###########################

  dg_list <- list()
  for(i in dg_dt[,unique(protein)]){
    dg_list[[i]] <- doubledeepms__allosteric_mutations_scatterplot(
      input_dt = copy(dg_dt)[protein==i],
      outpath = file.path(outpath, paste0(i, "_allosteric_mutations_scatter.pdf")),
      colour_scheme = colour_scheme)
  }
  dg_dt <- rbindlist(dg_list)

  ###########################
  ### Where are strongest binding effects?
  ###########################

  #Thresholded data
  plot_list <- list()
  for(i in c(1:40)/20){
    plot_list[[as.character(i)]] <- dg_dt[b_ddg_pred_conf==T & (b_ddg_pred-1.96*b_ddg_pred_sd)>i & f_ddg_pred_conf==T,]
    plot_list[[as.character(i)]][, b_ddg_pred_class := i]
    plot_list[[as.character(-i)]] <- dg_dt[b_ddg_pred_conf==T & (b_ddg_pred-1.96*b_ddg_pred_sd)<(-i) & f_ddg_pred_conf==T,]
    plot_list[[as.character(-i)]][, b_ddg_pred_class := (-i)]
  }
  plot_dt <- rbindlist(plot_list)
  
  #Plot
  plot_dt_all <- plot_dt[b_ddg_pred_conf==T,.(num_mutations = .N),.(b_ddg_pred_class, Pos_class, protein)]
  d <- ggplot2::ggplot(plot_dt_all[order(Pos_class)],ggplot2::aes(b_ddg_pred_class, num_mutations, color = Pos_class)) +
    ggplot2::geom_line(size = 1) +
    # ggplot2::geom_point(size = 0.5) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2) +
    ggplot2::facet_grid(protein~., scales = "free") +
    ggplot2::xlab(expression("Binding "*Delta*Delta*"G threshold")) +
    ggplot2::ylab("#Mutations") +
    ggplot2::labs(color = "Residue\nposition") +
    ggplot2::theme_bw()
  if(!is.null(colour_scheme)){
    d <- d + ggplot2::scale_colour_manual(values = unlist(colour_scheme[["shade 0"]][c(1, 3, 4)]))
  }
  ggplot2::ggsave(file.path(outpath, "binding_affinity_mutations_all.pdf"), d, width = 7, height = 5, useDingbats=FALSE)

  #Plot - only PSD95-PDZ3 and GRB2-SH3
  plot_dt_all <- plot_dt[b_ddg_pred_conf==T & protein!="GB1",.(num_mutations = .N),.(b_ddg_pred_class, Pos_class, protein)]
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
  ggplot2::ggsave(file.path(outpath, "binding_affinity_mutations.pdf"), d, width = 7, height = 3, useDingbats=FALSE)

  #Plot abundance unchanged
  plot_dt_all <- plot_dt[(abs(f_ddg_pred)-1.96*f_ddg_pred_sd)<0 & b_ddg_pred_conf==T,.(num_mutations = .N),.(b_ddg_pred_class, Pos_class, protein)]
  d <- ggplot2::ggplot(plot_dt_all[order(Pos_class)],ggplot2::aes(b_ddg_pred_class, num_mutations, color = Pos_class)) +
    ggplot2::geom_line(size = 1) +
    # ggplot2::geom_point(size = 0.5) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2) +
    ggplot2::facet_grid(protein~., scales = "free") +
    ggplot2::xlab(expression("Binding "*Delta*Delta*"G threshold")) +
    ggplot2::ylab("#Mutations") +
    ggplot2::labs(color = "Residue\nposition") +
    ggplot2::theme_bw()
  if(!is.null(colour_scheme)){
    d <- d + ggplot2::scale_colour_manual(values = unlist(colour_scheme[["shade 0"]][c(1, 3, 4)]))
  }
  ggplot2::ggsave(file.path(outpath, "binding_affinity_mutations_all_abundanceunchanged.pdf"), d, width = 7, height = 5, useDingbats=FALSE)

  ###########################
  ### Pleiotropy of allosteric mutations
  ###########################

  #Significant change in ddG binding
  dg_dt[b_ddg_pred_conf==T, b_ddg_pred_pvalue := doubledeepms__pvalue(b_ddg_pred, b_ddg_pred_sd)]
  dg_dt[b_ddg_pred_conf==T, b_ddg_pred_FDR := p.adjust(b_ddg_pred_pvalue, method = "BH"),protein]

  #Significant change in ddG folding
  dg_dt[f_ddg_pred_conf==T, f_ddg_pred_pvalue := doubledeepms__pvalue(f_ddg_pred, f_ddg_pred_sd)]
  dg_dt[f_ddg_pred_conf==T, f_ddg_pred_FDR := p.adjust(f_ddg_pred_pvalue, method = "BH"),protein]

  #Pleiotropy class
  dg_dt[, ddg_class := "Remainder"]
  # dg_dt[b_ddg_pred*f_ddg_pred>0 & b_ddg_pred_FDR<0.05 & f_ddg_pred_FDR<0.05, ddg_class := "Synergistic"]
  # dg_dt[b_ddg_pred*f_ddg_pred<0 & b_ddg_pred_FDR<0.05 & f_ddg_pred_FDR<0.05, ddg_class := "Antagonistic"]
  dg_dt[b_ddg_pred_conf==T & b_ddg_pred*f_ddg_pred>0, ddg_class := "Synergistic"]
  dg_dt[b_ddg_pred_conf==T & b_ddg_pred*f_ddg_pred<0, ddg_class := "Antagonistic"]
  # dg_dt[b_ddg_pred_FDR>=0.05 & f_ddg_pred_FDR<0.05, ddg_class := "Folding only"]
  # dg_dt[b_ddg_pred_FDR<0.05) & f_ddg_pred_FDR>=0.05, ddg_class := "Binding only"]

  dg_dt[b_ddg_pred_outlier==T & Pos_class!="binding_interface",.N,.(protein, ddg_class)][order(protein, ddg_class)]
  dg_dt[b_ddg_pred_outlier==T & Pos_class=="binding_interface",.N,.(protein, ddg_class)][order(protein, ddg_class)]

  fisher.test(matrix(c(
    dg_dt[b_ddg_pred_outlier==T & Pos_class=="binding_interface" & ddg_class=="Antagonistic",.N],
    dg_dt[b_ddg_pred_outlier==T & Pos_class=="binding_interface" & ddg_class!="Antagonistic",.N],
    dg_dt[b_ddg_pred_outlier==T & Pos_class!="binding_interface" & ddg_class=="Antagonistic",.N],
    dg_dt[b_ddg_pred_outlier==T & Pos_class!="binding_interface" & ddg_class!="Antagonistic",.N]), nrow = 2))

  # fisher.test(matrix(c(
  #   dg_dt[protein=="PSD95-PDZ3" & b_ddg_pred_outlier==T & Pos_class=="binding_interface" & ddg_class=="Antagonistic",.N],
  #   dg_dt[protein=="PSD95-PDZ3" & b_ddg_pred_outlier==T & Pos_class=="binding_interface" & ddg_class!="Antagonistic",.N],
  #   dg_dt[protein=="PSD95-PDZ3" & b_ddg_pred_outlier==T & Pos_class!="binding_interface" & ddg_class=="Antagonistic",.N],
  #   dg_dt[protein=="PSD95-PDZ3" & b_ddg_pred_outlier==T & Pos_class!="binding_interface" & ddg_class!="Antagonistic",.N]), nrow = 2))

  # fisher.test(matrix(c(
  #   dg_dt[protein=="GRB2-SH3" & b_ddg_pred_outlier==T & Pos_class=="binding_interface" & ddg_class=="Antagonistic",.N],
  #   dg_dt[protein=="GRB2-SH3" & b_ddg_pred_outlier==T & Pos_class=="binding_interface" & ddg_class!="Antagonistic",.N],
  #   dg_dt[protein=="GRB2-SH3" & b_ddg_pred_outlier==T & Pos_class!="binding_interface" & ddg_class=="Antagonistic",.N],
  #   dg_dt[protein=="GRB2-SH3" & b_ddg_pred_outlier==T & Pos_class!="binding_interface" & ddg_class!="Antagonistic",.N]), nrow = 2))

  # fisher.test(matrix(c(
  #   dg_dt[protein=="GB1" & b_ddg_pred_outlier==T & Pos_class=="binding_interface" & ddg_class=="Antagonistic",.N],
  #   dg_dt[protein=="GB1" & b_ddg_pred_outlier==T & Pos_class=="binding_interface" & ddg_class!="Antagonistic",.N],
  #   dg_dt[protein=="GB1" & b_ddg_pred_outlier==T & Pos_class!="binding_interface" & ddg_class=="Antagonistic",.N],
  #   dg_dt[protein=="GB1" & b_ddg_pred_outlier==T & Pos_class!="binding_interface" & ddg_class!="Antagonistic",.N]), nrow = 2))

  dg_dt[allosteric_mutation==T & Pos_class!="binding_interface",.N,.(protein, ddg_class)][order(protein, ddg_class)]
  dg_dt[allosteric_mutation==T & Pos_class=="binding_interface",.N,.(protein, ddg_class)][order(protein, ddg_class)]

  fisher.test(matrix(c(
    dg_dt[allosteric_mutation==T & Pos_class=="binding_interface" & ddg_class=="Antagonistic",.N],
    dg_dt[allosteric_mutation==T & Pos_class=="binding_interface" & ddg_class!="Antagonistic",.N],
    dg_dt[allosteric_mutation==T & Pos_class!="binding_interface" & ddg_class=="Antagonistic",.N],
    dg_dt[allosteric_mutation==T & Pos_class!="binding_interface" & ddg_class!="Antagonistic",.N]), nrow = 2))

  # dg_dt[Pos_class=="binding_interface" & b_ddg_pred_conf==T & f_ddg_pred_conf==T, .(cor = cor(f_ddg_pred, b_ddg_pred)),.(protein, Pos_ref)][,sum(cor<0, na.rm = T)/.N,protein]
  # dg_dt[Pos_class!="binding_interface" & b_ddg_pred_conf==T & f_ddg_pred_conf==T, .(cor = cor(f_ddg_pred, b_ddg_pred)),.(protein, Pos_ref)][,sum(cor<0, na.rm = T)/.N,protein]

  ###########################
  ### Free energy scatterplots
  ###########################

  #Free energy scatterplots by protein - all - conf
  plot_dt <- copy(dg_dt)[,.(protein, f_ddg_pred, b_ddg_pred, f_ddg_pred_conf, b_ddg_pred_conf, Pos_class, id, allosteric, allosteric_mutation)]
  plot_dt <- plot_dt[f_ddg_pred_conf==T & b_ddg_pred_conf==T,]
  plot_dt[, Pos_class_plot := "Remainder"]
  plot_dt[Pos_class=="binding_interface", Pos_class_plot := "binding_interface"]
  plot_dt[allosteric==T & Pos_class!="binding_interface", Pos_class_plot := "allosteric_site"]
  plot_dt[allosteric_mutation==T & Pos_class!="binding_interface", Pos_class_plot := "allosteric_mutation"]
  #Plot
  d <- ggplot2::ggplot(plot_dt[id!="-0-"],ggplot2::aes(f_ddg_pred, b_ddg_pred)) +
    ggplot2::geom_point(alpha = 0.5) +
    ggplot2::geom_smooth(method = "lm", formula = 'y~x') +
    ggplot2::geom_vline(xintercept = 0) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::geom_text(data = plot_dt[,.(label = paste("Pearson's r = ", round(cor(f_ddg_pred, b_ddg_pred, use = "pairwise.complete"), 2), sep="")),.(Pos_class_plot, protein)], ggplot2::aes(label=label, x=-Inf, y=Inf, hjust = 0, vjust = 1)) +
    ggplot2::facet_grid(Pos_class_plot~protein, scales = "free") +
    ggplot2::xlab(expression("Folding "*Delta*Delta*"G")) +
    ggplot2::ylab(expression("Binding "*Delta*Delta*"G")) +
    ggplot2::labs(color = "Residue\nposition") +
    ggplot2::theme_classic()
  if(!is.null(colour_scheme)){
    d <- d + ggplot2::scale_colour_manual(values = unlist(colour_scheme[["shade 0"]][c(1, 3, 4)]))
  }
  ggplot2::ggsave(file.path(outpath, "ddG_scatter_subsets_all.pdf"), d, width = 9, height = 9, useDingbats=FALSE)

}

