
#' doubledeepms_foldx_comparisons
#'
#' Plot FoldX comparisons.
#'
#' @param input_file path to input file (required)
#' @param foldx_file_list foldx output file list (required)
#' @param outpath output path for plots and saved objects (required)
#' @param colour_scheme colour scheme file (required)
#' @param execute whether or not to execute the analysis (default: TRUE)
#'
#' @return Nothing
#' @export
#' @import data.table
doubledeepms_foldx_comparisons <- function(
  input_file,
  foldx_file_list,
  outpath,
  colour_scheme,
  execute = TRUE
  ){

  #Return if analysis not executed
  if(!execute){
    return()
  }

  #Display status
  message(paste("\n\n*******", "running stage: doubledeepms_foldx_comparisons", "*******\n\n"))

  #Create output directory
  doubledeepms__create_dir(doubledeepms_dir = outpath)

  ### Load single mutant free energies
  ###########################

  #Load dg data
  dg_dt <- fread(input_file)

  ### Weighted mean folding free energy
  ###########################

  #Significant stabilising mutations
  dg_dt[f_ddg_pred_conf==T & id!="-0-", f_ddg_pred_fdr := p.adjust(doubledeepms__pvalue(f_ddg_pred, f_ddg_pred_sd), method = "BH"),.(protein)]

  #Stabilising mutations
  dg_dt[, f_ddg_pred_stab := 0]
  dg_dt[f_ddg_pred<0 & f_ddg_pred_fdr<0.05, f_ddg_pred_stab := 1]

  #Stabilising residues (at least 5)
  dg_dt[, f_ddg_pred_stab_res5 := F]
  for(i in names(foldx_file_list)){
    temp <- dg_dt[protein==i & f_ddg_pred_stab==1,table(Pos_ref)]
    stab_res <- as.integer(names(temp)[temp>=5])
    if(length(stab_res)==0){
      stab_res <- as.integer(names(temp))
      dg_dt[protein==i & Pos_ref %in% stab_res, f_ddg_pred_stab_res5 := T]
      # print(paste0("Destabilising residues (1+) for ", i, ": ", paste(dg_dt[protein==i & f_ddg_pred_stab_res5][!duplicated(Pos_ref),Pos_ref], collapse = ",")))
      # print(paste0("Surface destabilising (1+) residues for ", i, ": ", paste(dg_dt[protein==i & f_ddg_pred_stab_res5 & Pos_class=="surface"][!duplicated(Pos_ref),Pos_ref], collapse = ",")))
    }else{
      dg_dt[protein==i & Pos_ref %in% stab_res, f_ddg_pred_stab_res5 := T]
      # print(paste0("Destabilising residues (3+) for ", i, ": ", paste(dg_dt[protein==i & f_ddg_pred_stab_res5][!duplicated(Pos_ref),Pos_ref], collapse = ",")))
      # print(paste0("Surface destabilising (3+) residues for ", i, ": ", paste(dg_dt[protein==i & f_ddg_pred_stab_res5 & Pos_class=="surface"][!duplicated(Pos_ref),Pos_ref], collapse = ",")))
    }
  }

  ### Load foldx results
  ###########################

  result_list <- list()
  for(i in names(foldx_file_list)){
    #Load foldx results
    f_dt <- fread(foldx_file_list[[i]])
    names(f_dt) <- c("mut", "f_ddg_pred_foldx")
    f_dt[, id_mut := substr(mut, 5, nchar(mut))]
    #Mutant id
    dg_dt[protein==i, id_mut := substr(id_ref, 2, nchar(id_ref))]
    #Merge
    result_list[[i]] <- merge(dg_dt[protein==i], f_dt[!duplicated(id_mut),.(id_mut, f_ddg_pred_foldx)], by = "id_mut", all = T)
  }
  dg_dt <- rbindlist(result_list)

  ### Plot comparisons for protein and position class
  ###########################

  #Confident ddGs
  plot_dt <- dg_dt[f_ddg_pred_conf==T]
  #AA type
  plot_dt[, AA_type := "Remainder"]
  plot_dt[WT_AA %in% c("G", "A", "V"), AA_type := "Small"]
  plot_dt[WT_AA %in% c("H"), AA_type := "Histidine"]
  #Duplicated category of all residues
  plot_dt_all <- copy(plot_dt)
  plot_dt_all[, Pos_class := "all"]
  #Merge
  plot_dt <- rbind(plot_dt, plot_dt_all)
  #Duplicated category of surface destabilising residues
  plot_dt_sd <- plot_dt[f_ddg_pred_stab_res5==T & Pos_class=="surface"]
  plot_dt_sd[, Pos_class := "surface_destab"]
  #Merge
  plot_dt <- rbind(plot_dt, plot_dt_sd)
  
  #Plot - all
  cor_dt <- plot_dt[!is.na(Pos_class),.(cor = round(cor(f_ddg_pred_foldx, f_ddg_pred, use = "pairwise.complete"), 2), f_ddg_pred_foldx = Inf, f_ddg_pred = Inf),.(protein, Pos_class)]
  d <- ggplot2::ggplot(plot_dt[!is.na(Pos_class) & !is.na(f_ddg_pred_foldx)],ggplot2::aes(f_ddg_pred_foldx, f_ddg_pred)) +
    ggplot2::stat_binhex(bins = 30, size = 0.2, ggplot2::aes(color = Pos_class)) +
    ggplot2::scale_fill_gradientn(colours = c("white", "black"), trans = "log10") +
    ggplot2::geom_smooth(method = "lm", se = F, color = "black", formula = 'y ~ x', size = 0.5, linetype = 2) + 
    ggplot2::geom_hline(yintercept = 0, size = 0.5) + 
    ggplot2::geom_vline(xintercept = 0, size = 0.5) + 
    ggplot2::xlab(expression("Predicted Folding "*Delta*Delta*"G (FoldX)")) +
    ggplot2::ylab(expression("Inferred Folding "*Delta*Delta*"G (ddPCA)")) +
    ggplot2::facet_wrap(protein~Pos_class, scales = "free", nrow = 3) +
    ggplot2::geom_text(data = cor_dt, ggplot2::aes(label=paste("R = ", cor, sep="")), hjust = 1, vjust = 1) +
    ggplot2::theme_bw()
  if(!is.null(colour_scheme)){
    d <- d + ggplot2::scale_colour_manual(values = c("black", unlist(colour_scheme[["shade 0"]][c(1, 3, 4, 4)])))
  }
  ggplot2::ggsave(file.path(outpath, "foldx_vs_position_class.pdf"), d, width = 9, height = 6, useDingbats=FALSE)
  
  #Plot - excluding mutations at small/histidine WT residues
  cor_dt <- plot_dt[!is.na(Pos_class) & AA_type=="Remainder",.(cor = round(cor(f_ddg_pred_foldx, f_ddg_pred, use = "pairwise.complete"), 2), f_ddg_pred_foldx = Inf, f_ddg_pred = Inf),.(protein, Pos_class)]
  d <- ggplot2::ggplot(plot_dt[!is.na(Pos_class) & AA_type=="Remainder" & !is.na(f_ddg_pred_foldx)],ggplot2::aes(f_ddg_pred_foldx, f_ddg_pred)) +
    ggplot2::stat_binhex(bins = 30, size = 0.2, ggplot2::aes(color = Pos_class)) +
    ggplot2::scale_fill_gradientn(colours = c("white", "black"), trans = "log10") +
    ggplot2::geom_smooth(method = "lm", se = F, color = "black", formula = 'y ~ x', size = 0.5, linetype = 2) + 
    ggplot2::geom_hline(yintercept = 0, size = 0.5) + 
    ggplot2::geom_vline(xintercept = 0, size = 0.5) + 
    ggplot2::xlab(expression("Predicted Folding "*Delta*Delta*"G (FoldX)")) +
    ggplot2::ylab(expression("Inferred Folding "*Delta*Delta*"G (ddPCA)")) +
    ggplot2::facet_wrap(protein~Pos_class, scales = "free", nrow = 3) +
    ggplot2::geom_text(data = cor_dt, ggplot2::aes(label=paste("R = ", cor, sep="")), hjust = 1, vjust = 1) +
    ggplot2::theme_bw()
  if(!is.null(colour_scheme)){
    d <- d + ggplot2::scale_colour_manual(values = c("black", unlist(colour_scheme[["shade 0"]][c(1, 3, 4, 4)])))
  }
  ggplot2::ggsave(file.path(outpath, "foldx_vs_position_class_remainder.pdf"), d, width = 9, height = 6, useDingbats=FALSE)

  #Plot - excluding mutations at small/histidine WT residues
  cor_dt <- plot_dt[Pos_class=="all" & AA_type=="Remainder",.(cor = round(cor(f_ddg_pred_foldx, f_ddg_pred, use = "pairwise.complete"), 2), f_ddg_pred_foldx = Inf, f_ddg_pred = Inf),.(protein, Pos_class)]
  d <- ggplot2::ggplot(plot_dt[Pos_class=="all" & AA_type=="Remainder" & !is.na(f_ddg_pred_foldx)],ggplot2::aes(f_ddg_pred_foldx, f_ddg_pred)) +
    ggplot2::stat_binhex(bins = 30, size = 0.2, color = "grey") +
    ggplot2::scale_fill_gradientn(colours = c("white", "black"), trans = "log10") +
    ggplot2::geom_smooth(method = "lm", se = F, color = "black", formula = 'y ~ x', size = 0.5, linetype = 2) + 
    ggplot2::geom_hline(yintercept = 0, size = 0.5) + 
    ggplot2::geom_vline(xintercept = 0, size = 0.5) + 
    ggplot2::xlab(expression("Predicted Folding "*Delta*Delta*"G (FoldX)")) +
    ggplot2::ylab(expression("Inferred Folding "*Delta*Delta*"G (ddPCA)")) +
    ggplot2::facet_wrap(protein~., scales = "free") +
    ggplot2::geom_text(data = cor_dt, ggplot2::aes(label=paste("R = ", cor, sep="")), hjust = 1, vjust = 1) +
    ggplot2::theme_bw()
  ggplot2::ggsave(file.path(outpath, "foldx_vs_position_class_remainder_all.pdf"), d, width = 8, height = 3, useDingbats=FALSE)

  print(paste0("Correlation of FoldX ddGs with high confidence inferred folding ddGs (all): ", plot_dt[Pos_class=="all",round(cor(f_ddg_pred_foldx, f_ddg_pred, use = "pairwise.complete"), 2)]))
  print(paste0("Correlation of FoldX ddGs with high confidence inferred folding ddGs (excluding substitutions at HGAV): ", plot_dt[Pos_class=="all" & AA_type=="Remainder",round(cor(f_ddg_pred_foldx, f_ddg_pred, use = "pairwise.complete"), 2)]))

  print("Agreement with FoldX predictions for stabilising mutations at surface destabilising residues")
  print(plot_dt[Pos_class=="surface_destab" & f_ddg_pred<0 & f_ddg_pred_fdr<0.05,.(foldx_agree = sum(f_ddg_pred_foldx<0, na.rm = T), foldx_disagree = sum(f_ddg_pred_foldx>=0, na.rm = T)),protein])

  print("Agreement with FoldX predictions for surface stabilising mutations")
  print(plot_dt[Pos_class=="surface" & f_ddg_pred<0 & f_ddg_pred_fdr<0.05,.(foldx_agree = sum(f_ddg_pred_foldx<0, na.rm = T), foldx_disagree = sum(f_ddg_pred_foldx>=0, na.rm = T)),protein])

  print("Agreement with FoldX predictions for stabilising mutations")
  print(plot_dt[Pos_class=="all" & f_ddg_pred<0 & f_ddg_pred_fdr<0.05,.(foldx_agree = sum(f_ddg_pred_foldx<0, na.rm = T), foldx_disagree = sum(f_ddg_pred_foldx>=0, na.rm = T)),protein])



}

