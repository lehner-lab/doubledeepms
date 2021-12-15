
#' doubledeepms_polyphen2_comparisons
#'
#' Plot PolyPhen2 comparisons.
#'
#' @param input_file path to input file (required)
#' @param polyphen2_file polyphen2 output file (required)
#' @param position_offset position offset (required)
#' @param uniprot_id uniprot id (required)
#' @param outpath output path for plots and saved objects (required)
#' @param colour_scheme colour scheme file (required)
#' @param execute whether or not to execute the analysis (default: TRUE)
#'
#' @return Nothing
#' @export
#' @import data.table
doubledeepms_polyphen2_comparisons <- function(
  input_file,
  polyphen2_file,
  position_offset,
  uniprot_id,
  outpath,
  colour_scheme,
  execute = TRUE
  ){

  #Return if analysis not executed
  if(!execute){
    return()
  }

  #Display status
  message(paste("\n\n*******", "running stage: doubledeepms_polyphen2_comparisons", "*******\n\n"))

  #Create output directory
  doubledeepms__create_dir(doubledeepms_dir = outpath)

  ### Load single mutant free energies
  ###########################

  #Load dg data
  dg_dt <- fread(input_file)

  ### Format input for polyphen2
  ###########################

  polyphen2_list <- list()
  for(i in names(position_offset)){
    polyphen2_list[[i]] <- dg_dt[protein==i & id!="-0-",.(uniprot_id[[i]], Pos_ref+position_offset[[i]], WT_AA, Mut)]
  }
  polyphen2_dt <- rbindlist(polyphen2_list)

  #Save
  write.table(polyphen2_dt, file = file.path(outpath, "polyphen2_input.txt"), sep = " ", col.names = F, row.names = F, quote = F)

  ### Load polyphen2 results
  ###########################

  suppressWarnings(result_dt <- fread(polyphen2_file, header = T))
  names(result_dt) <- gsub("#", "", names(result_dt))
  for(i in names(position_offset)){
    #Protein
    result_dt[o_acc==uniprot_id[[i]], protein := i]
    #Mutant id
    result_dt[protein == i, id_ref := paste0(o_aa1, o_pos-position_offset[[i]], o_aa2)]
  }

  #Merge
  dg_dt <- merge(dg_dt, result_dt[,.(protein, id_ref, pph2_prob, pph2_FPR, pph2_TPR)], by = c("protein", "id_ref"), all.x = T)

  ### Plot comparisons for protein and position class
  ###########################

  #Confident ddGs
  plot_dt <- dg_dt[!is.na(pph2_prob) & f_ddg_pred_conf==T]
  #Duplicated category of all residues
  plot_dt_all <- copy(plot_dt)
  plot_dt_all[, Pos_class := "all"]
  #Merge
  plot_dt <- rbind(plot_dt, plot_dt_all)
  
  #Plot - all
  cor_dt <- plot_dt[!is.na(Pos_class),.(cor = round(cor(pph2_prob, f_ddg_pred, use = "pairwise.complete"), 2), pph2_prob = Inf, f_ddg_pred = Inf),.(protein, Pos_class)]
  d <- ggplot2::ggplot(plot_dt[!is.na(Pos_class) & !is.na(pph2_prob)],ggplot2::aes(pph2_prob, f_ddg_pred)) +
    ggplot2::stat_binhex(bins = 30, size = 0.2, ggplot2::aes(color = Pos_class)) +
    ggplot2::scale_fill_gradientn(colours = c("white", "black"), trans = "log10") +
    ggplot2::geom_smooth(method = "lm", se = F, color = "black", formula = 'y ~ x', size = 0.5, linetype = 2) + 
    ggplot2::geom_hline(yintercept = 0, size = 0.5) + 
    ggplot2::geom_vline(xintercept = 0, size = 0.5) + 
    ggplot2::xlab(expression("Predicted effect (PolyPhen2)")) +
    ggplot2::ylab(expression("Inferred Folding "*Delta*Delta*"G (ddPCA)")) +
    ggplot2::facet_wrap(protein~Pos_class, scales = "free", nrow = 2) +
    ggplot2::geom_text(data = cor_dt, ggplot2::aes(label=paste("R = ", cor, sep="")), hjust = 1, vjust = 1) +
    ggplot2::theme_bw()
  if(!is.null(colour_scheme)){
    d <- d + ggplot2::scale_colour_manual(values = c("black", unlist(colour_scheme[["shade 0"]][c(1, 3, 4, 4)])))
  }
  ggplot2::ggsave(file.path(outpath, "polyphen2_vs_position_class.pdf"), d, width = 8, height = 4, useDingbats=FALSE)
  
  #Plot - excluding mutations at small/histidine WT residues
  cor_dt <- plot_dt[Pos_class=="all" & !is.na(pph2_prob),.(cor = round(cor(pph2_prob, f_ddg_pred, use = "pairwise.complete"), 2), pph2_prob = Inf, f_ddg_pred = Inf),.(protein, Pos_class)]
  d <- ggplot2::ggplot(plot_dt[Pos_class=="all" & !is.na(pph2_prob)],ggplot2::aes(pph2_prob, f_ddg_pred)) +
    ggplot2::stat_binhex(bins = 30, size = 0.2, color = "grey") +
    ggplot2::scale_fill_gradientn(colours = c("white", "black"), trans = "log10") +
    ggplot2::geom_smooth(method = "lm", se = F, color = "black", formula = 'y ~ x', size = 0.5, linetype = 2) + 
    ggplot2::geom_hline(yintercept = 0, size = 0.5) + 
    ggplot2::geom_vline(xintercept = 0, size = 0.5) + 
    ggplot2::xlab(expression("Predicted effect (PolyPhen2)")) +
    ggplot2::ylab(expression("Inferred Folding "*Delta*Delta*"G (ddPCA)")) +
    ggplot2::facet_wrap(protein~., scales = "free") +
    ggplot2::geom_text(data = cor_dt, ggplot2::aes(label=paste("R = ", cor, sep="")), hjust = 1, vjust = 1) +
    ggplot2::theme_bw()
  ggplot2::ggsave(file.path(outpath, "polyphen2_vs_position_class_all.pdf"), d, width = 6, height = 3, useDingbats=FALSE)

}

