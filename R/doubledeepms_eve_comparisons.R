
#' doubledeepms_eve_comparisons
#'
#' Plot EVE comparisons.
#'
#' @param eve_file_list EVE output file list (required)
#' @param outpath output path for plots and saved objects (required)
#' @param colour_scheme colour scheme file (required)
#' @param execute whether or not to execute the analysis (default: TRUE)
#'
#' @return Nothing
#' @export
#' @import data.table
doubledeepms_eve_comparisons <- function(
  eve_file_list,
  outpath,
  colour_scheme,
  execute = TRUE
  ){

  #Return if analysis not executed
  if(!execute){
    return()
  }

  #Display status
  message(paste("\n\n*******", "running stage: doubledeepms_eve_comparisons", "*******\n\n"))

  #Create output directory
  doubledeepms__create_dir(doubledeepms_dir = outpath)

  ### Load EVE results
  ###########################

  #Load data for all domains
  dg_dt <- rbindlist(lapply(eve_file_list, fread))

  ### Plot EVE comparisons for protein and position class - Folding
  ###########################

  #Confident ddGs
  plot_dt <- dg_dt[!is.na(prediction_EVE_evol_indices) & f_ddg_pred_conf==T]
  #Duplicated category of all residues
  plot_dt_all <- copy(plot_dt)
  plot_dt_all[, Pos_class := "all"]
  #Merge
  plot_dt <- rbind(plot_dt, plot_dt_all)
  
  #Plot - all
  cor_dt <- plot_dt[!is.na(Pos_class),.(cor = round(cor(prediction_EVE_evol_indices, f_ddg_pred, use = "pairwise.complete"), 2), prediction_EVE_evol_indices = Inf, f_ddg_pred = Inf),.(protein, Pos_class)]
  d <- ggplot2::ggplot(plot_dt[!is.na(Pos_class) & !is.na(prediction_EVE_evol_indices)],ggplot2::aes(prediction_EVE_evol_indices, f_ddg_pred)) +
    ggplot2::stat_binhex(bins = 30, size = 0.2, ggplot2::aes(color = Pos_class)) +
    ggplot2::scale_fill_gradientn(colours = c("white", "black"), trans = "log10") +
    ggplot2::geom_smooth(method = "lm", se = F, color = "black", formula = 'y ~ x', size = 0.5, linetype = 2) + 
    ggplot2::geom_hline(yintercept = 0, size = 0.5) + 
    ggplot2::geom_vline(xintercept = 0, size = 0.5) + 
    ggplot2::xlab(expression("Predicted effect (EVE pathogenicity score)")) +
    ggplot2::ylab(expression("Inferred Folding "*Delta*Delta*"G (ddPCA)")) +
    ggplot2::facet_wrap(protein~Pos_class, scales = "free", nrow = 3) +
    ggplot2::geom_text(data = cor_dt, ggplot2::aes(label=paste("R = ", cor, sep="")), hjust = 1, vjust = 1) +
    ggplot2::theme_bw()
  if(!is.null(colour_scheme)){
    d <- d + ggplot2::scale_colour_manual(values = c("black", unlist(colour_scheme[["shade 0"]][c(1, 3, 4, 4)])))
  }
  ggplot2::ggsave(file.path(outpath, "eve_vs_f_ddg_position_class.pdf"), d, width = 8, height = 4, useDingbats=FALSE)
  
  #Plot - excluding mutations at small/histidine WT residues
  cor_dt <- plot_dt[Pos_class=="all" & !is.na(prediction_EVE_evol_indices),.(cor = round(cor(prediction_EVE_evol_indices, f_ddg_pred, use = "pairwise.complete"), 2), prediction_EVE_evol_indices = Inf, f_ddg_pred = Inf),.(protein, Pos_class)]
  d <- ggplot2::ggplot(plot_dt[Pos_class=="all" & !is.na(prediction_EVE_evol_indices)],ggplot2::aes(prediction_EVE_evol_indices, f_ddg_pred)) +
    ggplot2::stat_binhex(bins = 30, size = 0.2, color = "grey") +
    ggplot2::scale_fill_gradientn(colours = c("white", "black"), trans = "log10") +
    ggplot2::geom_smooth(method = "lm", se = F, color = "black", formula = 'y ~ x', size = 0.5, linetype = 2) + 
    ggplot2::geom_hline(yintercept = 0, size = 0.5) + 
    ggplot2::geom_vline(xintercept = 0, size = 0.5) + 
    ggplot2::xlab(expression("Predicted effect (EVE pathogenicity score)")) +
    ggplot2::ylab(expression("Inferred Folding "*Delta*Delta*"G (ddPCA)")) +
    ggplot2::facet_wrap(protein~., scales = "free") +
    ggplot2::geom_text(data = cor_dt, ggplot2::aes(label=paste("R = ", cor, sep="")), hjust = 1, vjust = 1) +
    ggplot2::theme_bw()
  ggplot2::ggsave(file.path(outpath, "eve_vs_f_ddg_position_class_all.pdf"), d, width = 8, height = 3, useDingbats=FALSE)

  ### Plot EVmutation comparisons for protein and position class - Folding
  ###########################

  #Confident ddGs
  plot_dt <- dg_dt[!is.na(prediction_epistatic) & f_ddg_pred_conf==T]
  #Duplicated category of all residues
  plot_dt_all <- copy(plot_dt)
  plot_dt_all[, Pos_class := "all"]
  #Merge
  plot_dt <- rbind(plot_dt, plot_dt_all)
  
  #Plot - all
  cor_dt <- plot_dt[!is.na(Pos_class),.(cor = round(cor(-prediction_epistatic, f_ddg_pred, use = "pairwise.complete"), 2), prediction_epistatic = Inf, f_ddg_pred = Inf),.(protein, Pos_class)]
  d <- ggplot2::ggplot(plot_dt[!is.na(Pos_class) & !is.na(prediction_epistatic)],ggplot2::aes(-prediction_epistatic, f_ddg_pred)) +
    ggplot2::stat_binhex(bins = 30, size = 0.2, ggplot2::aes(color = Pos_class)) +
    ggplot2::scale_fill_gradientn(colours = c("white", "black"), trans = "log10") +
    ggplot2::geom_smooth(method = "lm", se = F, color = "black", formula = 'y ~ x', size = 0.5, linetype = 2) + 
    ggplot2::geom_hline(yintercept = 0, size = 0.5) + 
    ggplot2::geom_vline(xintercept = 0, size = 0.5) + 
    ggplot2::xlab(expression("Predicted effect (-EVmutation score)")) +
    ggplot2::ylab(expression("Inferred Folding "*Delta*Delta*"G (ddPCA)")) +
    ggplot2::facet_wrap(protein~Pos_class, scales = "free", nrow = 3) +
    ggplot2::geom_text(data = cor_dt, ggplot2::aes(label=paste("R = ", cor, sep="")), hjust = 0, vjust = 1) +
    ggplot2::theme_bw()
  if(!is.null(colour_scheme)){
    d <- d + ggplot2::scale_colour_manual(values = c("black", unlist(colour_scheme[["shade 0"]][c(1, 3, 4, 4)])))
  }
  ggplot2::ggsave(file.path(outpath, "evemutation_vs_f_ddg_position_class.pdf"), d, width = 8, height = 4, useDingbats=FALSE)
  
  #Plot - excluding mutations at small/histidine WT residues
  cor_dt <- plot_dt[Pos_class=="all" & !is.na(prediction_epistatic),.(cor = round(cor(-prediction_epistatic, f_ddg_pred, use = "pairwise.complete"), 2), prediction_epistatic = Inf, f_ddg_pred = Inf),.(protein, Pos_class)]
  d <- ggplot2::ggplot(plot_dt[Pos_class=="all" & !is.na(prediction_epistatic)],ggplot2::aes(-prediction_epistatic, f_ddg_pred)) +
    ggplot2::stat_binhex(bins = 30, size = 0.2, color = "grey") +
    ggplot2::scale_fill_gradientn(colours = c("white", "black"), trans = "log10") +
    ggplot2::geom_smooth(method = "lm", se = F, color = "black", formula = 'y ~ x', size = 0.5, linetype = 2) + 
    ggplot2::geom_hline(yintercept = 0, size = 0.5) + 
    ggplot2::geom_vline(xintercept = 0, size = 0.5) + 
    ggplot2::xlab(expression("Predicted effect (-EVmutation score)")) +
    ggplot2::ylab(expression("Inferred Folding "*Delta*Delta*"G (ddPCA)")) +
    ggplot2::facet_wrap(protein~., scales = "free") +
    ggplot2::geom_text(data = cor_dt, ggplot2::aes(label=paste("R = ", cor, sep="")), hjust = 0, vjust = 1) +
    ggplot2::theme_bw()
  ggplot2::ggsave(file.path(outpath, "evemutation_vs_f_ddg_position_class_all.pdf"), d, width = 8, height = 3, useDingbats=FALSE)

  ### Plot EVE comparisons for protein and position class - Binding
  ###########################

  #Confident ddGs
  plot_dt <- dg_dt[!is.na(prediction_EVE_evol_indices) & b_ddg_pred_conf==T]
  #Duplicated category of all residues
  plot_dt_all <- copy(plot_dt)
  plot_dt_all[, Pos_class := "all"]
  #Merge
  plot_dt <- rbind(plot_dt, plot_dt_all)
  
  #Plot - all
  cor_dt <- plot_dt[!is.na(Pos_class),.(cor = round(cor(prediction_EVE_evol_indices, b_ddg_pred, use = "pairwise.complete"), 2), prediction_EVE_evol_indices = Inf, b_ddg_pred = Inf),.(protein, Pos_class)]
  d <- ggplot2::ggplot(plot_dt[!is.na(Pos_class) & !is.na(prediction_EVE_evol_indices)],ggplot2::aes(prediction_EVE_evol_indices, b_ddg_pred)) +
    ggplot2::stat_binhex(bins = 30, size = 0.2, ggplot2::aes(color = Pos_class)) +
    ggplot2::scale_fill_gradientn(colours = c("white", "black"), trans = "log10") +
    ggplot2::geom_smooth(method = "lm", se = F, color = "black", formula = 'y ~ x', size = 0.5, linetype = 2) + 
    ggplot2::geom_hline(yintercept = 0, size = 0.5) + 
    ggplot2::geom_vline(xintercept = 0, size = 0.5) + 
    ggplot2::xlab(expression("Predicted effect (EVE pathogenicity score)")) +
    ggplot2::ylab(expression("Inferred Folding "*Delta*Delta*"G (ddPCA)")) +
    ggplot2::facet_wrap(protein~Pos_class, scales = "free", nrow = 3) +
    ggplot2::geom_text(data = cor_dt, ggplot2::aes(label=paste("R = ", cor, sep="")), hjust = 1, vjust = 1) +
    ggplot2::theme_bw()
  if(!is.null(colour_scheme)){
    d <- d + ggplot2::scale_colour_manual(values = c("black", unlist(colour_scheme[["shade 0"]][c(1, 3, 4, 4)])))
  }
  ggplot2::ggsave(file.path(outpath, "eve_vs_b_ddg_position_class.pdf"), d, width = 8, height = 4, useDingbats=FALSE)
  
  #Plot - excluding mutations at small/histidine WT residues
  cor_dt <- plot_dt[Pos_class=="all" & !is.na(prediction_EVE_evol_indices),.(cor = round(cor(prediction_EVE_evol_indices, b_ddg_pred, use = "pairwise.complete"), 2), prediction_EVE_evol_indices = Inf, b_ddg_pred = Inf),.(protein, Pos_class)]
  d <- ggplot2::ggplot(plot_dt[Pos_class=="all" & !is.na(prediction_EVE_evol_indices)],ggplot2::aes(prediction_EVE_evol_indices, b_ddg_pred)) +
    ggplot2::stat_binhex(bins = 30, size = 0.2, color = "grey") +
    ggplot2::scale_fill_gradientn(colours = c("white", "black"), trans = "log10") +
    ggplot2::geom_smooth(method = "lm", se = F, color = "black", formula = 'y ~ x', size = 0.5, linetype = 2) + 
    ggplot2::geom_hline(yintercept = 0, size = 0.5) + 
    ggplot2::geom_vline(xintercept = 0, size = 0.5) + 
    ggplot2::xlab(expression("Predicted effect (EVE pathogenicity score)")) +
    ggplot2::ylab(expression("Inferred Folding "*Delta*Delta*"G (ddPCA)")) +
    ggplot2::facet_wrap(protein~., scales = "free") +
    ggplot2::geom_text(data = cor_dt, ggplot2::aes(label=paste("R = ", cor, sep="")), hjust = 1, vjust = 1) +
    ggplot2::theme_bw()
  ggplot2::ggsave(file.path(outpath, "eve_vs_b_ddg_position_class_all.pdf"), d, width = 8, height = 3, useDingbats=FALSE)

}

