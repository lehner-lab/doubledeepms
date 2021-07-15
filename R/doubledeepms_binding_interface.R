
#' doubledeepms_binding_interface
#'
#' Plot plots related to binding free energy and protein binding interface
#'
#' @param input_list path to MoCHI thermo model fit results (required)
#' @param outpath output path for plots and saved objects (required)
#' @param colour_scheme colour scheme file (required)
#' @param execute whether or not to execute the analysis (default: TRUE)
#'
#' @return Nothing
#' @export
#' @import data.table
doubledeepms_binding_interface <- function(
  input_list,
  outpath,
  colour_scheme,
  execute = TRUE
){
  
  #Return if analysis not executed
  if(!execute){
    return()
  }
  
  #Display status
  message(paste("\n\n*******", "running stage: doubledeepms_binding_interface_plots", "*******\n\n"))
  
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
  ### Position class violin plots
  ###########################
  
  plot_dt <- dg_dt[b_ddg_pred_conf==T & id!="-0-",.(protein, b_ddg_pred, Pos_class, RSASA, scHAmin_ligand, Pos_ref)]
  plot_dt[, Pos_class_plot := stringr::str_to_title(Pos_class)]
  plot_dt[Pos_class=="binding_interface", Pos_class_plot := "Binding\ninterface"]
  
  for(domain in names(input_list)){
    temp_test <- t.test(
      plot_dt[Pos_class=="binding_interface" & protein == domain, b_ddg_pred],
      plot_dt[!Pos_class=="binding_interface" & protein == domain, b_ddg_pred])
    print(paste0("Binding free energy change of mutations in binding interface higher than remainder (", i, "): p-value=", format(temp_test[["p.value"]], digits=2, scientific=T)))
  }
  
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(Pos_class_plot, b_ddg_pred, fill = Pos_class_plot)) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
    ggplot2::xlab("Residue position") +
    ggplot2::ylab(expression(Delta*Delta*"G of Binding")) +
    ggplot2::facet_grid(~protein, scales = "free") + 
    ggplot2::theme_classic() +
    ggplot2::coord_cartesian(ylim = c(-3, 8)) +
    ggplot2::labs(fill = "Residue\nposition")
  if(!is.null(colour_scheme)){
    d <- d + ggplot2::scale_fill_manual(values = unlist(colour_scheme[["shade 0"]][c(1, 3, 4)]))
  }
  suppressWarnings(ggplot2::ggsave(file.path(outpath, "position_violins_all.pdf"), d, width = 7, height = 3, useDingbats=FALSE))
  
  plot_dt <- dg_dt[b_ddg_pred_conf==T & id!="-0-" & protein!="GB1",.(protein, b_ddg_pred, Pos_class, RSASA, scHAmin_ligand)]
  plot_dt[, Pos_class_plot := stringr::str_to_title(Pos_class)]
  plot_dt[Pos_class=="binding_interface", Pos_class_plot := "Binding\ninterface"]
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(Pos_class_plot, b_ddg_pred, fill = Pos_class_plot)) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
    ggplot2::xlab("") +
    ggplot2::ylab(expression(Delta*Delta*"G of Binding")) +
    ggplot2::facet_grid(~protein, scales = "free") + 
    ggplot2::theme_classic() +
    ggplot2::coord_cartesian(ylim = c(-3, 8)) +
    ggplot2::labs(fill = "Residue\nposition") 
  if(!is.null(colour_scheme)){
    d <- d + ggplot2::scale_fill_manual(values = unlist(colour_scheme[["shade 0"]][c(1, 3, 4)]))
  }
  suppressWarnings(ggplot2::ggsave(file.path(outpath, "position_violins.pdf"), d, width = 5, height = 3, useDingbats=FALSE))
  
  
}

