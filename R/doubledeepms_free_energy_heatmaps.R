
#' doubledeepms_free_energy_heatmaps
#'
#' Plot free energy heatmaps.
#'
#' @param input_file path to MoCHI thermo model fit results (required)
#' @param domain_name domain name (required)
#' @param outpath output path for plots and saved objects (required)
#' @param colour_scheme colour scheme file (required)
#' @param plot_width heatmap plot width in inches (default:10)
#' @param plot_height heatmap plot height in inches (default:4)
#' @param execute whether or not to execute the analysis (default: TRUE)
#'
#' @return Nothing
#' @export
#' @import data.table
doubledeepms_free_energy_heatmaps <- function(
  input_file,
  domain_name,
  outpath,
  colour_scheme,
  plot_width = 10,
  plot_height = 4,
  execute = TRUE
  ){

  #Return if analysis not executed
  if(!execute){
    return()
  }

  #Display status
  message(paste("\n\n*******", paste("running stage: doubledeepms_free_energy_heatmaps for", domain_name), "*******\n\n"))

  #Create output directory
  doubledeepms__create_dir(doubledeepms_dir = outpath)

  #Load free energies
  dg_dt <- fread(input_file)[id!="-0-"]

  #Add WT and mutant AAs
  dg_dt[, WT_AA := substr(id, 1, 1)]
  dg_dt[, Mut := substr(id, nchar(id), nchar(id))]

  #Per residue metrics
  for(i in c("f_ddg", "b_ddg")){
    dg_dt[,paste0(i, "_posmeanabs") := mean(abs(.SD[[1]]), na.rm = T),Pos_ref,.SDcols = paste0(i, c("_pred"))]
    dg_dt[,paste0(i, "_posse") := sd(abs(.SD[[1]]), na.rm = T)/sqrt(sum(!is.na(.SD[[1]]))),Pos_ref,.SDcols = paste0(i, c("_pred"))]
    dg_dt[,paste0(i, "_wposmeanabs") := sum(abs(.SD[[1]])/.SD[[2]]^2, na.rm = T)/sum(1/.SD[[2]]^2, na.rm = T),Pos_ref,.SDcols = paste0(i, c("_pred", "_pred_sd"))]
    dg_dt[,paste0(i, "_wposse") := sqrt(1/sum(1/.SD[[2]]^2, na.rm = T)),Pos_ref,.SDcols = paste0(i, c("_pred", "_pred_sd"))]
  }

  #Per residue metrics - confident only
  for(i in c("f_ddg", "b_ddg")){
    dg_dt[get(paste0(i, "_pred_conf"))==TRUE,paste0(i, "_pred_filtered") := .SD[[1]],,.SDcols = paste0(i, "_pred")]
    dg_dt[get(paste0(i, "_pred_conf"))==TRUE,paste0(i, "_pred_sd_filtered") := .SD[[1]],,.SDcols = paste0(i, "_pred_sd")]
    dg_dt[,paste0(i, "_posmeanabs_conf") := mean(abs(.SD[[1]]), na.rm = T),Pos_ref,.SDcols = paste0(i, c("_pred_filtered"))]
    dg_dt[,paste0(i, "_posse_conf") := sd(abs(.SD[[1]]), na.rm = T)/sqrt(sum(!is.na(.SD[[1]]))),Pos_ref,.SDcols = paste0(i, c("_pred_filtered"))]
    dg_dt[,paste0(i, "_wposmeanabs_conf") := sum(abs(.SD[[1]])/.SD[[2]]^2, na.rm = T)/sum(1/.SD[[2]]^2, na.rm = T),Pos_ref,.SDcols = paste0(i, c("_pred_filtered", "_pred_sd_filtered"))]
    dg_dt[,paste0(i, "_wposse_conf") := sqrt(1/sum(1/.SD[[2]]^2, na.rm = T)),Pos_ref,.SDcols = paste0(i, c("_pred_filtered", "_pred_sd_filtered"))]
  }

  ###########################
  ### Folding heatmap
  ###########################

  doubledeepms__plot_heatmap(
    input_dt = dg_dt,
    variable_name = "f_ddg_pred",
    metric_name = "f_ddg_posmeanabs_conf",
    output_file = file.path(outpath, "folding_heatmap.pdf"),
    width = plot_width,
    height = plot_height,
    plot_title = paste0(domain_name, " amino acid position"),
    colour_clip = 2.5,
    colour_low = colour_scheme[["shade 0"]][[3]],
    colour_high = colour_scheme[["shade 0"]][[1]])

  ###########################
  ### Binding heatmap
  ###########################

  doubledeepms__plot_heatmap(
    input_dt = dg_dt,
    variable_name = "b_ddg_pred",
    metric_name = "b_ddg_posmeanabs_conf",
    output_file = file.path(outpath, "binding_heatmap.pdf"),
    width = plot_width,
    height = plot_height,
    plot_title = paste0(domain_name, " amino acid position"),
    colour_clip = 2.5,
    colour_low = colour_scheme[["shade 0"]][[3]],
    colour_high = colour_scheme[["shade 0"]][[1]])

}

