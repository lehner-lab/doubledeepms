
#' doubledeepms_fitness_heatmaps
#'
#' Plot fitness heatmaps.
#'
#' @param input_file path to MoCHI thermo model fit results (required)
#' @param input_file_fitness path to fitness data (required)
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
doubledeepms_fitness_heatmaps <- function(
  input_file,
  input_file_fitness,
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
  message(paste("\n\n*******", "running stage: doubledeepms_fitness_heatmaps", "*******\n\n"))

  #Create output directory
  doubledeepms__create_dir(doubledeepms_dir = outpath)

  #Load fitness and metrics data
  fit_list <- list()
  #Load free energies
  metrics_dt <- unique(fread(input_file)[id!="-0-",.(Pos, Pos_class, scHAmin_ligand, Pos_ref)])
  for(pca_type in c("Abundance", "Binding")){
    rdata_file <- list.files(file.path(input_file_fitness, pca_type))
    load(file.path(input_file_fitness, pca_type, rdata_file))
    all_variants[Nham_aa<=1 & STOP==F & STOP_readthrough==F, pca_type := pca_type]
    #Single mutant position
    wt_seq <- unlist(strsplit(all_variants[WT==T,aa_seq], ""))
    #Single mutant position class
    all_variants[Nham_aa==1, Pos := which(unlist(strsplit(aa_seq, "")) != wt_seq), aa_seq]
    all_variants[Nham_aa==1, WT_AA := wt_seq[Pos], aa_seq]
    all_variants[Nham_aa==1, Mut := substr(aa_seq, Pos, Pos), aa_seq]
    all_variants <- merge(all_variants, metrics_dt, by = "Pos", all = T)
    fit_list[[pca_type]] <- all_variants
  }
  fitness_dt <- rbindlist(fit_list, fill = T)

  ###########################
  ### Folding heatmap
  ###########################

  doubledeepms__plot_heatmap(
    input_dt = fitness_dt[pca_type=="Abundance" & Nham_aa==1,.(fitness = -fitness, fitness_conf = T, Pos_ref, scHAmin_ligand, WT_AA, Mut, Pos_class)],
    variable_name = "fitness",
    output_file = file.path(outpath, "abundance_heatmap.pdf"),
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
    input_dt = fitness_dt[pca_type=="Binding" & Nham_aa==1,.(fitness = -fitness, fitness_conf = T, Pos_ref, scHAmin_ligand, WT_AA, Mut, Pos_class)],
    variable_name = "fitness",
    output_file = file.path(outpath, "binding_heatmap.pdf"),
    width = plot_width,
    height = plot_height,
    plot_title = paste0(domain_name, " amino acid position"),
    colour_clip = 2.5,
    colour_low = colour_scheme[["shade 0"]][[3]],
    colour_high = colour_scheme[["shade 0"]][[1]])

}

