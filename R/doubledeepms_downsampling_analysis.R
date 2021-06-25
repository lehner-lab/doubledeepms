
#' doubledeepms_downsampling_analysis
#'
#' Downsampling analysis.
#'
#' @param mochi_outpath_prefix path prefix to MoCHI thermo model fit results (required)
#' @param sample_percentage sample percentages (required)
#' @param outpath_suffix output suffix (default:p)
#' @param literature_comparison_prefix literature comparison path prefix (default:empty string)
#' @param literature_comparison_list literature comparison file list (default:empty list)
#' @param temperature temperature in degrees celcuis (default:24)
#' @param outpath output path for plots and saved objects (required)
#' @param colour_scheme colour scheme file (required)
#' @param execute whether or not to execute the analysis (default: TRUE)
#'
#' @return Nothing
#' @export
#' @import data.table
doubledeepms_downsampling_analysis <- function(
  mochi_outpath_prefix,
  sample_percentage,
  outpath_suffix="p",
  literature_comparison_prefix="",
  literature_comparison_list=list(),
  temperature = 24,
  outpath,
  colour_scheme,
  execute = TRUE
  ){

  #Return if analysis not executed
  if(!execute){
    return()
  }

  #Domain name
  domain_name <- rev(unlist(strsplit(basename(outpath), "_")))[1]

  #Display status
  message(paste("\n\n*******", paste("running stage: doubledeepms_downsampling_analysis for", domain_name), "*******\n\n"))

  #Create output directory
  doubledeepms__create_dir(doubledeepms_dir = outpath)

  #Constants
  gas_constant <- 0.001987
  RT <- gas_constant*(273+temperature)

  #Mochi output list
  mochi_outpath_list <- as.list(gsub("_subsample100p", "", paste0(mochi_outpath_prefix, sample_percentage, outpath_suffix)))
  names(mochi_outpath_list) <- as.character(sample_percentage)

  #Literature comparison list
  if(length(literature_comparison_list)!=0){
    literature_comparison_list_folding <- as.list(gsub("subsample100p", "", file.path(paste0(literature_comparison_prefix, sample_percentage, outpath_suffix), literature_comparison_list[["folding"]])))
    names(literature_comparison_list_folding) <- as.character(sample_percentage)
    literature_comparison_list_binding <- as.list(gsub("subsample100p", "", file.path(paste0(literature_comparison_prefix, sample_percentage, outpath_suffix), literature_comparison_list[["binding"]])))
    names(literature_comparison_list_binding) <- as.character(sample_percentage)
    literature_comparison_list_all <- list(
      "folding" = literature_comparison_list_folding,
      "binding" = literature_comparison_list_binding)
  }else{
    literature_comparison_list_all <- list()
  }

  #Load model data
  fitness_dt <- fread(file.path(mochi_outpath_list[["100"]], "model_data.txt"))

  #Load model results
  pred_dt <- doubledeepms__get_model_results(
    input_folder = mochi_outpath_list[["100"]], 
    input_dt = fitness_dt, 
    RT = RT)[mut_order==2]

  #Mutation ids
  pred_dt[, mut1 := sapply(strsplit(id, ","), '[', 1)]
  pred_dt[, mut2 := sapply(strsplit(id, ","), '[', 2)]

  #Add predicted fitness for all subsamples
  for(i in names(mochi_outpath_list)){
    if(i=="100"){
      next
    }

    #Load coefficients
    coef_dt <- fread(file.path(mochi_outpath_list[[i]], "model_weights_0.txt"))

    #Training set ids
    fold_training_set_ids <- fread(file.path(mochi_outpath_list[[i]], "model_data.txt"))[training_set==1 & dataset_binding==0,id]
    bind_training_set_ids <- fread(file.path(mochi_outpath_list[[i]], "model_data.txt"))[training_set==1 & dataset_binding==1,id]

    #Dictionaries
    folding_coef_list <- coef_dt[,folding_coefficient]
    names(folding_coef_list) <- coef_dt[,id]
    binding_coef_list <- coef_dt[,binding_coefficient]
    names(binding_coef_list) <- coef_dt[,id]

    #Energies
    pred_dt[, f_dg_pred_ss := apply(data.table(
      dg_wt = folding_coef_list[["WT"]], 
      dg_mut1 = folding_coef_list[mut1], 
      dg_mut2 = folding_coef_list[mut2]), 1, sum, na.rm = T)*RT]
    pred_dt[, b_dg_pred_ss := apply(data.table(
      dg_wt = binding_coef_list[["WT"]], 
      dg_mut1 = binding_coef_list[mut1], 
      dg_mut2 = binding_coef_list[mut2]), 1, sum, na.rm = T)*RT]

    #Fitness - folding
    pred_dt[dataset_binding==0, paste0("predicted_fitness_", i) := doubledeepms__predict_fitness(
      mochi_outpath = mochi_outpath_list[[i]], 
      folding_energy = f_dg_pred_ss,
      binding_energy = b_dg_pred_ss,
      RT = RT)[["fitness_folding"]]]
    pred_dt[dataset_binding==0 & id %in% fold_training_set_ids, paste0("training_set_", i) := 1]
    pred_dt[dataset_binding==0 & !id %in% fold_training_set_ids, paste0("training_set_", i) := 0]
    #Fitness - binding
    pred_dt[dataset_binding==1, paste0("predicted_fitness_", i) := doubledeepms__predict_fitness(
      mochi_outpath = mochi_outpath_list[[i]], 
      folding_energy = f_dg_pred_ss,
      binding_energy = b_dg_pred_ss,
      RT = RT)[["fitness_binding"]]]
    pred_dt[dataset_binding==1 & id %in% bind_training_set_ids, paste0("training_set_", i) := 1]
    pred_dt[dataset_binding==1 & !id %in% bind_training_set_ids, paste0("training_set_", i) := 0]
  }

  #Plot - subsample correlation with observed fitness
  r2_list <- list()
  for(i in names(mochi_outpath_list)){
    if(i=="100"){
      cname <- "predicted_fitness"
    }else{
      cname <- paste0("predicted_fitness_", i)
    }
    r2_list[[paste0("r2_", i, "_abundance_train")]] <- pred_dt[dataset_binding==0 & get(gsub("predicted_fitness", "training_set", cname))==1,cor(.SD[[1]], .SD[[2]])^2,,.SDcols = c("observed_fitness", cname)]
    r2_list[[paste0("r2_", i, "_abundance_val")]] <- pred_dt[dataset_binding==0 & get(gsub("predicted_fitness", "training_set", cname))==0,cor(.SD[[1]], .SD[[2]])^2,,.SDcols = c("observed_fitness", cname)]
    r2_list[[paste0("r2_", i, "_binding_train")]] <- pred_dt[dataset_binding==1 & get(gsub("predicted_fitness", "training_set", cname))==1,cor(.SD[[1]], .SD[[2]])^2,,.SDcols = c("observed_fitness", cname)]
    r2_list[[paste0("r2_", i, "_binding_val")]] <- pred_dt[dataset_binding==1 & get(gsub("predicted_fitness", "training_set", cname))==0,cor(.SD[[1]], .SD[[2]])^2,,.SDcols = c("observed_fitness", cname)]
  }
  plot_dt <- data.table(r2 = unlist(r2_list), name = names(r2_list))
  plot_dt[, phenotype := sapply(strsplit(name, "_"), '[', 3)]
  plot_dt[, training_set := sapply(lapply(strsplit(name, "_"), rev), '[', 1)=="train"]
  plot_dt[, sample_size := sapply(strsplit(name, "_"), '[', 2)]
  plot_dt[, sample_size_num := as.numeric(sapply(strsplit(name, "_"), '[', 2))]
  plot_dt[, sample_size_plot := factor(sample_size, levels = unique(sample_size))]
  plot_cols <- unlist(colour_scheme[["shade 0"]][c(1,3)])
  names(plot_cols) <- c("binding", "abundance")
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(sample_size_num, r2, color = phenotype, shape = training_set)) +
    ggplot2::geom_point() +
    ggplot2::geom_line(ggplot2::aes(linetype = training_set)) +
    ggplot2::xlab("%Retained double AA mutants") +
    ggplot2::ylab("%Fitness variance explained") +
    ggplot2::scale_colour_manual(values=plot_cols) +
    ggplot2::theme_classic()
  ggplot2::ggsave(file.path(outpath, "subsample_performance_folding_binding.pdf"), d, width = 4, height = 3, useDingbats=FALSE)

  #Plot - subsample correlation with validations
  if(length(literature_comparison_list_all)!=0){
    cor_list <- list()
    for(i in names(literature_comparison_list_all[["folding"]])){
      val_dt <- fread(literature_comparison_list_all[["folding"]][[i]])
      cor_list[[paste0("cor_", i, "_folding")]] <- val_dt[,cor(col_lit, col_data)]
    }
    for(i in names(literature_comparison_list_all[["binding"]])){
      val_dt <- fread(literature_comparison_list_all[["binding"]][[i]])
      cor_list[[paste0("cor_", i, "_binding")]] <- val_dt[,cor(col_lit, col_data)]
    }
    plot_dt <- data.table(cor = unlist(cor_list), name = names(cor_list))
    plot_dt[, phenotype := sapply(strsplit(name, "_"), '[', 3)]
    plot_dt[, sample_size := sapply(strsplit(name, "_"), '[', 2)]
    plot_dt[, sample_size_num := as.numeric(sapply(strsplit(name, "_"), '[', 2))]
    plot_dt[, sample_size_plot := factor(sample_size, levels = unique(sample_size))]
    plot_cols <- unlist(colour_scheme[["shade 0"]][c(1,3)])
    names(plot_cols) <- c("binding", "folding")
    d <- ggplot2::ggplot(plot_dt,ggplot2::aes(sample_size_num, cor, color = phenotype)) +
      ggplot2::geom_point() +
      ggplot2::geom_line() +
      ggplot2::xlab("%Retained double AA mutants") +
      ggplot2::ylab(expression("Correlation with literature "*Delta*Delta*"Gs")) +
      ggplot2::scale_colour_manual(values=plot_cols) +
      ggplot2::theme_classic()
    ggplot2::ggsave(file.path(outpath, "subsample_validation_cor_folding_binding.pdf"), d, width = 4, height = 3, useDingbats=FALSE)
  }

}

