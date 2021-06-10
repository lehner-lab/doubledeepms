
#' doubledeepms_allostery_scatterplots
#'
#' Plot free energy scatterplots for allosteric mutations
#'
#' @param input_file path to input file (required)
#' @param temperature temperature in degrees celcuis (default:24)
#' @param fitness_list list of folder paths to fitness data (required)
#' @param mochi_outpath_list list of paths to MoCHI thermo model fit results (required)
#' @param outpath output path for plots and saved objects (required)
#' @param colour_scheme colour scheme file (required)
#' @param execute whether or not to execute the analysis (default: TRUE)
#'
#' @return Nothing
#' @export
#' @import data.table
doubledeepms_allostery_scatterplots <- function(
  input_file,
  temperature = 24,
  fitness_list,
  mochi_outpath_list,
  outpath,
  colour_scheme,
  execute = TRUE
  ){

  #Return if analysis not executed
  if(!execute){
    return()
  }

  #Display status
  message(paste("\n\n*******", "running stage: doubledeepms_allostery_scatterplots", "*******\n\n"))

  #Create output directory
  doubledeepms__create_dir(doubledeepms_dir = outpath)

  ### Load single mutant free energies
  ###########################

  #Constants
  gas_constant <- 0.001987
  RT <- gas_constant*(273+temperature)

  #Load dg data
  dg_dt <- fread(input_file)

  #Load fitness and metrics data
  fit_list <- list()
  for(protein in names(fitness_list)){
    fit_list[[protein]] <- list()
    for(pca_type in c("Abundance", "Binding")){
      rdata_file <- list.files(file.path(fitness_list[[protein]], pca_type))
      load(file.path(fitness_list[[protein]], pca_type, rdata_file))
      all_variants[, protein := protein]
      all_variants[, pca_type := pca_type]
      #Single mutant position
      wt_seq <- unlist(strsplit(all_variants[WT==T,aa_seq], ""))
      #Single mutant position class
      all_variants[Nham_aa==1, Pos := which(unlist(strsplit(aa_seq, "")) != wt_seq), aa_seq]
      all_variants[Nham_aa==1, WT_AA := wt_seq[Pos], aa_seq]
      all_variants[Nham_aa==1, Mut := substr(aa_seq, Pos, Pos), aa_seq]
      all_variants[Nham_aa==1, id := paste0(WT_AA, Pos, Mut), aa_seq]
      all_variants[Nham_aa==0, id := "-0-"]
      fit_list[[protein]][[pca_type]] <- all_variants
    }
  }
  fitness_dt <- rbindlist(unlist(fit_list, recursive = FALSE), fill = T)
  #Add allostery information
  fitness_dt <- merge(
    fitness_dt[pca_type=="Abundance" & (Nham_aa==1 | WT==T) & STOP==F & STOP_readthrough==F,.(id, fitness_abundance = fitness, protein, sigma_abundance = sigma)], 
    fitness_dt[pca_type=="Binding" & (Nham_aa==1 | WT==T) & STOP==F & STOP_readthrough==F,.(id, fitness_binding = fitness, protein, sigma_binding = sigma)], by = c("id", "protein"))
  fitness_dt <- merge(
    fitness_dt, 
    dg_dt[!duplicated(id),.(id, allosteric, allosteric_mutation, Pos_ref, Pos_class)], by = "id", all.x = T)

  #Add zero binding model fitness
  for(i in fitness_dt[,unique(protein)]){
    #Load model
    modpar <- fread(file.path(mochi_outpath_list[[i]], "model_parameters_0.txt"), header = F)[,V1]
    modpar_list <- as.list(as.numeric(modpar[seq(2, length(modpar), 2)]))
    names(modpar_list) <- modpar[seq(1, length(modpar), 2)]
    #Zero binding ddG model mapping Binding fitness to folding fitness
    fitness_dt[protein==i, fitness_binding_zbm := doubledeepms__predict_binding_fitness_from_folding_fitness(
      folding_fitness = unlist(fitness_abundance),
      folding_linear_kernel = modpar_list[["folding_linear_kernel"]],
      folding_linear_bias = modpar_list[["folding_linear_bias"]],
      binding_linear_kernel = modpar_list[["binding_linear_kernel"]],
      binding_linear_bias = modpar_list[["binding_linear_bias"]],
      b_dg_wt = dg_dt[id=="-0-" & protein==i,b_dg_pred][1],
      RT = RT)]
  }

  ###########################
  ### Free energy scatterplots for allosteric sites
  ###########################

  #Free energy scatterplots by protein - all - conf
  for(i in dg_dt[,unique(protein)]){
    plot_dt <- copy(dg_dt)[protein==i][,.(f_dg_pred, f_ddg_pred_sd, b_dg_pred, b_ddg_pred_sd, f_ddg_pred_conf, b_ddg_pred_conf, Pos_class, allosteric, id, Pos_ref)]
    plot_dt <- plot_dt[f_ddg_pred_conf==T & b_ddg_pred_conf==T,]
    plot_dt[, Pos_ref_plot := factor(Pos_ref)]
    #Plot
    d <- ggplot2::ggplot(plot_dt[id!="-0-"],ggplot2::aes(f_dg_pred, b_dg_pred)) +
      ggplot2::geom_point(alpha = 0.5, size = 1, color = "lightgrey") +
      ggplot2::geom_point(data = plot_dt[allosteric==T & Pos_class !="binding_interface",], alpha = 0.5, ggplot2::aes(color = Pos_ref_plot), size = 2) +
      ggplot2::geom_linerange(data = plot_dt[allosteric==T & Pos_class !="binding_interface",], alpha = 0.25, ggplot2::aes(ymin = b_dg_pred-b_ddg_pred_sd*1.96, ymax = b_dg_pred+b_ddg_pred_sd*1.96, color = Pos_ref_plot)) +
      ggplot2::geom_linerange(data = plot_dt[allosteric==T & Pos_class !="binding_interface",], alpha = 0.25, ggplot2::aes(xmin = f_dg_pred-f_ddg_pred_sd*1.96, xmax = f_dg_pred+f_ddg_pred_sd*1.96, color = Pos_ref_plot)) +
      ggplot2::geom_vline(xintercept = 0) +
      ggplot2::geom_vline(data = plot_dt[id=="-0-",], ggplot2::aes(xintercept = f_dg_pred), linetype = 2) +
      ggplot2::geom_hline(yintercept = 0) +
      ggplot2::geom_hline(data = plot_dt[id=="-0-",], ggplot2::aes(yintercept = b_dg_pred), linetype = 2) +
      ggplot2::xlab(expression("Folding "*Delta*"G")) +
      ggplot2::ylab(expression("Binding "*Delta*"G")) +
      ggplot2::labs(color = "Allosteric\nsite") +
      ggplot2::theme_classic()
    ggplot2::ggsave(file.path(outpath, paste0("dG_scatter_allosteric_sites_", i, ".pdf")), d, width = 4, height = 3, useDingbats=FALSE)
    #Plot - xylim
    d <- ggplot2::ggplot(plot_dt[id!="-0-"],ggplot2::aes(f_dg_pred, b_dg_pred)) +
      ggplot2::geom_point(alpha = 0.5, size = 1, color = "lightgrey") +
      ggplot2::geom_point(data = plot_dt[allosteric==T & Pos_class !="binding_interface",], alpha = 0.5, ggplot2::aes(color = Pos_ref_plot), size = 2) +
      ggplot2::geom_linerange(data = plot_dt[allosteric==T & Pos_class !="binding_interface",], alpha = 0.25, ggplot2::aes(ymin = b_dg_pred-b_ddg_pred_sd*1.96, ymax = b_dg_pred+b_ddg_pred_sd*1.96, color = Pos_ref_plot)) +
      ggplot2::geom_linerange(data = plot_dt[allosteric==T & Pos_class !="binding_interface",], alpha = 0.25, ggplot2::aes(xmin = f_dg_pred-f_ddg_pred_sd*1.96, xmax = f_dg_pred+f_ddg_pred_sd*1.96, color = Pos_ref_plot)) +
      ggplot2::geom_vline(xintercept = 0) +
      ggplot2::geom_vline(data = plot_dt[id=="-0-",], ggplot2::aes(xintercept = f_dg_pred), linetype = 2) +
      ggplot2::geom_hline(yintercept = 0) +
      ggplot2::geom_hline(data = plot_dt[id=="-0-",], ggplot2::aes(yintercept = b_dg_pred), linetype = 2) +
      ggplot2::xlab(expression("Folding "*Delta*"G")) +
      ggplot2::ylab(expression("Binding "*Delta*"G")) +
      ggplot2::labs(color = "Allosteric\nsite") +
      ggplot2::theme_classic()
    if(i!="GB1"){
      d <- d + ggplot2::coord_cartesian(ylim = c(-2.5, 2.5), xlim = c(-2.5,2.5))
    }
    ggplot2::ggsave(file.path(outpath, paste0("dG_scatter_allosteric_sites_", i, "_xylim.pdf")), d, width = 4, height = 3, useDingbats=FALSE)
  }

  ###########################
  ### Free energy scatterplots for allosteric mutations
  ###########################

  #Free energy scatterplots by protein - all - conf
  for(i in dg_dt[,unique(protein)]){
    plot_dt <- copy(dg_dt)[protein==i][,.(f_dg_pred, f_ddg_pred_sd, b_dg_pred, b_ddg_pred_sd, f_ddg_pred_conf, b_ddg_pred_conf, Pos_class, allosteric, allosteric_mutation, id, Pos_ref)]
    plot_dt <- plot_dt[f_ddg_pred_conf==T & b_ddg_pred_conf==T,]
    plot_dt[, Pos_ref_plot := "Remainder"]
    plot_dt[allosteric==T & Pos_class != "binding_interface", Pos_ref_plot := "Site"]
    plot_dt[allosteric_mutation==T & Pos_class != "binding_interface", Pos_ref_plot := "Mutation"]
    plot_cols = c("grey", colour_scheme[["shade 0"]][c(2,4)])
    names(plot_cols) <- c("Remainder", "Site", "Mutation")
    #Plot
    d <- ggplot2::ggplot(plot_dt[id!="-0-"],ggplot2::aes(f_dg_pred, b_dg_pred)) +
      ggplot2::geom_point(data = plot_dt[Pos_ref_plot=="Remainder"], alpha = 0.5, size = 1, color = "lightgrey") +
      ggplot2::geom_point(data = plot_dt[Pos_ref_plot!="Remainder"], alpha = 0.5, ggplot2::aes(color = Pos_ref_plot), size = 2) +
      ggplot2::geom_linerange(data = plot_dt[Pos_ref_plot!="Remainder"], alpha = 0.25, ggplot2::aes(ymin = b_dg_pred-b_ddg_pred_sd*1.96, ymax = b_dg_pred+b_ddg_pred_sd*1.96, color = Pos_ref_plot)) +
      ggplot2::geom_linerange(data = plot_dt[Pos_ref_plot!="Remainder"], alpha = 0.25, ggplot2::aes(xmin = f_dg_pred-f_ddg_pred_sd*1.96, xmax = f_dg_pred+f_ddg_pred_sd*1.96, color = Pos_ref_plot)) +
      ggplot2::geom_vline(xintercept = 0) +
      ggplot2::geom_vline(data = plot_dt[id=="-0-",], ggplot2::aes(xintercept = f_dg_pred), linetype = 2) +
      ggplot2::geom_hline(yintercept = 0) +
      ggplot2::geom_hline(data = plot_dt[id=="-0-",], ggplot2::aes(yintercept = b_dg_pred), linetype = 2) +
      ggplot2::xlab(expression("Folding "*Delta*"G")) +
      ggplot2::ylab(expression("Binding "*Delta*"G")) +
      ggplot2::labs(color = "Allosteric") +
      ggplot2::scale_colour_manual(values=plot_cols) +
      ggplot2::theme_classic()
    ggplot2::ggsave(file.path(outpath, paste0("dG_scatter_allosteric_mutations_", i, ".pdf")), d, width = 4, height = 3, useDingbats=FALSE)
    #Plot - xylim
    d <- ggplot2::ggplot(plot_dt[id!="-0-"],ggplot2::aes(f_dg_pred, b_dg_pred)) +
      ggplot2::geom_point(data = plot_dt[Pos_ref_plot=="Remainder"], alpha = 0.5, size = 1, color = "lightgrey") +
      ggplot2::geom_point(data = plot_dt[Pos_ref_plot!="Remainder"], alpha = 0.5, ggplot2::aes(color = Pos_ref_plot), size = 2) +
      ggplot2::geom_linerange(data = plot_dt[Pos_ref_plot!="Remainder"], alpha = 0.25, ggplot2::aes(ymin = b_dg_pred-b_ddg_pred_sd*1.96, ymax = b_dg_pred+b_ddg_pred_sd*1.96, color = Pos_ref_plot)) +
      ggplot2::geom_linerange(data = plot_dt[Pos_ref_plot!="Remainder"], alpha = 0.25, ggplot2::aes(xmin = f_dg_pred-f_ddg_pred_sd*1.96, xmax = f_dg_pred+f_ddg_pred_sd*1.96, color = Pos_ref_plot)) +
      ggplot2::geom_vline(xintercept = 0) +
      ggplot2::geom_vline(data = plot_dt[id=="-0-",], ggplot2::aes(xintercept = f_dg_pred), linetype = 2) +
      ggplot2::geom_hline(yintercept = 0) +
      ggplot2::geom_hline(data = plot_dt[id=="-0-",], ggplot2::aes(yintercept = b_dg_pred), linetype = 2) +
      ggplot2::xlab(expression("Folding "*Delta*"G")) +
      ggplot2::ylab(expression("Binding "*Delta*"G")) +
      ggplot2::labs(color = "Allosteric") +
      ggplot2::scale_colour_manual(values=plot_cols) +
      ggplot2::theme_classic()
    if(i!="GB1"){
      d <- d + ggplot2::coord_cartesian(ylim = c(-2.5, 2.5), xlim = c(-2.5,2.5))
    }
    ggplot2::ggsave(file.path(outpath, paste0("dG_scatter_allosteric_mutations_", i, "_xylim.pdf")), d, width = 4, height = 3, useDingbats=FALSE)
  }

  ###########################
  ### Fitness scatterplots for allosteric sites
  ###########################

  #Fitness scatterplots by protein - all
  for(i in fitness_dt[,unique(protein)]){
    plot_dt <- copy(fitness_dt)[protein==i][,.(fitness_abundance, sigma_abundance, fitness_binding, fitness_binding_zbm, sigma_binding, Pos_class, allosteric, id, Pos_ref)]
    plot_dt[, Pos_ref_plot := factor(Pos_ref)]
    #Plot
    d <- ggplot2::ggplot(plot_dt[id!="-0-"],ggplot2::aes(fitness_abundance, fitness_binding)) +
      ggplot2::geom_point(alpha = 0.5, size = 1, color = "lightgrey") +
      ggplot2::geom_point(data = plot_dt[allosteric==T & Pos_class !="binding_interface",], alpha = 0.5, ggplot2::aes(color = Pos_ref_plot), size = 2) +
      ggplot2::geom_linerange(data = plot_dt[allosteric==T & Pos_class !="binding_interface",], alpha = 0.25, ggplot2::aes(ymin = fitness_binding-sigma_binding*1.96, ymax = fitness_binding+sigma_binding*1.96, color = Pos_ref_plot)) +
      ggplot2::geom_linerange(data = plot_dt[allosteric==T & Pos_class !="binding_interface",], alpha = 0.25, ggplot2::aes(xmin = fitness_abundance-sigma_abundance*1.96, xmax = fitness_abundance+sigma_abundance*1.96, color = Pos_ref_plot)) +
      ggplot2::geom_line(data = plot_dt[fitness_binding_zbm>=min(fitness_binding)], ggplot2::aes(fitness_abundance, fitness_binding_zbm), color = colour_scheme[["shade 0"]][[1]]) +
      ggplot2::geom_vline(xintercept = 0) +
      ggplot2::geom_abline(linetype = 2) +
      ggplot2::geom_hline(yintercept = 0) +
      ggplot2::xlab("Fitness (Abundance)") +
      ggplot2::ylab("Fitness (Binding)") +
      ggplot2::labs(color = "Allosteric\nsite") +
      ggplot2::theme_classic()
    ggplot2::ggsave(file.path(outpath, paste0("fitness_scatter_allosteric_sites_", i, ".pdf")), d, width = 4, height = 3, useDingbats=FALSE)
  }

  ###########################
  ### Fitness scatterplots for allosteric mutations
  ###########################

  #Fitness scatterplots by protein - all
  for(i in fitness_dt[,unique(protein)]){
    plot_dt <- copy(fitness_dt)[protein==i][,.(fitness_abundance, sigma_abundance, fitness_binding, fitness_binding_zbm, sigma_binding, Pos_class, allosteric, allosteric_mutation, id, Pos_ref)]
    plot_dt[, Pos_ref_plot := "Remainder"]
    plot_dt[allosteric==T & Pos_class != "binding_interface", Pos_ref_plot := "Site"]
    plot_dt[allosteric_mutation==T & Pos_class != "binding_interface", Pos_ref_plot := "Mutation"]
    plot_cols = c("grey", colour_scheme[["shade 0"]][c(2,4)])
    names(plot_cols) <- c("Remainder", "Site", "Mutation")
    #Plot
    d <- ggplot2::ggplot(plot_dt[id!="-0-"],ggplot2::aes(fitness_abundance, fitness_binding)) +
      ggplot2::geom_point(data = plot_dt[Pos_ref_plot=="Remainder"], alpha = 0.5, size = 1, color = "lightgrey") +
      ggplot2::geom_point(data = plot_dt[Pos_ref_plot!="Remainder"], alpha = 0.5, ggplot2::aes(color = Pos_ref_plot), size = 2) +
      ggplot2::geom_linerange(data = plot_dt[Pos_ref_plot!="Remainder"], alpha = 0.25, ggplot2::aes(ymin = fitness_binding-sigma_binding*1.96, ymax = fitness_binding+sigma_binding*1.96, color = Pos_ref_plot)) +
      ggplot2::geom_linerange(data = plot_dt[Pos_ref_plot!="Remainder"], alpha = 0.25, ggplot2::aes(xmin = fitness_abundance-sigma_abundance*1.96, xmax = fitness_abundance+sigma_abundance*1.96, color = Pos_ref_plot)) +
      ggplot2::geom_line(data = plot_dt[fitness_binding_zbm>=min(fitness_binding)], ggplot2::aes(fitness_abundance, fitness_binding_zbm), color = colour_scheme[["shade 0"]][[1]]) +
      ggplot2::geom_vline(xintercept = 0) +
      ggplot2::geom_abline(linetype = 2) +
      ggplot2::geom_hline(yintercept = 0) +
      ggplot2::xlab("Fitness (Abundance)") +
      ggplot2::ylab("Fitness (Binding)") +
      ggplot2::labs(color = "Allosteric") +
      ggplot2::scale_colour_manual(values=plot_cols) +
      ggplot2::theme_classic()
    ggplot2::ggsave(file.path(outpath, paste0("fitness_scatter_allosteric_mutations_", i, ".pdf")), d, width = 4, height = 3, useDingbats=FALSE)
  }

}

