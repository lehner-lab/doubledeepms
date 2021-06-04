
#' doubledeepms__allosteric_mutations_scatterplot
#'
#' Allosteric mutations scatterplot.
#'
#' @param input_dt input data table (required)
#' @param outpath plot output path (required)
#' @param colour_scheme colour scheme file (required)
#'
#' @return Allosteric mutations
#' @export
doubledeepms__allosteric_mutations_scatterplot <- function(
  input_dt,
  outpath,
  colour_scheme
  ){

  #Outlier changes in binding free energy
  reg_threshold <- input_dt[b_ddg_pred_conf==T & Pos_class=="binding_interface"][!duplicated(Pos_ref)][order(b_ddg_posmeanabs_conf, decreasing = T)][5][,b_ddg_posmeanabs_conf]
  input_dt[b_ddg_pred_conf==T, b_ddg_pred_outlier := p.adjust(doubledeepms__pvalue(abs(b_ddg_pred)-reg_threshold, b_ddg_pred_sd), method = "BH")<0.05 & (abs(b_ddg_pred)-reg_threshold)>0]

  # #Outlier changes in binding free energy
  # reg_threshold <- input_dt[b_ddg_pred_conf==T, mean(b_ddg_pred)]
  # input_dt[b_ddg_pred_conf==T, b_ddg_pred_outlier := p.adjust(doubledeepms__pvalue(b_ddg_pred-reg_threshold, b_ddg_pred_sd), method = "BH")<0.05]

  #Allosteric mutations (not within binding interface and not at allosteric site)
  input_dt[b_ddg_pred_conf==T & is.na(allosteric), allosteric_mutation := b_ddg_pred_outlier]

  #Scatterplot
  plot_dt <- input_dt[b_ddg_pred_conf==T]
  plot_dt[, Pos_class_plot := "Remainder"]
  plot_dt[allosteric_mutation & Pos_class=="binding_interface", Pos_class_plot := "Orthosteric mutation"]
  plot_dt[allosteric_mutation & Pos_class=="core", Pos_class_plot := "Core allosteric mutation"]
  plot_dt[allosteric_mutation & Pos_class=="surface", Pos_class_plot := "Surface allosteric mutation"]
  plot_dt[!is.na(allosteric) & Pos_class=="binding_interface", Pos_class_plot := "Orthosteric site"]
  plot_dt[!is.na(allosteric) & Pos_class!="binding_interface", Pos_class_plot := "Allosteric site"]
  plot_dt[, Pos_class_plot := factor(Pos_class_plot, levels = c("Remainder", "Orthosteric site", "Orthosteric mutation", "Allosteric site", "Core allosteric mutation", "Surface allosteric mutation"))]
  plot_cols = c("grey", colour_scheme[["shade 0"]][c(1,1:4)])
  names(plot_cols) <- c("Remainder", "Orthosteric site", "Orthosteric mutation", "Allosteric site", "Core allosteric mutation", "Surface allosteric mutation")
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(Pos_ref, b_ddg_pred, color = Pos_class_plot)) +
    ggplot2::geom_point(data = plot_dt[Pos_class_plot=="Remainder"], color = "grey", size = 0.25) +
    ggplot2::geom_point(data = plot_dt[Pos_class_plot %in% c("Allosteric site", "Orthosteric site")], size = 0.75, shape = 21) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2) +
    ggplot2::geom_point(data = plot_dt[grepl("mutation", Pos_class_plot)], size = 1.5) +
    ggrepel::geom_text_repel(data = plot_dt[allosteric_mutation==T & Pos_class=="surface"], ggplot2::aes(label = id_ref, color = Pos_class_plot), show.legend = F, 
      max.overlaps = Inf) +
    ggplot2::xlab("Amino acid position") +
    ggplot2::ylab(expression("Binding "*Delta*Delta*"G")) +
    ggplot2::theme_classic() +
    ggplot2::scale_colour_manual(values=plot_cols) +
    ggplot2::labs(color = "Mutation\ntype")   
  suppressWarnings(ggplot2::ggsave(outpath, d, width = 6, height = 3, useDingbats=FALSE))

  return(input_dt)

  # #Mutant delta charges
  # resi_codes <- unlist(strsplit("GAVLMIFYWKRHDESTCNQP", ""))
  # resi_charge <- rep(0, 20)
  # names(resi_charge) <- resi_codes
  # resi_charge[resi_codes %in% c("R", "H", "K")] <- 1
  # resi_charge[resi_codes %in% c("D", "E")] <- -1
  # input_dt[, WT_AA_charge := resi_charge[WT_AA]]
  # input_dt[, Mut_charge := resi_charge[Mut]]
  # input_dt[, delta_charge := Mut_charge-WT_AA_charge]
  # input_dt[, abs_delta_charge := abs(Mut_charge)-abs(WT_AA_charge)]

  # #Enrichment of allosteric mutations in binding interface
  # temp_in <- input_dt[b_ddg_pred_conf==T & b_ddg_pred_outlier,Pos_class=="binding_interface"]
  # temp_out <- input_dt[b_ddg_pred_conf==T & !b_ddg_pred_outlier,Pos_class=="binding_interface"]
  # fisher.test(matrix(c(sum(temp_in), sum(!temp_in), sum(temp_out), sum(!temp_out)), nrow = 2))

  # #Enrichment of allosteric mutations in allosteric sites
  # temp_in <- input_dt[b_ddg_pred_conf==T & b_ddg_pred_outlier,!is.na(allosteric)]
  # temp_out <- input_dt[b_ddg_pred_conf==T & !b_ddg_pred_outlier,!is.na(allosteric)]
  # fisher.test(matrix(c(sum(temp_in), sum(!temp_in), sum(temp_out), sum(!temp_out)), nrow = 2))

  # #Enrichment of allosteric mutations (not within binding interface and not at allosteric site) in core versus surface
  # temp_in <- input_dt[Pos_class!="binding_interface" & allosteric_mutation,Pos_class=="core"]
  # temp_out <- input_dt[Pos_class!="binding_interface" & !allosteric_mutation,Pos_class=="core"]
  # fisher.test(matrix(c(sum(temp_in), sum(!temp_in), sum(temp_out), sum(!temp_out)), nrow = 2))

  # #Enrichment of allosteric mutations (not within binding interface and not at allosteric site) in proximity to ligand
  # temp_in <- input_dt[Pos_class!="binding_interface" & allosteric_mutation,scHAmin_ligand]
  # temp_out <- input_dt[Pos_class!="binding_interface" & !allosteric_mutation,scHAmin_ligand]
  # doubledeepms__mann_whitney_U_wrapper(temp_out, temp_in)

  # #Enrichment of allosteric mutations at surface that involve change in charge
  # temp_in <- input_dt[allosteric_mutation & Pos_class=="surface",delta_charge!=0]
  # temp_out <- input_dt[b_ddg_pred_conf==T & is.na(allosteric) & Pos_class=="surface",delta_charge!=0]
  # fisher.test(matrix(c(sum(temp_in), sum(!temp_in), sum(temp_out), sum(!temp_out)), nrow = 2))
  # print(sum(temp_in)/length(temp_in))

  # #Test change in mass
  # temp_in <- input_dt[b_ddg_pred_conf==T & allosteric_mutation][,delta_mass]
  # temp_out <- input_dt[b_ddg_pred_conf==T & !allosteric_mutation][,delta_mass]
  # fisher.test(matrix(c(sum(temp_in>0), sum(temp_in<=0), sum(temp_out>0), sum(temp_out<=0)), nrow = 2))
  # doubledeepms__mann_whitney_U_wrapper(temp_in, temp_out)

  # #Test change in mass excluding charged WT AAs
  # temp_in <- input_dt[b_ddg_pred_conf==T & is.na(allosteric) & allosteric_mutation][,delta_mass]
  # temp_out <- input_dt[b_ddg_pred_conf==T & is.na(allosteric) & !allosteric_mutation][,delta_mass]
  # fisher.test(matrix(c(sum(temp_in>0), sum(temp_in<=0), sum(temp_out>0), sum(temp_out<=0)), nrow = 2))
  # doubledeepms__mann_whitney_U_wrapper(temp_in, temp_out)

  # #Test change in mass excluding charged WT AAs and excluding glycine
  # temp_in <- input_dt[b_ddg_pred_conf==T & WT_AA_charge==0 & WT_AA!="G" & allosteric_mutation][,delta_mass]
  # temp_out <- input_dt[b_ddg_pred_conf==T & WT_AA_charge==0 & WT_AA!="G" & !allosteric_mutation][,delta_mass]
  # fisher.test(matrix(c(sum(temp_in>0), sum(temp_in<=0), sum(temp_out>0), sum(temp_out<=0)), nrow = 2))
  # doubledeepms__mann_whitney_U_wrapper(temp_in, temp_out)

  # #Mutant delta masses
  # resi_codes <- unlist(strsplit("GAVLMIFYWKRHDESTCNQP", ""))
  # resi_names <- bio3d::aa123(resi_codes)
  # resi_masses <- bio3d::aa2mass(resi_names, addter=FALSE)
  # names(resi_masses) <- resi_codes
  # input_dt[, WT_AA_mass := resi_masses[WT_AA]]
  # input_dt[, Mut_mass := resi_masses[Mut]]
  # input_dt[, delta_mass := Mut_mass-WT_AA_mass]
}
