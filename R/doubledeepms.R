
#' doubledeepms
#'
#' Main analysis script.
#'
#' @param startStage Start at a specified analysis stage (default:1)
#' @param stopStage Stop at a specified analysis stage (default:0 i.e. no stop condition)
#' @param base_dir Base directory for all output file (default:NB private CRG server path; change accordingly)
#' @param rerun_raw re-run raw data preprocessing (default:F)
#' @param Ncores # cores available to parallelize leave-one-out cross validation workflow (stage 4)
#'
#' @return Nothing
#' @export
doubledeepms <- function(
  startStage=1,
  stopStage=0,
  base_dir = "/users/project/prj004631/afaure/DMS/Results/doubledeepms_proj",
  rerun_raw = F,
  Ncores = 1
  ){

  # startStage=1
  # stopStage=0
  # base_dir = "/users/project/prj004631/afaure/DMS/Results/doubledeepms_proj"
  # rerun_raw = F
  # Ncores = 1

  colour_scheme <- list(
    "shade 0" = list(
      "#F4270C",
      "#F4AD0C",
      "#1B38A6",
      "#09B636"),
    "shade 1" = list(
      "#FFB0A5",
      "#FFE4A5",
      "#9DACE3",
      "#97E9AD"),
    "shade 2" = list(
      "#FF6A56",
      "#FFCB56",
      "#4C63B7",
      "#43C766"),
    "shade 3" = list(
      "#A31300",
      "#A37200",
      "#0C226F",
      "#007A20"),
    "shade 4" = list(
      "#410800",
      "#412D00",
      "#020B2C",
      "#00300D"),
    "shade 5" = list(
      "#CCCCCC",
      "#FF991F",
      "#5CB8FF",
      "#B22222"
    ))

  #First and last analysis stages
  first_stage <- startStage
  last_stage <- stopStage

  ### Fit thermo model
  ###########################

  ### Evaluate thermo model results
  ###########################

  stagenum <- 1
  #GB1
  doubledeepms_thermo_model_results(
    mochi_outpath = file.path(base_dir, "Data", "mochi", "GB1", "mochi__fit_tmodel_3state"),
    literature_free_energies = file.path(base_dir, "Data", "in_vitro", "GB1_literature_free_energies.txt"),
    position_offset = 1,
    outpath = doubledeepms__format_dir(dir_suffix="_doubledeepms_thermo_model_results_GB1", stagenum=stagenum, base_dir=base_dir),
    colour_scheme = colour_scheme,
    execute = (first_stage <= stagenum & (last_stage == 0 | last_stage >= stagenum)))
  #PSD95-PDZ3
  doubledeepms_thermo_model_results(
    mochi_outpath = file.path(base_dir, "Data", "mochi", "PSD95-PDZ3", "mochi__fit_tmodel_3state"),
    literature_free_energies = file.path(base_dir, "Data", "in_vitro", "PDZ_literature_free_energies.txt"),
    position_offset = 310,
    folding_energy_max_sd = 0.37,
    binding_energy_max_sd = 0.5,
    outpath = doubledeepms__format_dir(dir_suffix="_doubledeepms_thermo_model_results_PSD95-PDZ3", stagenum=stagenum, base_dir=base_dir),
    colour_scheme = colour_scheme,
    execute = (first_stage <= stagenum & (last_stage == 0 | last_stage >= stagenum)))
  #GRB2-SH3
  doubledeepms_thermo_model_results(
    mochi_outpath = file.path(base_dir, "Data", "mochi", "GRB2-SH3", "mochi__fit_tmodel_3state"),
    position_offset = 0,
    outpath = doubledeepms__format_dir(dir_suffix="_doubledeepms_thermo_model_results_GRB2-SH3", stagenum=stagenum, base_dir=base_dir),
    colour_scheme = colour_scheme,
    execute = (first_stage <= stagenum & (last_stage == 0 | last_stage >= stagenum)))

  ### Add 3D structure metrics
  ###########################

  stagenum <- 2
  #GB1
  doubledeepms_structure_metrics(
    input_file = file.path(base_dir, paste0("001", "_doubledeepms_thermo_model_results_GB1"), "dg_singles.txt"),
    outpath = doubledeepms__format_dir(dir_suffix="_doubledeepms_structure_metrics_GB1", stagenum=stagenum, base_dir=base_dir),
    pdb_file = file.path(base_dir, "Data", "pdb", "1fcc.pdb"),
    pdb_chain_query = "C",
    pdb_chain_target = "A",
    execute = (first_stage <= stagenum & (last_stage == 0 | last_stage >= stagenum)))
  #PSD95-PDZ3
  doubledeepms_structure_metrics(
    input_file = file.path(base_dir, paste0("001", "_doubledeepms_thermo_model_results_PSD95-PDZ3"), "dg_singles.txt"),
    outpath = doubledeepms__format_dir(dir_suffix="_doubledeepms_structure_metrics_PSD95-PDZ3", stagenum=stagenum, base_dir=base_dir),
    pdb_file = file.path(base_dir, "Data", "pdb", "1be9.pdb"),
    execute = (first_stage <= stagenum & (last_stage == 0 | last_stage >= stagenum)))
  #GRB2-SH3
  doubledeepms_structure_metrics(
    input_file = file.path(base_dir, paste0("001", "_doubledeepms_thermo_model_results_GRB2-SH3"), "dg_singles.txt"),
    outpath = doubledeepms__format_dir(dir_suffix="_doubledeepms_structure_metrics_GRB2-SH3", stagenum=stagenum, base_dir=base_dir),
    pdb_file = file.path(base_dir, "Data", "pdb", "2vwf.pdb"),
    execute = (first_stage <= stagenum & (last_stage == 0 | last_stage >= stagenum)))

  ### Fitness plots
  ###########################

  stagenum <- 3
  #PSD95-PDZ3 and GRB2-SH3
  doubledeepms_fitness_plots(
    fitness_list = list(
      "PSD95-PDZ3" = file.path(base_dir, "Data", "fitness", "PSD95-PDZ3"),
      "GRB2-SH3" = file.path(base_dir, "Data", "fitness", "GRB2-SH3")),
    structure_metrics_list = list(
      "PSD95-PDZ3" = file.path(base_dir, paste0("002", "_doubledeepms_structure_metrics_PSD95-PDZ3"), "dg_singles.txt"),
      "GRB2-SH3" = file.path(base_dir, paste0("002", "_doubledeepms_structure_metrics_GRB2-SH3"), "dg_singles.txt")),
    scatter_examples_list = list(
      "GRB2-SH3" = c("Abundance", "singles", "rep2", "rep1"),
      "PSD95-PDZ3" = c("Binding", "doubles", "rep2", "rep3")),
    val_inpath = file.path(base_dir, "Data", "experimental_validations", "GRB2-SH3", "experimenal_validations_ODs.txt"),
    outpath = doubledeepms__format_dir(dir_suffix="_doubledeepms_fitness_plots", stagenum=stagenum, base_dir=base_dir),
    colour_scheme = colour_scheme,
    execute = (first_stage <= stagenum & (last_stage == 0 | last_stage >= stagenum)))
  
  ### Plot free fitness heatmaps
  ###########################

  stagenum <- 4
  #PSD95-PDZ3
  doubledeepms_fitness_heatmaps(
    input_file = file.path(base_dir, paste0("002", "_doubledeepms_structure_metrics_PSD95-PDZ3"), "dg_singles.txt"),
    input_file_fitness = file.path(base_dir, "Data", "fitness", "PSD95-PDZ3"),
    domain_name = "PSD95 PDZ3",
    outpath = doubledeepms__format_dir(dir_suffix="_doubledeepms_fitness_heatmaps_PSD95-PDZ3", stagenum=stagenum, base_dir=base_dir),
    colour_scheme = colour_scheme,
    plot_width = 15,
    execute = (first_stage <= stagenum & (last_stage == 0 | last_stage >= stagenum)))
  #GRB2-SH3
  doubledeepms_fitness_heatmaps(
    input_file = file.path(base_dir, paste0("002", "_doubledeepms_structure_metrics_GRB2-SH3"), "dg_singles.txt"),
    input_file_fitness = file.path(base_dir, "Data", "fitness", "GRB2-SH3"),
    domain_name = "GRB2 SH3",
    outpath = doubledeepms__format_dir(dir_suffix="_doubledeepms_fitness_heatmaps_GRB2-SH3", stagenum=stagenum, base_dir=base_dir),
    colour_scheme = colour_scheme,
    plot_width = 11,
    execute = (first_stage <= stagenum & (last_stage == 0 | last_stage >= stagenum)))

  ### Plot free energy scatterplots
  ###########################

  stagenum <- 5
  #All domains
  doubledeepms_free_energy_scatterplots(
    input_list = list(
      "GB1" = file.path(base_dir, paste0("002", "_doubledeepms_structure_metrics_GB1"), "dg_singles.txt"),
      "PSD95-PDZ3" = file.path(base_dir, paste0("002", "_doubledeepms_structure_metrics_PSD95-PDZ3"), "dg_singles.txt"),
      "GRB2-SH3" = file.path(base_dir, paste0("002", "_doubledeepms_structure_metrics_GRB2-SH3"), "dg_singles.txt")),
    outpath = doubledeepms__format_dir(dir_suffix="_doubledeepms_free_energy_scatterplots", stagenum=stagenum, base_dir=base_dir),
    colour_scheme = colour_scheme,
    execute = (first_stage <= stagenum & (last_stage == 0 | last_stage >= stagenum)))

  ### Plot free energy heatmaps
  ###########################

  stagenum <- 6
  #GB1
  doubledeepms_free_energy_heatmaps(
    input_file = file.path(base_dir, paste0("002", "_doubledeepms_structure_metrics_GB1"), "dg_singles.txt"),
    domain_name = "GB1",
    outpath = doubledeepms__format_dir(dir_suffix="_doubledeepms_free_energy_heatmaps_GB1", stagenum=stagenum, base_dir=base_dir),
    colour_scheme = colour_scheme,
    plot_width = 11,
    execute = (first_stage <= stagenum & (last_stage == 0 | last_stage >= stagenum)))
  #PSD95-PDZ3
  doubledeepms_free_energy_heatmaps(
    input_file = file.path(base_dir, paste0("002", "_doubledeepms_structure_metrics_PSD95-PDZ3"), "dg_singles.txt"),
    domain_name = "PSD95 PDZ3",
    outpath = doubledeepms__format_dir(dir_suffix="_doubledeepms_free_energy_heatmaps_PSD95-PDZ3", stagenum=stagenum, base_dir=base_dir),
    colour_scheme = colour_scheme,
    plot_width = 15,
    execute = (first_stage <= stagenum & (last_stage == 0 | last_stage >= stagenum)))
  #GRB2-SH3
  doubledeepms_free_energy_heatmaps(
    input_file = file.path(base_dir, paste0("002", "_doubledeepms_structure_metrics_GRB2-SH3"), "dg_singles.txt"),
    domain_name = "GRB2 SH3",
    outpath = doubledeepms__format_dir(dir_suffix="_doubledeepms_free_energy_heatmaps_GRB2-SH3", stagenum=stagenum, base_dir=base_dir),
    colour_scheme = colour_scheme,
    plot_width = 11,
    execute = (first_stage <= stagenum & (last_stage == 0 | last_stage >= stagenum)))

  ### Protein stability plots
  ###########################

  stagenum <- 7
  #All domains
  doubledeepms_protein_stability_plots(
    input_list = list(
      "GB1" = file.path(base_dir, paste0("002", "_doubledeepms_structure_metrics_GB1"), "dg_singles.txt"),
      "PSD95-PDZ3" = file.path(base_dir, paste0("002", "_doubledeepms_structure_metrics_PSD95-PDZ3"), "dg_singles.txt"),
      "GRB2-SH3" = file.path(base_dir, paste0("002", "_doubledeepms_structure_metrics_GRB2-SH3"), "dg_singles.txt")),
    pdb_file_list = list(
      "GB1" = file.path(base_dir, "Data", "pdb", "1fcc.pdb"),
      "PSD95-PDZ3" = file.path(base_dir, "Data", "pdb", "1be9.pdb"),
      "GRB2-SH3" = file.path(base_dir, "Data", "pdb", "2vwf.pdb")),
    pdb_chain_query_list = list(
      "GB1" = "C",
      "PSD95-PDZ3" = "A",
      "GRB2-SH3" = "A"),
    aaprop_file = file.path(base_dir, "Data", "amino_acid_properties", "amino_acid_properties_annotated_supplementary.txt"),
    aaprop_file_selected = file.path(base_dir, "Data", "amino_acid_properties", "selected.amino_acid_properties.txt"),
    outpath = doubledeepms__format_dir(dir_suffix="_doubledeepms_protein_stability_plots", stagenum=stagenum, base_dir=base_dir),
    colour_scheme = colour_scheme,
    execute = (first_stage <= stagenum & (last_stage == 0 | last_stage >= stagenum)))

  ### Binding interface plots (Julia to add)
  ###########################

  stagenum <- 8
  # GRB2-SH3
  doubledeepms_interface_mechanisms(
    base_dir = base_dir, 
    domain_name = "GRB2-SH3", 
    outpath = doubledeepms__format_dir(dir_suffix="_doubledeepms_ligand_interface_plots", stagenum=stagenum, base_dir=base_dir), 
    colour_scheme = colour_scheme,
    mut_subset_list = list(
      c("T", "F","Q","H", "S", "V", "I", "C"),
      c("T", "S", "L", "V", "A", "K", "R")),
    pos_subset_list = list(
      c(7, 51),
      c(46)),
    ligand_pos_list = list(
      c(3,4),
      c(11)),
    pdb_id = "2vwf",
    pdb_view = "(0.837024331,0.261427671,0.480670929,0.209657252,0.658192515,-0.723066926,-0.505404055,0.706001341,0.496113181,0.000000000,0.000000000,-141.782043457,-2.309909821,13.775757790,0.139801025,111.782043457,171.782043457,-20.000000000 )",
    )

  ### Allostery plots
  ###########################

  stagenum <- 9
  #All domains
  doubledeepms_allostery_plots(
    input_list = list(
      "GB1" = file.path(base_dir, paste0("002", "_doubledeepms_structure_metrics_GB1"), "dg_singles.txt"),
      "PSD95-PDZ3" = file.path(base_dir, paste0("002", "_doubledeepms_structure_metrics_PSD95-PDZ3"), "dg_singles.txt"),
      "GRB2-SH3" = file.path(base_dir, paste0("002", "_doubledeepms_structure_metrics_GRB2-SH3"), "dg_singles.txt")),
    pdb_file_list = list(
      "GB1" = file.path(base_dir, "Data", "pdb", "1fcc.pdb"),
      "PSD95-PDZ3" = file.path(base_dir, "Data", "pdb", "1be9.pdb"),
      "GRB2-SH3" = file.path(base_dir, "Data", "pdb", "2vwf.pdb")),
    pdb_chain_query_list = list(
      "GB1" = "C",
      "PSD95-PDZ3" = "A",
      "GRB2-SH3" = "A"),
    literature_list = list(
      "GB1" = c(),
      "PSD95-PDZ3" = c(322, 327, 329, 330, 336, 362, 372, 375, 379),
      "GRB2-SH3" = c()),
    aaprop_file = file.path(base_dir, "Data", "amino_acid_properties", "amino_acid_properties_annotated_supplementary.txt"),
    aaprop_file_selected = file.path(base_dir, "Data", "amino_acid_properties", "selected.amino_acid_properties.txt"),
    outpath = doubledeepms__format_dir(dir_suffix="_doubledeepms_allostery_plots", stagenum=stagenum, base_dir=base_dir),
    colour_scheme = colour_scheme,
    execute = (first_stage <= stagenum & (last_stage == 0 | last_stage >= stagenum)))

}



