
#' doubledeepms
#'
#' Main analysis script.
#'
#' @param startStage Start at a specified analysis stage (default:1)
#' @param stopStage Stop at a specified analysis stage (default:15)
#' @param base_dir Base directory for all output files (default:NB private CRG server path; change accordingly)
#' @param tmodel_job_number Thermodynamic model fit job number: 1:final model, 2-11: monte carlo simluations for confidence intervals of model-inferred free energies (default:1)
#' @param tmodel_grid_search Thermodynamic model fit grid search to determine optimal hyperparameters (default:FALSE)
#' @param tmodel_protein Thermodynamic model fit proteins: comma-separated list of names or 'all' (default:'GRB2-SH3')
#' @param tmodel_subset Thermodynamic model fit data subset: either an integer percentage of subsampled doubles (1-100) or "binding_only" or "singles_only" (default:100)
#'
#' @return Nothing
#' @export
doubledeepms <- function(
  startStage = 1,
  stopStage = 15,
  base_dir = "/users/project/prj004631/afaure/DMS/Results/doubledeepms_proj",
  tmodel_job_number = 1,
  tmodel_grid_search = FALSE,
  tmodel_protein = "GRB2-SH3",
  tmodel_subset = 100
  ){

  # startStage=1
  # stopStage=15
  # base_dir = "/users/project/prj004631/afaure/DMS/Results/doubledeepms_proj"
  # rerun_raw = F

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

  ### Fit thermo models
  ###########################

  stagenum <- 0
  doubledeepms_fit_thermo_model(
    base_dir = base_dir,
    tmodel_job_number = tmodel_job_number,
    tmodel_grid_search = tmodel_grid_search,
    tmodel_protein = tmodel_protein,
    tmodel_subset = as.character(tmodel_subset),
    tmodel_hyperparameters = file.path(base_dir, "Data", "mochi", "hyperparameters.txt"),
    execute = (first_stage <= stagenum & last_stage >= stagenum))

  ### Evaluate thermo model results
  ###########################

  stagenum <- 1
  #GB1
  doubledeepms_thermo_model_results(
    mochi_outpath = file.path(base_dir, "Data", "mochi", "GB1", "mochi__fit_tmodel_3state_sparse_dimsum128"),
    literature_free_energies = file.path(base_dir, "Data", "in_vitro", "GB1_literature_free_energies.txt"),
    position_offset = 1,
    outpath = doubledeepms__format_dir(dir_suffix="_doubledeepms_thermo_model_results_GB1", stagenum=stagenum, base_dir=base_dir),
    colour_scheme = colour_scheme,
    execute = (first_stage <= stagenum & last_stage >= stagenum))
  #GB1 - subsample doubles 2%
  doubledeepms_thermo_model_results(
    mochi_outpath = file.path(base_dir, "Data", "mochi", "GB1", "mochi__fit_tmodel_3state_sparse_dimsum128_subsample2p"),
    literature_free_energies = file.path(base_dir, "Data", "in_vitro", "GB1_literature_free_energies.txt"),
    position_offset = 1,
    outpath = doubledeepms__format_dir(dir_suffix="_doubledeepms_thermo_model_results_GB1subsample2p", stagenum=stagenum, base_dir=base_dir),
    colour_scheme = colour_scheme,
    execute = (first_stage <= stagenum & last_stage >= stagenum))
  #GB1 - subsample doubles 5%
  doubledeepms_thermo_model_results(
    mochi_outpath = file.path(base_dir, "Data", "mochi", "GB1", "mochi__fit_tmodel_3state_sparse_dimsum128_subsample5p"),
    literature_free_energies = file.path(base_dir, "Data", "in_vitro", "GB1_literature_free_energies.txt"),
    position_offset = 1,
    outpath = doubledeepms__format_dir(dir_suffix="_doubledeepms_thermo_model_results_GB1subsample5p", stagenum=stagenum, base_dir=base_dir),
    colour_scheme = colour_scheme,
    execute = (first_stage <= stagenum & last_stage >= stagenum))
  #GB1 - subsample doubles 10%
  doubledeepms_thermo_model_results(
    mochi_outpath = file.path(base_dir, "Data", "mochi", "GB1", "mochi__fit_tmodel_3state_sparse_dimsum128_subsample10p"),
    literature_free_energies = file.path(base_dir, "Data", "in_vitro", "GB1_literature_free_energies.txt"),
    position_offset = 1,
    outpath = doubledeepms__format_dir(dir_suffix="_doubledeepms_thermo_model_results_GB1subsample10p", stagenum=stagenum, base_dir=base_dir),
    colour_scheme = colour_scheme,
    execute = (first_stage <= stagenum & last_stage >= stagenum))
  #GB1 - subsample doubles 20%
  doubledeepms_thermo_model_results(
    mochi_outpath = file.path(base_dir, "Data", "mochi", "GB1", "mochi__fit_tmodel_3state_sparse_dimsum128_subsample20p"),
    literature_free_energies = file.path(base_dir, "Data", "in_vitro", "GB1_literature_free_energies.txt"),
    position_offset = 1,
    outpath = doubledeepms__format_dir(dir_suffix="_doubledeepms_thermo_model_results_GB1subsample20p", stagenum=stagenum, base_dir=base_dir),
    colour_scheme = colour_scheme,
    execute = (first_stage <= stagenum & last_stage >= stagenum))
  #GB1 - subsample doubles 50%
  doubledeepms_thermo_model_results(
    mochi_outpath = file.path(base_dir, "Data", "mochi", "GB1", "mochi__fit_tmodel_3state_sparse_dimsum128_subsample50p"),
    literature_free_energies = file.path(base_dir, "Data", "in_vitro", "GB1_literature_free_energies.txt"),
    position_offset = 1,
    outpath = doubledeepms__format_dir(dir_suffix="_doubledeepms_thermo_model_results_GB1subsample50p", stagenum=stagenum, base_dir=base_dir),
    colour_scheme = colour_scheme,
    execute = (first_stage <= stagenum & last_stage >= stagenum))
  #PSD95-PDZ3
  doubledeepms_thermo_model_results(
    mochi_outpath = file.path(base_dir, "Data", "mochi", "PSD95-PDZ3", "mochi__fit_tmodel_3state_sparse_dimsum128"),
    literature_free_energies = file.path(base_dir, "Data", "in_vitro", "PDZ_literature_free_energies.txt"),
    position_offset = 310,
    outpath = doubledeepms__format_dir(dir_suffix="_doubledeepms_thermo_model_results_PSD95-PDZ3", stagenum=stagenum, base_dir=base_dir),
    colour_scheme = colour_scheme,
    execute = (first_stage <= stagenum & last_stage >= stagenum))
  #PSD95-PDZ3 - binding only
  doubledeepms_thermo_model_results(
    mochi_outpath = file.path(base_dir, "Data", "mochi", "PSD95-PDZ3", "mochi__fit_tmodel_3state_sparse_dimsum128_bindingonly"),
    literature_free_energies = file.path(base_dir, "Data", "in_vitro", "PDZ_literature_free_energies.txt"),
    position_offset = 310,
    outpath = doubledeepms__format_dir(dir_suffix="_doubledeepms_thermo_model_results_PSD95-PDZ3bindingonly", stagenum=stagenum, base_dir=base_dir),
    colour_scheme = colour_scheme,
    execute = (first_stage <= stagenum & last_stage >= stagenum))
  #PSD95-PDZ3 - singles only
  doubledeepms_thermo_model_results(
    mochi_outpath = file.path(base_dir, "Data", "mochi", "PSD95-PDZ3", "mochi__fit_tmodel_3state_sparse_dimsum128_singlesonly"),
    literature_free_energies = file.path(base_dir, "Data", "in_vitro", "PDZ_literature_free_energies.txt"),
    position_offset = 310,
    outpath = doubledeepms__format_dir(dir_suffix="_doubledeepms_thermo_model_results_PSD95-PDZ3singlesonly", stagenum=stagenum, base_dir=base_dir),
    colour_scheme = colour_scheme,
    execute = (first_stage <= stagenum & last_stage >= stagenum))
  #PSD95-PDZ3 - subsample doubles 2%
  doubledeepms_thermo_model_results(
    mochi_outpath = file.path(base_dir, "Data", "mochi", "PSD95-PDZ3", "mochi__fit_tmodel_3state_sparse_dimsum128_subsample2p"),
    literature_free_energies = file.path(base_dir, "Data", "in_vitro", "PDZ_literature_free_energies.txt"),
    position_offset = 310,
    outpath = doubledeepms__format_dir(dir_suffix="_doubledeepms_thermo_model_results_PSD95-PDZ3subsample2p", stagenum=stagenum, base_dir=base_dir),
    colour_scheme = colour_scheme,
    execute = (first_stage <= stagenum & last_stage >= stagenum))
  #PSD95-PDZ3 - subsample doubles 5%
  doubledeepms_thermo_model_results(
    mochi_outpath = file.path(base_dir, "Data", "mochi", "PSD95-PDZ3", "mochi__fit_tmodel_3state_sparse_dimsum128_subsample5p"),
    literature_free_energies = file.path(base_dir, "Data", "in_vitro", "PDZ_literature_free_energies.txt"),
    position_offset = 310,
    outpath = doubledeepms__format_dir(dir_suffix="_doubledeepms_thermo_model_results_PSD95-PDZ3subsample5p", stagenum=stagenum, base_dir=base_dir),
    colour_scheme = colour_scheme,
    execute = (first_stage <= stagenum & last_stage >= stagenum))
  #PSD95-PDZ3 - subsample doubles 10%
  doubledeepms_thermo_model_results(
    mochi_outpath = file.path(base_dir, "Data", "mochi", "PSD95-PDZ3", "mochi__fit_tmodel_3state_sparse_dimsum128_subsample10p"),
    literature_free_energies = file.path(base_dir, "Data", "in_vitro", "PDZ_literature_free_energies.txt"),
    position_offset = 310,
    outpath = doubledeepms__format_dir(dir_suffix="_doubledeepms_thermo_model_results_PSD95-PDZ3subsample10p", stagenum=stagenum, base_dir=base_dir),
    colour_scheme = colour_scheme,
    execute = (first_stage <= stagenum & last_stage >= stagenum))
  #PSD95-PDZ3 - subsample doubles 20%
  doubledeepms_thermo_model_results(
    mochi_outpath = file.path(base_dir, "Data", "mochi", "PSD95-PDZ3", "mochi__fit_tmodel_3state_sparse_dimsum128_subsample20p"),
    literature_free_energies = file.path(base_dir, "Data", "in_vitro", "PDZ_literature_free_energies.txt"),
    position_offset = 310,
    outpath = doubledeepms__format_dir(dir_suffix="_doubledeepms_thermo_model_results_PSD95-PDZ3subsample20p", stagenum=stagenum, base_dir=base_dir),
    colour_scheme = colour_scheme,
    execute = (first_stage <= stagenum & last_stage >= stagenum))
  #PSD95-PDZ3 - subsample doubles 50%
  doubledeepms_thermo_model_results(
    mochi_outpath = file.path(base_dir, "Data", "mochi", "PSD95-PDZ3", "mochi__fit_tmodel_3state_sparse_dimsum128_subsample50p"),
    literature_free_energies = file.path(base_dir, "Data", "in_vitro", "PDZ_literature_free_energies.txt"),
    position_offset = 310,
    outpath = doubledeepms__format_dir(dir_suffix="_doubledeepms_thermo_model_results_PSD95-PDZ3subsample50p", stagenum=stagenum, base_dir=base_dir),
    colour_scheme = colour_scheme,
    execute = (first_stage <= stagenum & last_stage >= stagenum))
  #GRB2-SH3
  doubledeepms_thermo_model_results(
    mochi_outpath = file.path(base_dir, "Data", "mochi", "GRB2-SH3", "mochi__fit_tmodel_3state_sparse_dimsum128"),
    literature_free_energies = file.path(base_dir, "Data", "in_vitro", "GRB2_literature_free_energies.txt"),
    position_offset = 0,
    outpath = doubledeepms__format_dir(dir_suffix="_doubledeepms_thermo_model_results_GRB2-SH3", stagenum=stagenum, base_dir=base_dir),
    colour_scheme = colour_scheme,
    execute = (first_stage <= stagenum & last_stage >= stagenum))
  #GRB2-SH3 - binding only
  doubledeepms_thermo_model_results(
    mochi_outpath = file.path(base_dir, "Data", "mochi", "GRB2-SH3", "mochi__fit_tmodel_3state_sparse_dimsum128_bindingonly"),
    literature_free_energies = file.path(base_dir, "Data", "in_vitro", "GRB2_literature_free_energies.txt"),
    position_offset = 0,
    outpath = doubledeepms__format_dir(dir_suffix="_doubledeepms_thermo_model_results_GRB2-SH3bindingonly", stagenum=stagenum, base_dir=base_dir),
    colour_scheme = colour_scheme,
    execute = (first_stage <= stagenum & last_stage >= stagenum))
  #GRB2-SH3 - singles only
  doubledeepms_thermo_model_results(
    mochi_outpath = file.path(base_dir, "Data", "mochi", "GRB2-SH3", "mochi__fit_tmodel_3state_sparse_dimsum128_singlesonly"),
    literature_free_energies = file.path(base_dir, "Data", "in_vitro", "GRB2_literature_free_energies.txt"),
    position_offset = 0,
    outpath = doubledeepms__format_dir(dir_suffix="_doubledeepms_thermo_model_results_GRB2-SH3singlesonly", stagenum=stagenum, base_dir=base_dir),
    colour_scheme = colour_scheme,
    execute = (first_stage <= stagenum & last_stage >= stagenum))
  #GRB2-SH3 - subsample doubles 2%
  doubledeepms_thermo_model_results(
    mochi_outpath = file.path(base_dir, "Data", "mochi", "GRB2-SH3", "mochi__fit_tmodel_3state_sparse_dimsum128_subsample2p"),
    literature_free_energies = file.path(base_dir, "Data", "in_vitro", "GRB2_literature_free_energies.txt"),
    position_offset = 0,
    outpath = doubledeepms__format_dir(dir_suffix="_doubledeepms_thermo_model_results_GRB2-SH3subsample2p", stagenum=stagenum, base_dir=base_dir),
    colour_scheme = colour_scheme,
    execute = (first_stage <= stagenum & last_stage >= stagenum))
  #GRB2-SH3 - subsample doubles 5%
  doubledeepms_thermo_model_results(
    mochi_outpath = file.path(base_dir, "Data", "mochi", "GRB2-SH3", "mochi__fit_tmodel_3state_sparse_dimsum128_subsample5p"),
    literature_free_energies = file.path(base_dir, "Data", "in_vitro", "GRB2_literature_free_energies.txt"),
    position_offset = 0,
    outpath = doubledeepms__format_dir(dir_suffix="_doubledeepms_thermo_model_results_GRB2-SH3subsample5p", stagenum=stagenum, base_dir=base_dir),
    colour_scheme = colour_scheme,
    execute = (first_stage <= stagenum & last_stage >= stagenum))
  #GRB2-SH3 - subsample doubles 10%
  doubledeepms_thermo_model_results(
    mochi_outpath = file.path(base_dir, "Data", "mochi", "GRB2-SH3", "mochi__fit_tmodel_3state_sparse_dimsum128_subsample10p"),
    literature_free_energies = file.path(base_dir, "Data", "in_vitro", "GRB2_literature_free_energies.txt"),
    position_offset = 0,
    outpath = doubledeepms__format_dir(dir_suffix="_doubledeepms_thermo_model_results_GRB2-SH3subsample10p", stagenum=stagenum, base_dir=base_dir),
    colour_scheme = colour_scheme,
    execute = (first_stage <= stagenum & last_stage >= stagenum))
  #GRB2-SH3 - subsample doubles 20%
  doubledeepms_thermo_model_results(
    mochi_outpath = file.path(base_dir, "Data", "mochi", "GRB2-SH3", "mochi__fit_tmodel_3state_sparse_dimsum128_subsample20p"),
    literature_free_energies = file.path(base_dir, "Data", "in_vitro", "GRB2_literature_free_energies.txt"),
    position_offset = 0,
    outpath = doubledeepms__format_dir(dir_suffix="_doubledeepms_thermo_model_results_GRB2-SH3subsample20p", stagenum=stagenum, base_dir=base_dir),
    colour_scheme = colour_scheme,
    execute = (first_stage <= stagenum & last_stage >= stagenum))
  #GRB2-SH3 - subsample doubles 50%
  doubledeepms_thermo_model_results(
    mochi_outpath = file.path(base_dir, "Data", "mochi", "GRB2-SH3", "mochi__fit_tmodel_3state_sparse_dimsum128_subsample50p"),
    literature_free_energies = file.path(base_dir, "Data", "in_vitro", "GRB2_literature_free_energies.txt"),
    position_offset = 0,
    outpath = doubledeepms__format_dir(dir_suffix="_doubledeepms_thermo_model_results_GRB2-SH3subsample50p", stagenum=stagenum, base_dir=base_dir),
    colour_scheme = colour_scheme,
    execute = (first_stage <= stagenum & last_stage >= stagenum))

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
    execute = (first_stage <= stagenum & last_stage >= stagenum))
  #PSD95-PDZ3
  doubledeepms_structure_metrics(
    input_file = file.path(base_dir, paste0("001", "_doubledeepms_thermo_model_results_PSD95-PDZ3"), "dg_singles.txt"),
    outpath = doubledeepms__format_dir(dir_suffix="_doubledeepms_structure_metrics_PSD95-PDZ3", stagenum=stagenum, base_dir=base_dir),
    pdb_file = file.path(base_dir, "Data", "pdb", "1be9.pdb"),
    execute = (first_stage <= stagenum & last_stage >= stagenum))
  #GRB2-SH3
  doubledeepms_structure_metrics(
    input_file = file.path(base_dir, paste0("001", "_doubledeepms_thermo_model_results_GRB2-SH3"), "dg_singles.txt"),
    outpath = doubledeepms__format_dir(dir_suffix="_doubledeepms_structure_metrics_GRB2-SH3", stagenum=stagenum, base_dir=base_dir),
    pdb_file = file.path(base_dir, "Data", "pdb", "2vwf.pdb"),
    execute = (first_stage <= stagenum & last_stage >= stagenum))

  ### Fitness plots
  ###########################

  stagenum <- 3
  #PSD95-PDZ3 and GRB2-SH3
  doubledeepms_fitness_plots(
    fitness_list = list(
      "GB1" = file.path(base_dir, "Data", "fitness", "GB1"),
      "PSD95-PDZ3" = file.path(base_dir, "Data", "fitness", "PSD95-PDZ3"),
      "GRB2-SH3" = file.path(base_dir, "Data", "fitness", "GRB2-SH3")),
    structure_metrics_list = list(
      "GB1" = file.path(base_dir, paste0("002", "_doubledeepms_structure_metrics_GB1"), "dg_singles.txt"),
      "PSD95-PDZ3" = file.path(base_dir, paste0("002", "_doubledeepms_structure_metrics_PSD95-PDZ3"), "dg_singles.txt"),
      "GRB2-SH3" = file.path(base_dir, paste0("002", "_doubledeepms_structure_metrics_GRB2-SH3"), "dg_singles.txt")),
    val_inpath = file.path(base_dir, "Data", "experimental_validations", "GRB2-SH3", "experimenal_validations_ODs.txt"),
    outpath = doubledeepms__format_dir(dir_suffix="_doubledeepms_fitness_plots", stagenum=stagenum, base_dir=base_dir),
    colour_scheme = colour_scheme,
    execute = (first_stage <= stagenum & last_stage >= stagenum))
  
  ### Plot fitness heatmaps
  ###########################

  stagenum <- 4
  #GB1
  doubledeepms_fitness_heatmaps(
    input_file = file.path(base_dir, paste0("002", "_doubledeepms_structure_metrics_GB1"), "dg_singles.txt"),
    input_file_fitness = file.path(base_dir, "Data", "fitness", "GB1"),
    domain_name = "GB1",
    outpath = doubledeepms__format_dir(dir_suffix="_doubledeepms_fitness_heatmaps_GB1", stagenum=stagenum, base_dir=base_dir),
    colour_scheme = colour_scheme,
    plot_width = 11,
    plot_traits = "Binding",
    execute = (first_stage <= stagenum & last_stage >= stagenum))
  #PSD95-PDZ3
  doubledeepms_fitness_heatmaps(
    input_file = file.path(base_dir, paste0("002", "_doubledeepms_structure_metrics_PSD95-PDZ3"), "dg_singles.txt"),
    input_file_fitness = file.path(base_dir, "Data", "fitness", "PSD95-PDZ3"),
    input_file_MSA = file.path(base_dir, "Data", "MSA", "PSD95-PDZ3", "frequencies.csv"),
    domain_name = "PSD95 PDZ3",
    outpath = doubledeepms__format_dir(dir_suffix="_doubledeepms_fitness_heatmaps_PSD95-PDZ3", stagenum=stagenum, base_dir=base_dir),
    colour_scheme = colour_scheme,
    plot_width = 15,
    execute = (first_stage <= stagenum & last_stage >= stagenum))
  #GRB2-SH3
  doubledeepms_fitness_heatmaps(
    input_file = file.path(base_dir, paste0("002", "_doubledeepms_structure_metrics_GRB2-SH3"), "dg_singles.txt"),
    input_file_fitness = file.path(base_dir, "Data", "fitness", "GRB2-SH3"),
    input_file_MSA = file.path(base_dir, "Data", "MSA", "GRB2-SH3", "frequencies.csv"),
    domain_name = "GRB2 SH3",
    outpath = doubledeepms__format_dir(dir_suffix="_doubledeepms_fitness_heatmaps_GRB2-SH3", stagenum=stagenum, base_dir=base_dir),
    colour_scheme = colour_scheme,
    plot_width = 11,
    execute = (first_stage <= stagenum & last_stage >= stagenum))

  ### Plot free energy scatterplots
  ###########################

  stagenum <- 5
  #All domains
  doubledeepms_free_energy_scatterplots(
    input_list = list(
      "GB1" = file.path(base_dir, paste0("002", "_doubledeepms_structure_metrics_GB1"), "dg_singles.txt"),
      "PSD95-PDZ3" = file.path(base_dir, paste0("002", "_doubledeepms_structure_metrics_PSD95-PDZ3"), "dg_singles.txt"),
      "GRB2-SH3" = file.path(base_dir, paste0("002", "_doubledeepms_structure_metrics_GRB2-SH3"), "dg_singles.txt")),
    input_MSA_list = list(
      "GB1" = file.path(base_dir, "Data", "MSA", "GB1", "frequencies.csv"),
      "PSD95-PDZ3" = file.path(base_dir, "Data", "MSA", "PSD95-PDZ3", "frequencies.csv"),
      "GRB2-SH3" = file.path(base_dir, "Data", "MSA", "GRB2-SH3", "frequencies.csv")
    ),
    outpath = doubledeepms__format_dir(dir_suffix="_doubledeepms_free_energy_scatterplots", stagenum=stagenum, base_dir=base_dir),
    colour_scheme = colour_scheme,
    execute = (first_stage <= stagenum & last_stage >= stagenum))

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
    execute = (first_stage <= stagenum & last_stage >= stagenum))
  #PSD95-PDZ3
  doubledeepms_free_energy_heatmaps(
    input_file = file.path(base_dir, paste0("002", "_doubledeepms_structure_metrics_PSD95-PDZ3"), "dg_singles.txt"),
    domain_name = "PSD95 PDZ3",
    outpath = doubledeepms__format_dir(dir_suffix="_doubledeepms_free_energy_heatmaps_PSD95-PDZ3", stagenum=stagenum, base_dir=base_dir),
    colour_scheme = colour_scheme,
    plot_width = 15,
    execute = (first_stage <= stagenum & last_stage >= stagenum))
  #GRB2-SH3
  doubledeepms_free_energy_heatmaps(
    input_file = file.path(base_dir, paste0("002", "_doubledeepms_structure_metrics_GRB2-SH3"), "dg_singles.txt"),
    domain_name = "GRB2 SH3",
    outpath = doubledeepms__format_dir(dir_suffix="_doubledeepms_free_energy_heatmaps_GRB2-SH3", stagenum=stagenum, base_dir=base_dir),
    colour_scheme = colour_scheme,
    plot_width = 11,
    execute = (first_stage <= stagenum & last_stage >= stagenum))

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
    input_MSA_list = list(
      "GB1" = file.path(base_dir, "Data", "MSA", "GB1", "frequencies.csv"),
      "PSD95-PDZ3" = file.path(base_dir, "Data", "MSA", "PSD95-PDZ3", "frequencies.csv"),
      "GRB2-SH3" = file.path(base_dir, "Data", "MSA", "GRB2-SH3", "frequencies.csv")
    ),
    outpath = doubledeepms__format_dir(dir_suffix="_doubledeepms_protein_stability_plots", stagenum=stagenum, base_dir=base_dir),
    colour_scheme = colour_scheme,
    execute = (first_stage <= stagenum & last_stage >= stagenum))

  ### Binding interface plots
  ###########################

  stagenum <- 8
  # GRB2-SH3
  doubledeepms_interface_mechanisms(
    base_dir = base_dir, 
    domain_name = "GRB2-SH3", 
    outpath = doubledeepms__format_dir(dir_suffix="_doubledeepms_ligand_interface_plots", stagenum=stagenum, base_dir=base_dir), 
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
    colour_scheme = colour_scheme,
    execute = (first_stage <= stagenum & last_stage >= stagenum))
  # All domains
  doubledeepms_binding_interface(
    input_list = list(
      "GB1" = file.path(base_dir, paste0("002", "_doubledeepms_structure_metrics_GB1"), "dg_singles.txt"),
      "PSD95-PDZ3" = file.path(base_dir, paste0("002", "_doubledeepms_structure_metrics_PSD95-PDZ3"), "dg_singles.txt"),
      "GRB2-SH3" = file.path(base_dir, paste0("002", "_doubledeepms_structure_metrics_GRB2-SH3"), "dg_singles.txt")),
    outpath = doubledeepms__format_dir(dir_suffix="_doubledeepms_binding_interface", stagenum=stagenum, base_dir=base_dir),
    colour_scheme = colour_scheme,
    execute = (first_stage <= stagenum & last_stage >= stagenum))
  

  ### Allostery plots
  ###########################

  stagenum <- 9
  #All domains
  doubledeepms_allostery_plots(
    input_file = file.path(base_dir, paste0("007", "_doubledeepms_protein_stability_plots"), "dg_singles.txt"),
    pdb_file_list = list(
      "GB1" = file.path(base_dir, "Data", "pdb", "1fcc.pdb"),
      "PSD95-PDZ3" = file.path(base_dir, "Data", "pdb", "1be9.pdb"),
      "GRB2-SH3" = file.path(base_dir, "Data", "pdb", "2vwf.pdb")),
    pdb_chain_query_list = list(
      "GB1" = "C",
      "PSD95-PDZ3" = "A",
      "GRB2-SH3" = "A"),
    annotation_list = list(
      "PSD95-PDZ3" = file.path(base_dir, "Data", "annotations", "PSD95-PDZ3", "PSD95-PDZ3_annotations.txt")),
    ohm_file_list = list(
      "GB1" = file.path(base_dir, "Data", "ohm", "ohm_2gb1_ddPCA_GB1only.txt"), 
      "PSD95-PDZ3" = file.path(base_dir, "Data", "ohm", "ohm_1be9_ddPCA_PDZonly.txt"),
      "GRB2-SH3" = file.path(base_dir, "Data", "ohm", "ohm_2vwf_ddPCA_GRB2only.txt")),
    outpath = doubledeepms__format_dir(dir_suffix="_doubledeepms_allostery_plots", stagenum=stagenum, base_dir=base_dir),
    colour_scheme = colour_scheme,
    execute = (first_stage <= stagenum & last_stage >= stagenum))

  ### Allostery scatterplots
  ###########################

  stagenum <- 10
  #All domains
  doubledeepms_allostery_scatterplots(
    input_file = file.path(base_dir, paste0("009", "_doubledeepms_allostery_plots"), "dg_singles.txt"),
    fitness_list = list(
      "PSD95-PDZ3" = file.path(base_dir, "Data", "fitness", "PSD95-PDZ3"),
      "GRB2-SH3" = file.path(base_dir, "Data", "fitness", "GRB2-SH3")),
    mochi_outpath_list = list(
      "PSD95-PDZ3" = file.path(base_dir, "Data", "mochi", "PSD95-PDZ3", "mochi__fit_tmodel_3state_sparse_dimsum128"),
      "GRB2-SH3" = file.path(base_dir, "Data", "mochi", "GRB2-SH3", "mochi__fit_tmodel_3state_sparse_dimsum128")),
    outpath = doubledeepms__format_dir(dir_suffix="_doubledeepms_allostery_scatterplots", stagenum=stagenum, base_dir=base_dir),
    colour_scheme = colour_scheme,
    execute = (first_stage <= stagenum & last_stage >= stagenum))

  ### Downsampling analysis
  ###########################

  stagenum <- 11
  #GB1
  doubledeepms_downsampling_analysis(
    mochi_outpath_prefix = file.path(base_dir, "Data", "mochi", "GB1", "mochi__fit_tmodel_3state_sparse_dimsum128_subsample"),
    sample_percentage = c(2, 5, 10, 20, 50, 100),
    literature_comparison_prefix = file.path(base_dir, paste0("001", "_doubledeepms_thermo_model_results_GB1subsample")),
    literature_comparison_list = list(
      "folding" = "validation_scatter_(mAvg)_mean_f_ddg_conf.txt"),
    outpath = doubledeepms__format_dir(dir_suffix="_doubledeepms_downsampling_analysis_GB1", stagenum=stagenum, base_dir=base_dir),
    colour_scheme = colour_scheme,
    execute = (first_stage <= stagenum & last_stage >= stagenum))
  #GRB2-SH3
  doubledeepms_downsampling_analysis(
    mochi_outpath_prefix = file.path(base_dir, "Data", "mochi", "GRB2-SH3", "mochi__fit_tmodel_3state_sparse_dimsum128_subsample"),
    sample_percentage = c(2, 5, 10, 20, 50, 100),
    outpath = doubledeepms__format_dir(dir_suffix="_doubledeepms_downsampling_analysis_GRB2-SH3", stagenum=stagenum, base_dir=base_dir),
    colour_scheme = colour_scheme,
    execute = (first_stage <= stagenum & last_stage >= stagenum))
  #PSD95-PDZ3
  doubledeepms_downsampling_analysis(
    mochi_outpath_prefix = file.path(base_dir, "Data", "mochi", "PSD95-PDZ3", "mochi__fit_tmodel_3state_sparse_dimsum128_subsample"),
    sample_percentage = c(2, 5, 10, 20, 50, 100),
    literature_comparison_prefix = file.path(base_dir, paste0("001", "_doubledeepms_thermo_model_results_PSD95-PDZ3subsample")),
    literature_comparison_list = list(
      "folding" = "validation_scatter_Calosci_2008_f_ddg_conf.txt",
      "binding" = "validation_scatter_Laursen_2020_b_ddg_conf.txt"),
    outpath = doubledeepms__format_dir(dir_suffix="_doubledeepms_downsampling_analysis_PSD95-PDZ3", stagenum=stagenum, base_dir=base_dir),
    colour_scheme = colour_scheme,
    execute = (first_stage <= stagenum & last_stage >= stagenum))

  ### Foldx comparisons
  ###########################

  stagenum <- 12
  #All domains
  doubledeepms_foldx_comparisons(
    input_file = file.path(base_dir, paste0("009", "_doubledeepms_allostery_plots"), "dg_singles.txt"),
    foldx_file_list = list(
      "GB1" = file.path(base_dir, "Data", "foldx", "GB1", "PS_2gb1_nowater_noligand_restrict_scanning_output.txt"), 
      "PSD95-PDZ3" = file.path(base_dir, "Data", "foldx", "PSD95-PDZ3", "PS_1be9_nowater_noligand_restrict_scanning_output.txt"),
      "GRB2-SH3" = file.path(base_dir, "Data", "foldx", "GRB2-SH3", "PS_2vwf_nowater_noligand_restrict_scanning_output.txt")),
    outpath = doubledeepms__format_dir(dir_suffix="_doubledeepms_foldx_comparisons", stagenum=stagenum, base_dir=base_dir),
    colour_scheme = colour_scheme,
    execute = (first_stage <= stagenum & last_stage >= stagenum))

  ### PolyPhen2 comparisons
  ###########################

  stagenum <- 13
  #All domains
  doubledeepms_polyphen2_comparisons(
    input_file = file.path(base_dir, paste0("009", "_doubledeepms_allostery_plots"), "dg_singles.txt"),
    polyphen2_file = file.path(base_dir, "Data", "polyphen2", "poyphen2_short.txt"),
    position_offset = list(
      # "GB1" = 0, 
      "PSD95-PDZ3" = 0,
      "GRB2-SH3" = 158),
    uniprot_id = list(
      # "GB1" = "", 
      "PSD95-PDZ3" = "DLG4_HUMAN",
      "GRB2-SH3" = "GRB2_HUMAN"),
    outpath = doubledeepms__format_dir(dir_suffix="_doubledeepms_polyphen2_comparisons", stagenum=stagenum, base_dir=base_dir),
    colour_scheme = colour_scheme,
    execute = (first_stage <= stagenum & last_stage >= stagenum))

  ### 3did comparisons
  ###########################

  stagenum <- 14
  #All domains
  doubledeepms_3did_comparisons(
    input_file = file.path(base_dir, paste0("009", "_doubledeepms_allostery_plots"), "dg_singles.txt"),
    threedid_file = file.path(base_dir, "Data", "3did", "3did_interface_flat"),
    alignment_file_list = list(
      "PSD95-PDZ3" = file.path(base_dir, "Data", "3did", "PSD95-PDZ3", "PF00595_seed_rat_human_align.fasta"),
      "GRB2-SH3" = file.path(base_dir, "Data", "3did", "GRB2-SH3", "PF00018_seed_chicken_human_align.fasta")),
    threedid_domain_name_list = list(
      "PSD95-PDZ3" = "PDZ",
      "GRB2-SH3" = "SH3_1"),
    position_offset_list = list(
      "PSD95-PDZ3" = 310,
      "GRB2-SH3" = 0),
    outpath = doubledeepms__format_dir(dir_suffix="_doubledeepms_3did_comparisons", stagenum=stagenum, base_dir=base_dir),
    colour_scheme = colour_scheme,
    execute = (first_stage <= stagenum & last_stage >= stagenum))

  ### EVE comparisons
  ###########################

  stagenum <- 15
  #All domains
  doubledeepms_eve_comparisons(
    eve_file_list = list(
      "GB1" = file.path(base_dir, "Data", "EVE", "GB1", "SPG1_STRSG_218-292_11-26-2021_b05_single_mutant_matrix.csv"),
      "PSD95-PDZ3" = file.path(base_dir, "Data", "EVE", "PSD95-PDZ3", "DLG4_HUMAN_301-404_11-26-2021_b04_single_mutant_matrix.csv"),
      "GRB2-SH3" = file.path(base_dir, "Data", "EVE", "GRB2-SH3", "GRB2_HUMAN_149-217_11-26-2021_b02_single_mutant_matrix.csv")),
    outpath = doubledeepms__format_dir(dir_suffix="_doubledeepms_eve_comparisons", stagenum=stagenum, base_dir=base_dir),
    colour_scheme = colour_scheme,
    execute = (first_stage <= stagenum & last_stage >= stagenum))

}



