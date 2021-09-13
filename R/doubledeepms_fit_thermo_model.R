
#' doubledeepms_fit_thermo_model
#'
#' Fit thermo models
#'
#' @param base_dir Base directory for all output files
#' @param mochi_output_name temperature in degrees celcuis (default:30)
#' @param tmodel_job_number Thermodynamic model fit job number: 1:final model, 2-11: monte carlo simluations for confidence intervals of model-inferred free energies (default:1)
#' @param tmodel_grid_search Thermodynamic model fit grid search to determine optimal hyperparameters (default:FALSE)
#' @param tmodel_protein Thermodynamic model fit proteins: comma-separated list of names or 'all' (default:'all')
#' @param tmodel_subset Thermodynamic model fit data subset: either an integer percentage of subsampled doubles (1-100) or "binding_only" or "singles_only" (default:"100")
#' @param tmodel_hyperparameters path table of hyperparameters
#' @param execute whether or not to execute the analysis (default: TRUE)
#'
#' @return Nothing
#' @export
#' @import data.table
doubledeepms_fit_thermo_model <- function(
  base_dir,
  mochi_output_name = "mochi__fit_tmodel_3state_sparse_dimsum128",
  tmodel_job_number = 1,
  tmodel_grid_search = FALSE,
  tmodel_protein = "all",
  tmodel_subset = "100",
  tmodel_hyperparameters,
  execute = TRUE
  ){

  #Return if analysis not executed
  if(!execute){
    return()
  }

  #Display status
  message(paste("\n\n*******", paste("running stage: doubledeepms_fit_thermo_model for", tmodel_protein), "*******\n\n"))

  #Create output directory
  mochi_base_dir <- file.path(base_dir, "Data", "mochi")
  doubledeepms__create_dir(doubledeepms_dir = mochi_base_dir, overwrite_dir = FALSE)

  #Domain name
  domain_names <- unlist(strsplit(tmodel_protein, ","))
  if(domain_names == "all"){ 
    domain_names = c("GB1", "GRB2-SH3", "PSD95-PDZ3")
  }

  #Loop over all proteins
  for(domain_name in domain_names){
    #Fitness data
    fitness_folding_inpath <- NULL
    fitness_binding_inpath <- NULL
    fitness_dir <- file.path(base_dir, "Data", "fitness", domain_name)
    if(!file.exists(fitness_dir)){
      stop(paste0("No fitness data for ", domain_name), call. = FALSE)
    }

    #Abundance fitness data (optional)
    pca_type <- "Abundance"
    fitness_dir <- file.path(base_dir, "Data", "fitness", domain_name, pca_type)
    if(file.exists(fitness_dir)){
      fitness_folding_inpath <- file.path(fitness_dir, list.files(fitness_dir))
    }

    #Binding fitness data (required)
    pca_type <- "Binding"
    fitness_dir <- file.path(base_dir, "Data", "fitness", domain_name, pca_type)
    if(file.exists(fitness_dir)){
      fitness_binding_inpath <- file.path(fitness_dir, list.files(fitness_dir))
    }else{
      stop(paste0("No binding fitness data for ", domain_name), call. = FALSE)
    }

    #Subset
    fit_type <- "fb"
    mut_order_subset_train <- "1,2"
    mut_order_subset_valid <- "2"
    mut_order_subset_obs <- "union"
    tmodel_subsample_prop <- 1
    if(!tmodel_subset %in% c(as.character(1:100), "binding_only", "singles_only")){
      stop(paste0("Invalid '", "tmodel_subset", "' argument (either an integer % of subsampled doubles (1-100) or 'binding_only' or 'singles_only')"), call. = FALSE)
    }
    if(tmodel_subset=="binding_only"){
      mochi_output_name <- paste0(mochi_output_name, "_bindingonly")
      fit_type <- "b"
    }else if(tmodel_subset=="singles_only"){
      mochi_output_name <- paste0(mochi_output_name, "_singlesonly")
      mut_order_subset_train <- "1"
      mut_order_subset_valid <- "1"
      mut_order_subset_obs <- "1,2"
    }else{
      tmodel_subsample_prop <- as.integer(tmodel_subset)/100
      mochi_output_name <- paste0(mochi_output_name, "_subsample", tmodel_subset, "p")
    }

    #Hyperparameters
    #Turn off grid search for monte carlo simulations
    if(tmodel_job_number!=1){
      tmodel_grid_search <- FALSE
    }
    learning_rate <- "0.001"
    if(tmodel_grid_search==TRUE){
      num_samples <- "128,256,512,1024"
    }else{
      hyper_dt <- fread(tmodel_hyperparameters)
      if(hyper_dt[protein==domain_name & subsample_prop==as.integer(tmodel_subsample_prop*100) & binding_only==(tmodel_subset=="binding_only") & singles_only==(tmodel_subset=="singles_only"),.N]!=1){
        stop(paste0("Hyperparameters not found. Rerun with tmodel_grid_search=TRUE and tmodel_job_number=1."), call. = FALSE)
      }
      num_samples <- hyper_dt[protein==domain_name & subsample_prop==as.integer(tmodel_subsample_prop*100) & binding_only==(tmodel_subset=="binding_only") & singles_only==(tmodel_subset=="singles_only"),as.character(num_samples)]
    }

    #Output folder
    doubledeepms__create_dir(doubledeepms_dir = file.path(mochi_base_dir, domain_name), overwrite_dir = FALSE)
    mochi_outpath <- file.path(mochi_base_dir, domain_name, mochi_output_name)

    #Run
    doubledeepms__mochi__fit_tmodel_3state(
      fitness_folding_inpath = fitness_folding_inpath,
      fitness_binding_inpath = fitness_binding_inpath,
      mut_order_subset_train = mut_order_subset_train,
      mut_order_subset_valid = mut_order_subset_valid,
      mut_order_subset_obs = mut_order_subset_obs,
      fit_type = fit_type,
      subsample_prop = tmodel_subsample_prop,
      mochi_script = system.file("python", "mochi__fit_tmodel_3state_doubledeepms.py", package = "doubledeepms"),
      mochi_outpath = mochi_outpath,
      num_samples = num_samples,
      learning_rate = learning_rate,
      job_number = tmodel_job_number)
  }
}
