
#' doubledeepms__mochi__fit_tmodel_3state
#'
#' Fit 3-state thermodynamic model to ddPCA data.
#'
#' @param fitness_folding_inpath Folding fitness DiMSum RData file (default:NULL)
#' @param fitness_binding_inpath Folding fitness DiMSum RData file (default:NULL)
#' @param id_subset_folding_inpath Folding coefficient id subset file (default:all)
#' @param id_subset_binding_inpath Folding coefficient id subset file (default:all)
#' @param mut_order_subset_train Mutation order subset for training data, in addition to WT (default:1,2)
#' @param mut_order_subset_valid Mutation order subset for validation data, in addition to WT (default:2)
#' @param mut_order_subset_obs Mutation order subset for all observed data, in addition to WT (default:union of mutOrderSubsetTrain and mutOrderSubsetValid)
#' @param min_multiplicity_valid Minimum multiplicity for validation data (default:2)
#' @param fit_type Fit type, either fb=folding+binding, f=folding, b=binding (default:fb)
#' @param sequence_type Sequence type (default:aminoacid)
#' @param subsample_prop Doubles subsample proportion (default:1)
#' @param mochi_script MoCHI script (required)
#' @param mochi_outpath Output folder (required)
#' @param num_epochs Number of epochs (default:1000)
#' @param num_epochs_grid Number of epochs for gridsearch (default:100)
#' @param validation_set_proportion Validation set proportion (default:0.1)
#' @param num_resamplings Number of random resamplings of target values (default:10)
#' @param num_samples Number of samples per gradient update (default:128,256,512,1024)
#' @param learning_rate Learning rate (default:0.0001,0.001,0.01,0.1)
#' @param job_number Job number (default:1)
#'
#' @return Nothing
#' @export
#' @import data.table
doubledeepms__mochi__fit_tmodel_3state <- function(
  fitness_folding_inpath=NULL,
  fitness_binding_inpath=NULL,
  id_subset_folding_inpath="all",
  id_subset_binding_inpath="all",
  mut_order_subset_train="1,2",
  mut_order_subset_valid="2",
  mut_order_subset_obs="union",
  min_multiplicity_valid=2,
  fit_type="fb",
  sequence_type="aminoacid",
  subsample_prop=1,
  mochi_script,
  mochi_outpath,
  num_epochs=1000,
  num_epochs_grid=100,
  validation_set_proportion=0.1,
  num_resamplings=10,
  num_samples="128,256,512,1024",
  learning_rate="0.0001,0.001,0.01,0.1",
  job_number=1
  ){

  ###########################
  ### Globals
  ###########################

  job_number <- job_number-1

  ###########################
  ### Load data
  ###########################

  #Create output directory
  dir.create(mochi_outpath, recursive = T)

  #Organise mutation order subsets
  mut_order_subset_valid <- unique(c(0, as.integer(unlist(strsplit(mut_order_subset_valid, ",")))))
  #Training mutation order subset must include validation mutation order subset
  mut_order_subset_train <- unique(c(0, mut_order_subset_valid, as.integer(unlist(strsplit(mut_order_subset_train, ",")))))
  #Observed mutation order subset must include validation and training mutation order subsets
  if(mut_order_subset_obs!="union"){
    mut_order_subset_obs <- unique(c(0, as.integer(unlist(strsplit(mut_order_subset_obs, ",")))))
  }else{
    mut_order_subset_obs <- c(0)
  }
  mut_order_subset_obs <- unique(c(mut_order_subset_obs, mut_order_subset_train, mut_order_subset_valid))

  #Load coefficient subsets
  id_subset_folding <- NULL
  id_subset_binding <- NULL
  if(id_subset_folding_inpath!="all"){id_subset_folding <- fread(id_subset_folding_inpath)[,id]}
  if(id_subset_binding_inpath!="all"){id_subset_binding <- fread(id_subset_binding_inpath)[,id]}

  #Load folding data
  fitness_folding_dt <- data.table()
  if(!is.null(fitness_folding_inpath)){
    fitness_folding_dt <- doubledeepms__mochi__load_fitness_for_model(
      dimsum_RData_file = fitness_folding_inpath, 
      order_subset = mut_order_subset_obs,
      sequence_type = sequence_type)
    fitness_folding_dt[, dataset_folding := 1]
    fitness_folding_dt[, dataset_binding := 0]
  }

  #Load binding data
  fitness_binding_dt <- data.table()
  if(!is.null(fitness_binding_inpath)){
    fitness_binding_dt <- doubledeepms__mochi__load_fitness_for_model(
      dimsum_RData_file = fitness_binding_inpath,
      order_subset = mut_order_subset_obs,
      sequence_type = sequence_type)
    fitness_binding_dt[, dataset_folding := 0]
    fitness_binding_dt[, dataset_binding := 1]
  }

  #Merge datasets
  fitness_dt <- rbind(fitness_folding_dt, fitness_binding_dt, fill = T)

  #Term order
  fitness_dt[, mut_order := sapply(strsplit(id, ","), length)]
  fitness_dt[id=="-0-", mut_order := 0]

  #Subsample doubles
  if(subsample_prop!=1){
    fitness_dt[mut_order==2, subsample := ceiling((1:.N)/(.N*subsample_prop))]
    set.seed(1)
    fitness_dt[mut_order==2, subsample := sample(subsample, .N, replace = F)]
    fitness_dt <- fitness_dt[which(mut_order!=2 | subsample==1),.SD,,.SDcols = names(fitness_dt)[names(fitness_dt)!="subsample"]]
  }

  #Remove constant regions
  fitness_dt[, var := doubledeepms__mochi__remove_constant_positions(var)]

  #Mutation ids
  fitness_dt[, mut1 := sapply(strsplit(id, ","), '[', 1)]
  fitness_dt[, mut2 := sapply(strsplit(id, ","), '[', 2)]

  #Folding and binding multiplicity
  all_mult <- table(unlist(fitness_dt[id!="-0-",.(mut1, mut2)]))
  fold_mult <- table(unlist(fitness_dt[dataset_folding==1 & id!="-0-",.(mut1, mut2)]))
  bind_mult <- table(unlist(fitness_dt[dataset_binding==1 & id!="-0-",.(mut1, mut2)]))
  all_mult_dict <- as.list(as.integer(all_mult))
  names(all_mult_dict) <- names(all_mult)
  fold_mult_dict <- as.list(as.integer(fold_mult))
  names(fold_mult_dict) <- names(fold_mult)
  bind_mult_dict <- as.list(as.integer(bind_mult))
  names(bind_mult_dict) <- names(bind_mult)

  ###########################
  ### Coefficients that can be fit
  ###########################

  #Coefficients that can be fit
  fold_coef_permitted <- names(all_mult_dict)[unlist(all_mult_dict)>0]
  bind_coef_permitted <- names(bind_mult_dict)[unlist(bind_mult_dict)>0]

  #Subset data to that associated with coefficients that can be fit
  fitness_dt <- fitness_dt[id=="-0-" | dataset_binding==1 | (dataset_folding==1 & mut1 %in% fold_coef_permitted & (is.na(mut2) | mut2 %in% fold_coef_permitted))]
  fitness_dt <- fitness_dt[id=="-0-" | dataset_folding==1 | (dataset_binding==1 & mut1 %in% bind_coef_permitted & (is.na(mut2) | mut2 %in% bind_coef_permitted))]

  #Folding coefficients for model
  if(is.null(id_subset_folding)){
    id_subset_folding <- fold_coef_permitted
  }else{
    id_subset_folding <- id_subset_folding[id_subset_folding %in% fold_coef_permitted]
  }

  #Binding coefficients for model
  if(is.null(id_subset_binding)){
    id_subset_binding <- bind_coef_permitted
  }else{
    id_subset_binding <- id_subset_binding[id_subset_binding %in% bind_coef_permitted]
  }

  ###########################
  ### Data permitted in validation set 
  ###########################

  #Mutations that are permitted in validation set (WT and expendible doubles)
  fold_val_permitted <- names(all_mult_dict)[unlist(all_mult_dict)>=min_multiplicity_valid]
  bind_val_permitted <- names(bind_mult_dict)[unlist(bind_mult_dict)>=min_multiplicity_valid]
  fitness_dt[, validation_permitted := F]
  # fitness_dt[dataset_folding==1 & mut_order==2, validation_permitted := mut1 %in% fold_val_permitted & mut2 %in% fold_val_permitted]
  # fitness_dt[dataset_binding==1 & mut_order==2, validation_permitted := mut1 %in% bind_val_permitted & mut2 %in% bind_val_permitted]
  fitness_dt[dataset_folding==1 & mut_order %in% mut_order_subset_valid & fit_type %in% c("f", "fb"), validation_permitted := mut_order==0 | (mut1 %in% fold_val_permitted & mut2 %in% fold_val_permitted)]
  fitness_dt[dataset_binding==1 & mut_order %in% mut_order_subset_valid & fit_type %in% c("b", "fb"), validation_permitted := mut_order==0 | (mut1 %in% bind_val_permitted & mut2 %in% bind_val_permitted)]

  ###########################
  ### Duplicate WT
  ###########################

  #Duplicate WT to represent 1% of training data
  fitness_dt <- rbind(fitness_dt[rep(which(id=="-0-" & dataset_folding==1), fitness_dt[dataset_folding==1,as.integer(.N/100)]),], fitness_dt[id!="-0-" | dataset_folding!=1])
  fitness_dt <- rbind(fitness_dt[rep(which(id=="-0-" & dataset_binding==1), fitness_dt[dataset_binding==1,as.integer(.N/100)]),], fitness_dt[id!="-0-" | dataset_binding!=1])

  ###########################
  ### Assign training/validation data
  ###########################

  if(validation_set_proportion!=0){
    #Assign training/validation data
    fitness_dt[validation_permitted==T, training_fold := ceiling((1:.N)/(.N*validation_set_proportion))]
    #Split all validation data into 10 equally sized folds
    if(job_number==0){
      #Choose orthogonal random validation set
      set.seed(2)
    }else{
      set.seed(1)
    }
    fitness_dt[validation_permitted==T, training_fold := sample(training_fold, .N, replace = F)]
    fitness_dt[validation_permitted==T, training_set := training_fold]
    fitness_dt[training_set==job_number | (job_number==0 & training_set==1), training_set := 0]
    #Training data or not in this validation fold
    fitness_dt[(dataset_folding==1 & mut_order %in% mut_order_subset_train & fit_type %in% c("f", "fb") & validation_permitted==F) | training_set!=0, training_set := 1]
    fitness_dt[(dataset_binding==1 & mut_order %in% mut_order_subset_train & fit_type %in% c("b", "fb") & validation_permitted==F) | training_set!=0, training_set := 1]
  }else{
    #Training data
    fitness_dt[dataset_folding==1 & mut_order %in% mut_order_subset_train & fit_type %in% c("f", "fb"), training_set := 1]
    fitness_dt[dataset_binding==1 & mut_order %in% mut_order_subset_train & fit_type %in% c("b", "fb"), training_set := 1]
  }

  ###########################
  ### Shuffle
  ###########################

  #Shuffle data
  set.seed(100)
  fitness_dt <- fitness_dt[sample(1:.N, .N, replace = F)]

  ###########################
  ### Save data
  ###########################

  #Run name
  run_name <- ""
  if(job_number!=0){
    run_name <- file.path("bootstrap", paste0("boot_", job_number))
  }

  #Output directory
  dir.create(file.path(mochi_outpath, run_name), recursive = T)

  #Save model data
  if(job_number==0){
    write.table(fitness_dt, file = file.path(mochi_outpath, run_name, "model_data.txt"), sep = "\t", row.names = F, col.names = T, quote = F)
  }

  #Feature matrix
  input_matrix_binding <- doubledeepms__mochi__construct_feature_matrix(fitness_dt, prefix = "bind", mut_ids = id_subset_binding)
  input_matrix_folding <- doubledeepms__mochi__construct_feature_matrix(fitness_dt, prefix = "fold", mut_ids = id_subset_folding)

  #Selection columns
  selection_matrix <- fitness_dt[,.(dataset_folding, dataset_binding)]

  #Target variable, target variable SE and training_set 
  target_matrix <- fitness_dt[,.(fitness, fitness_sd = SE, training_set, variant_sequence = var)]

  #All input data
  input_matrix <- cbind(selection_matrix, input_matrix_binding, input_matrix_folding[,.SD,,.SDcols = !grepl("WT", names(input_matrix_folding))], target_matrix)

  #Save - all
  write.table(input_matrix, file = file.path(mochi_outpath, run_name, "dataset_all.txt"), sep = ",", row.names = F, quote = F)
  #Save - train
  write.table(input_matrix[training_set==1], file = file.path(mochi_outpath, run_name, "dataset_train.txt"), sep = ",", row.names = F, quote = F)
  #Save - validate
  write.table(input_matrix[training_set==0 | (training_set==1 & validation_set_proportion==0)], file = file.path(mochi_outpath, run_name, "dataset_valid.txt"), sep = ",", row.names = F, quote = F)

  ###########################
  ### Delete objects to free up memory
  ###########################

  rm(list = ls()[grep("matrix|fitness", ls())])
  gc()

  ###########################
  ### Fit model
  ###########################

  #Run mochi python script on command-line
  system(paste0(
    "python ", mochi_script,
    " --data_train ",
    file.path(mochi_outpath, run_name, "dataset_train.txt"),
    " --data_valid ",
    file.path(mochi_outpath, run_name, "dataset_valid.txt"),
    " --data_obs ",
    file.path(mochi_outpath, run_name, "dataset_all.txt"),
    " -o ",
    file.path(mochi_outpath, run_name),
    " -e ",
    num_epochs_grid,
    " --l1_regularization_factor ",
    "0",
    " --l2_regularization_factor ",
    "0",
    " -p ",
    num_epochs,
    " --num_samples ",
    num_samples,
    " --learning_rate ",
    learning_rate,
    " --num_resamplings ",
    num_resamplings,
    " --num_models ",
    1,
    " --random_seed ",
    job_number+1))

}

