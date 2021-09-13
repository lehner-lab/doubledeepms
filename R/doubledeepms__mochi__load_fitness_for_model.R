
#' doubledeepms__mochi__load_fitness_for_model
#'
#' Load fitness to be used for model fitting.
#'
#' @param dimsum_RData_file DiMSum .RData file (required)
#' @param order_subset Subset of mutation orders (default:0:2)
#' @param sequence_type Sequence type: 'nucleotide' or 'aminoacid' (default:aminoacid)
#' @param mean_input_count_threshold Mean input count threshold (default:0)
#'
#' @return Data.table of fitness
#' @export
#' @import data.table
doubledeepms__mochi__load_fitness_for_model <- function(
  dimsum_RData_file,
  order_subset = 0:2,
  sequence_type = "aminoacid",
  mean_input_count_threshold = 0
  ){

  #Load and encode
  load(dimsum_RData_file)

  #Subset variants
  if(sequence_type == "aminoacid"){
    all_variants <- all_variants[STOP==F & STOP_readthrough==F & Nham_aa %in% order_subset & Nham_aa==Nmut_codons & mean_count>=mean_input_count_threshold & !is.na(fitness)]
  }
  if(sequence_type == "nucleotide"){
    all_variants <- all_variants[Nham_nt %in% order_subset & mean_count>=mean_input_count_threshold & !is.na(fitness)]
  }
  fitness_dt <- doubledeepms__mochi__encode_genotypes(all_variants, sequenceType = sequence_type, remove_constant = F)[["all"]]

  #Add columns
  fitness_dt[, aa_seq := all_variants[, aa_seq]]
  fitness_dt[, mean_count := all_variants[, mean_count]]
  fitness_dt[, Nham_nt := all_variants[, Nham_nt]]
  fitness_dt[, Nham_aa := all_variants[, Nham_aa]]

  return(fitness_dt)
}
