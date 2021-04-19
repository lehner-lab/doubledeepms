
#' doubledeepms__format_dir
#'
#' Plot fitness versus growth rate.
#'
#' @param dir_suffix directory suffix (required)
#' @param stagenum stage number (required)
#' @param base_dir base directory (required)
#' @param no_digits number of digits to use for stage number
#'
#' @return Absolute directory character string
#' @export
doubledeepms__format_dir <- function(
  dir_suffix, 
  stagenum,
  base_dir,
  no_digits=3
  ){
  stage_id_zeroes <- paste0(rep("0", no_digits-nchar(as.character(stagenum))), collapse = "")
  stage_path <- file.path(base_dir, paste0(stage_id_zeroes, stagenum, dir_suffix))
  return(stage_path)
}
