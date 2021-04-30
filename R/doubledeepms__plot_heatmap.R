
#' doubledeepms__plot_heatmap
#'
#' Plot free energy heatmap.
#'
#' @param input_dt input data.table (required)
#' @param variable_name variable name for heatmap cells (required)
#' @param metric_name input data.table (required)
#' @param output_file plot output file path (required)
#' @param width plot width in inches (default:10)
#' @param height plot height in inches (default:4)
#' @param plot_title main title for plot (default:empty)
#' @param colour_clip maximum absolute value of colour scale (default:4)
#' @param colour_low colour scale lower limit colour (default:blue)
#' @param colour_high colour scale upper limit colour (default:red)
#' @param colour_mid colour scale zero colour (default:lightgrey)
#' @param colour_mid_text text mid colour (default:grey)
#'
#' @return Nothing
#' @export
#' @import data.table
doubledeepms__plot_heatmap <- function(
  input_dt,
  variable_name,
  metric_name,
  output_file,
  width = 10,
  height = 4,
  plot_title = "",
  colour_clip = 4,
  colour_low = "blue",
  colour_high = "red",
  colour_mid = "lightgrey",
  colour_mid_text = "grey"
  ){

  #All amino acids
  aa_obj <- Biostrings::AAString("GAVLMIFYWKRHDESTCNQP")
  aa_list <- Biostrings::AMINO_ACID_CODE[strsplit(as.character(aa_obj), NULL)[[1]]]
  # aa_list["*"] <- "X"

  #WT sequence
  wt_seq <- unlist(unique(input_dt[order(input_dt[,"Pos_ref"]),c("WT_AA", "Pos_ref")])[,"WT_AA"])

  #Construct heatmap matrix
  heat_mat <- matrix(nrow = length(aa_list), ncol = input_dt[,length(unique(Pos_ref))])
  rownames(heat_mat) <- names(aa_list)
  colnames(heat_mat) <- input_dt[order(Pos_ref),unique(Pos_ref)]
  for(aa_pos in input_dt[order(Pos_ref),unique(Pos_ref)]){
    for(aa_id in names(aa_list)){
      temp_index <- which(input_dt[,Pos_ref]==aa_pos & input_dt[,Mut]==aa_id)
      if(length(temp_index)==1){
        heat_mat[aa_id,as.character(aa_pos)] <- input_dt[temp_index,.SD[[1]],,.SDcols = variable_name]
      }
    }
  }

  #Add distance and metric
  misc_mat <- matrix(c(
    # input_dt[order(Pos_ref)][!duplicated(Pos_ref),colour_clip*.SD[[1]]/max(.SD[[1]]),,.SDcols = metric_name],
    input_dt[order(Pos_ref)][!duplicated(Pos_ref),colour_clip-2*colour_clip*scHAmin_ligand/max(scHAmin_ligand)]), nrow = 1, byrow = T)
  # rownames(misc_mat) <- c("Effect", "Proximity")
  rownames(misc_mat) <- c("Prox.")
  colnames(misc_mat) <- colnames(heat_mat)

  #Merge
  heat_mat <- rbind(heat_mat, misc_mat)

  #PLot low confidence estimates as points
  heat_mat_point <- matrix(F, nrow = dim(heat_mat)[1], ncol = dim(heat_mat)[2])
  rownames(heat_mat_point) <- rownames(heat_mat)
  colnames(heat_mat_point) <- colnames(heat_mat)
  for(aa_pos in input_dt[order(Pos_ref),unique(Pos_ref)]){
    for(aa_id in names(aa_list)){
      temp_index <- which(input_dt[,Pos_ref==aa_pos & Mut==aa_id & .SD[[1]]==F,,.SDcols = c(paste0(variable_name, "_conf"), variable_name)])
      if(length(temp_index)==1){
        heat_mat_point[aa_id,as.character(aa_pos)] <- T
      }
    }
  }
  if(sum(heat_mat_point)==0){
    heat_mat_point <- NULL
  }

  #Add WT sequence to column names
  colnames(heat_mat) <- paste0(colnames(heat_mat), "\n", wt_seq)

  #Matrix text - add binding interface residues
  input_matrix_text <- matrix("", nrow = dim(heat_mat)[1], ncol = dim(heat_mat)[2])
  bi_residues <- which(input_dt[order(Pos_ref)][!duplicated(Pos_ref),Pos_class=="binding_interface"])
  input_matrix_text[dim(input_matrix_text)[1],bi_residues] <- "B"
  rownames(input_matrix_text) <- rownames(heat_mat)
  colnames(input_matrix_text) <- colnames(heat_mat)

  #Plot
  doubledeepms__tile_heatmap_wrapper(
    heat_mat, 
    text_size = 2.5, 
    output_file = output_file, 
    colour_low = colour_low, 
    colour_high = colour_high, 
    colour_mid = colour_mid, 
    xlab = plot_title, 
    ylab = "Mutant AA", 
    cluster =  "none", 
    width = width, 
    height = height, 
    xaxis_angle = 0, 
    xaxis_hjust = 0.5, 
    xaxis_size = 5, 
    na_colour = "white",
    colour_clip = colour_clip,
    input_matrix_point = heat_mat_point,
    input_matrix_text = input_matrix_text)
}
