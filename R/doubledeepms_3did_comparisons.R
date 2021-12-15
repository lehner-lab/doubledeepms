
#' doubledeepms_3did_comparisons
#'
#' Plot free energy heatmaps.
#'
#' @param input_file path to input file (required)
#' @param threedid_file 3did flat file (required)
#' @param alignment_file_list alignment file list (required)
#' @param threedid_domain_name_list 3did domain name list (required)
#' @param position_offset_list position offset list (required)
#' @param outpath output path for plots and saved objects (required)
#' @param colour_scheme colour scheme file (required)
#' @param execute whether or not to execute the analysis (default: TRUE)
#'
#' @return Nothing
#' @export
#' @import data.table
doubledeepms_3did_comparisons <- function(
  input_file,
  threedid_file,
  alignment_file_list,
  threedid_domain_name_list,
  position_offset_list,
  outpath,
  colour_scheme,
  execute = TRUE
  ){

  #Return if analysis not executed
  if(!execute){
    return()
  }

  #Display status
  message(paste("\n\n*******", "running stage: doubledeepms_3did_comparisons", "*******\n\n"))

  #Create output directory
  doubledeepms__create_dir(doubledeepms_dir = outpath)

  ### Load single mutant free energies
  ###########################

  #Load dg data
  dg_dt <- fread(input_file)

  ### Load 3did file
  ###########################

  #Load dg data
  tdid_dt <- fread(threedid_file, header = F, fill = T)
  #rows of interest
  idrows <- which(tdid_dt[,1]=="#=ID")
  idrows_oi <- which(tdid_dt[,V1=="#=ID" & V2 %in% unlist(threedid_domain_name_list)])
  #Subset to rows of interest
  for(i in idrows_oi){
    tdid_dt[(i+1):(idrows[which(idrows==i)+1]-1), did1 := unlist(tdid_dt[i,2])] 
    tdid_dt[(i+1):(idrows[which(idrows==i)+1]-1), did2 := unlist(tdid_dt[i,3])] 
  }
  #Reformat
  tdid_dt <- tdid_dt[!is.na(did1) & V1!="//"]
  tdid_dt[, V2 := sapply(strsplit(V2, ": "), '[', 2)]
  #Format positions to comma separated list instead of ranges
  for(i in 1:tdid_dt[,.N]){
    tdid_dt[i, Pos_3did := paste0(sapply(strsplit(tdid_dt[i,V2], " ")[[1]], function(x){paste0(unlist(strsplit(x, "-"))[1]:unlist(strsplit(x, "-"))[2], collapse=",")}), collapse = ",")]
  }

  ### Add 3did data to variants 
  ###########################

  tdid_list <- list()
  for(i in names(alignment_file_list)){
    #Convert positions to ddPCA residues positions
    #Alignment positions data.table
    astring_dt <- data.table(
      alignment = unlist(strsplit(unlist(fread(alignment_file_list[[i]], header = F)[2,1]), "")))
    astring_dt[alignment!="-", Pos_alignment := 1:.N]
    #Mapping between alignment and motif/3did positions
    apos_dt <- fread(gsub(".fasta", "_pos.txt", alignment_file_list[[i]]))
    apos_list <- as.list(unlist(apos_dt[,Pos_3did]))
    names(apos_list) <- unlist(apos_dt[,Pos_alignment])
    #Add to data.table
    astring_dt[as.character(Pos_alignment) %in% names(apos_list), Pos_3did := unlist(apos_list[as.character(Pos_alignment)])]
    #Conservation of domain residues
    aseed_dt <- as.data.table(do.call("rbind", strsplit(unlist(fread(gsub("_seed.*", "_seed.txt", alignment_file_list[[i]]), header = F)[,2]), "")))
    aseed_cons <- aseed_dt[,apply(.SD, 2, function(x){sum(x==names(table(x))[which(table(x)==max(table(x)))[1]])})/.N]
    #Add to data.table
    astring_dt[alignment!="-", conservation := aseed_cons]
    #Remove gaps
    astring_dt <- astring_dt[alignment!="."]
    #Add ddPCA positions
    astring_dt[, Pos_ref := 1:.N + position_offset_list[[i]]]
    #Translate 3did interaction positions to ddPCA positions
    pos_dict <- as.list(unlist(astring_dt[!is.na(Pos_3did),Pos_ref]))
    names(pos_dict) <- astring_dt[!is.na(Pos_3did),Pos_3did]
    tdid_dt[did1==threedid_domain_name_list[[i]], Pos_ref := sapply(strsplit(Pos_3did, ","), function(x){paste0(pos_dict[as.character(x)], collapse=",")})]
    #Add 3did data to ddG data.table (ignore ligands)
    tdid_list[[i]] <- unlist(strsplit(tdid_dt[did1==threedid_domain_name_list[[i]] & !grepl(paste0(threedid_domain_name_list[[i]], "_LIG_"), did2),Pos_ref], ","))
    tdid_list[[i]] <- as.numeric(tdid_list[[i]][tdid_list[[i]]!="NULL"])
    dg_dt[protein==i, threedid_intsum := sum(tdid_list[[i]]==Pos_ref),Pos_ref]
    dg_dt[protein==i & Pos_ref %in% astring_dt[is.na(Pos_3did),Pos_ref], threedid_intsum := NA]
    #Add 3did data to ddG data.table (ignore ligands) - homotypic interactions only
    tdid_list[[i]] <- unlist(strsplit(tdid_dt[did1==threedid_domain_name_list[[i]] & did2==threedid_domain_name_list[[i]],Pos_ref], ","))
    tdid_list[[i]] <- as.numeric(tdid_list[[i]][tdid_list[[i]]!="NULL"])
    dg_dt[protein==i, threedid_intsum_homo := sum(tdid_list[[i]]==Pos_ref),Pos_ref]
    dg_dt[protein==i & Pos_ref %in% astring_dt[is.na(Pos_3did),Pos_ref], threedid_intsum_homo := NA]
    #Add 3did conservation data to ddG data.table
    cons_list <- as.list(unlist(astring_dt[!is.na(conservation),conservation]))
    names(cons_list) <- unlist(astring_dt[!is.na(conservation),Pos_ref])
    dg_dt[protein==i & as.character(Pos_ref) %in% names(cons_list), domain_conservation := unlist(cons_list[as.character(Pos_ref)])]
  }

  ### Enrichment of allosteric mutations at interaction interfaces
  ###########################

  plot_dt <- dg_dt[!is.na(allosteric_mutation) & !is.na(threedid_intsum),.(prop_mut = sum(allosteric_mutation)/length(allosteric_mutation)*100, scHAmin_ligand = unique(scHAmin_ligand), Pos_class=unique(Pos_class), allosteric=unique(allosteric), threedid_intsum = unique(threedid_intsum), domain_conservation = unique(domain_conservation)),.(Pos_ref, protein)]
  cor_dt <- plot_dt[,.(cor = round(cor(threedid_intsum, prop_mut, use = "pairwise.complete", method = "spearman"), 2), threedid_intsum = Inf, prop_mut = Inf),protein]
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(threedid_intsum, prop_mut)) +
    ggplot2::geom_smooth(method = "lm", se = F, color = "grey", linetype = 2, formula = 'y ~ x') + 
    ggplot2::geom_point(ggplot2::aes(shape = !is.na(allosteric), color = Pos_class, size = domain_conservation)) +
    ggplot2::xlab(expression("Number of distinct interaction interfaces for homologs (3did)")) +
    ggplot2::ylab("%Allosteric mutations per residue") +
    ggplot2::facet_wrap(protein~., scales = "free", ncol = 2) +
    ggplot2::labs(color = "Residue\nposition", shape = "Major allosteric\nsite", size = "Residue\nconservation") +
    ggrepel::geom_text_repel(ggplot2::aes(label = Pos_ref, color = Pos_class), show.legend = F, max.overlaps = 5) +
    # ggplot2::geom_text(ggplot2::aes(label = Pos_ref), colour = "black", size = 1) +
    ggplot2::geom_text(data = cor_dt, ggplot2::aes(label=paste("rho = ", cor, sep="")), hjust = 1, vjust = 1) +
    ggplot2::theme_bw()
  if(!is.null(colour_scheme)){
    d <- d + ggplot2::scale_colour_manual(values = unlist(colour_scheme[["shade 0"]][c(3, 4)]))
  }
  ggplot2::ggsave(file.path(outpath, "allosteric_mutations_vs_3did.pdf"), d, width = 8, height = 4, useDingbats=FALSE)

  ### Enrichment of allosteric mutations at interaction interfaces - homotypic
  ###########################

  plot_dt <- dg_dt[!is.na(allosteric_mutation) & !is.na(threedid_intsum_homo),.(prop_mut = sum(allosteric_mutation)/length(allosteric_mutation)*100, scHAmin_ligand = unique(scHAmin_ligand), Pos_class=unique(Pos_class), allosteric=unique(allosteric), threedid_intsum_homo = unique(threedid_intsum_homo), domain_conservation = unique(domain_conservation)),.(Pos_ref, protein)]
  cor_dt <- plot_dt[,.(cor = round(cor(threedid_intsum_homo, prop_mut, use = "pairwise.complete", method = "spearman"), 2), threedid_intsum_homo = Inf, prop_mut = Inf),protein]
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(threedid_intsum_homo, prop_mut)) +
    ggplot2::geom_smooth(method = "lm", se = F, color = "grey", linetype = 2, formula = 'y ~ x') + 
    ggplot2::geom_point(ggplot2::aes(shape = !is.na(allosteric), color = Pos_class, size = domain_conservation)) +
    ggplot2::xlab(expression("Number of distinct interaction interfaces for homologs (3did)")) +
    ggplot2::ylab("%Allosteric mutations per residue") +
    ggplot2::facet_wrap(protein~., scales = "free", ncol = 2) +
    ggplot2::labs(color = "Residue\nposition", shape = "Major allosteric\nsite", size = "Residue\nconservation") +
    ggrepel::geom_text_repel(ggplot2::aes(label = Pos_ref, color = Pos_class), show.legend = F, max.overlaps = 5) +
    # ggplot2::geom_text(ggplot2::aes(label = Pos_ref), colour = "black", size = 1) +
    ggplot2::geom_text(data = cor_dt, ggplot2::aes(label=paste("rho = ", cor, sep="")), hjust = 1, vjust = 1) +
    ggplot2::theme_bw()
  if(!is.null(colour_scheme)){
    d <- d + ggplot2::scale_colour_manual(values = unlist(colour_scheme[["shade 0"]][c(3, 4)]))
  }
  ggplot2::ggsave(file.path(outpath, "allosteric_mutations_vs_3did_homo.pdf"), d, width = 8, height = 4, useDingbats=FALSE)

}

