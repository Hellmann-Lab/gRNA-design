# Libraries
library(tidyverse)
library(plyranges)
library(knitr)
library(kableExtra)
library(biomaRt)
library(GenomeInfoDb)
library(universalmotif)
library(wesanderson)
library(patchwork)
library(Gviz)
library(rtracklayer)
library(data.table)
library(ggtext)
library(lmerTest)
library(ggplotify)
library(gridExtra)
library(yaml)
library(reticulate)

######################################
## UDFs
######################################
annotate_nanopore_transcripts <- function(np, ref, set_seqlevel_style = "UCSC", 
                                          min_overlap = .6) {
  #' Annotates nanopore transcripts
  #'
  #'Uses a reference genome range to annotate nanopore transcripts. Can be used 
  #'for similar genome ranges.
  #'
  #'@param np Nanopore transcripts
  #'@param ref Reference genome range
  #'@param set_seqlevel_style seqlevelStyle for nanopore transcript
  #'@param min_overlap Minimum overlap of sequence between the reference and input transcripts. Values between 0 and 1.
  
  seqlevelsStyle(np) <- set_seqlevel_style
  seqlevelsStyle(ref) <- set_seqlevel_style
  
  ref <- ref %>% mutate(gencode_strand = strand) %>% 
    group_by(gene_name,gencode_strand) %>% 
    reduce_ranges()
  
  np_strand<- np %>% plyranges::filter(type == "exon") %>% as_tibble %>% 
    group_by(transcript_id) %>% 
    mutate( orig.size= sum(width)) %>% 
    as_granges() %>% 
    join_overlap_intersect( ref  ) %>% 
    as_tibble %>% 
    group_by( transcript_id, gene_name) %>% 
    mutate( bp_ol = sum(width),
            rel_ol =sum(width)/orig.size,
            strand = ifelse(sum(strand=="+")>sum(strand=="-"),"+","-") ) %>%
    group_by(transcript_id) %>% 
    filter( bp_ol == max(bp_ol) ) %>%
    filter( rel_ol >=min_overlap) %>% #, gencode_strand == strand ) %>% 
    transmute(strand=gencode_strand,transcript_id,gene_name, bp_ol,rel_ol)
  
  transcripts <- np %>% as_tibble %>%  
    filter(type == "mRNA") %>% 
    dplyr::select(-strand) %>% 
    inner_join(np_strand) %>%
    distinct_all %>% as_granges() 
  
  return(transcripts)
  
}

## Helper function - common TFs

calculate_motif_info_content<-function( pfm_file =jaspar_pfm_file, type = "jaspar"){
  #'Retrieves motif info from JASPAR or Homer file
  #'
  #'Retrieves motif info from JASPAR or Homer file
  #'
  #'@param pfm_file Position Frequency Matrix file
  #'@param type Type of PFM file
  
  if( type == "jaspar" ){
    x <- universalmotif::read_jaspar(pfm_file)
  }else if (type == "homer"){
    x<- universalmotif::read_homer(pfm_file)
  }else{
    return("type must be homer or jaspar\n")
  }
  
  tibble(motif_id = sapply( 1:length(x), function(i){ x[[i]]@name }),
         IC = sapply( 1:length(x), function(i){ x[[i]]@icscore }) )
}

# Get rank from tss_id
get_rank <- function(tss_id){
  #'Get rank from TSS_ID
  #'
  #'Splits the existing tss_id to retrieve rank as number (Internal Function)
  #'
  #'@param tss_id TSS ID
  return(as.integer(strsplit(tss_id, '_')[[1]][-1]))
}

# Select tss_id with minimum rank
tss_names_sort <- function(tss_ids){
  #'Handles multi-named tss_ids
  #'
  #'Any tss_id with multiple ids concatenated handled by this internal function (Internal Function)
  #'
  #'@param tss_ids TSS IDs
  list_ids <- strsplit(tss_ids, ',')[[1]]
  list_ranks <- sapply(list_ids, get_rank)
  
  return(list_ids[which.min(list_ranks)])
  
}

# distance to gene
TSS_dist2gene<- function(tss, gene_model, max.up = 1e5){
  #'Handles multi-named tss_ids
  #'
  #'Any transcript_id with multiple ids concatenated handled by this internal function (Internal Function)
  #'
  #'@param tss GRange object for TSSs
  #'@param gene_model Gene model prepared from GTF file
  #'@param max.up Maximum distance to consider upstream of gene
  
  gm_tss<- gene_model %>% 
    as_granges()  %>% 
    anchor_5p() %>% 
    mutate( gene_size = width,
            width =1)
  #could also be directed
  dd<-bind_ranges(tss %>%
                    as_granges() %>% 
                    join_nearest_upstream( gm_tss, suffix=c("_x",""), distance=T) %>% 
                    filter(gene_name_x == gene_name & distance < gene_size ) %>% 
                    plyranges::select(-gene_name_x) %>% 
                    mutate(dir = "downstream"), 
                  tss %>% 
                    as_granges() %>%
                    join_nearest_downstream(gm_tss ,suffix=c("_x",""),distance=T) %>% 
                    filter(gene_name_x == gene_name & distance < max.up) %>% 
                    plyranges::select(-gene_name_x)  %>% 
                    mutate(dir="upstream") ) %>% 
    as_tibble %>% 
    dplyr::group_by( transcript_id, gene_name ) %>% 
    dplyr::filter( distance == min(distance)) %>% 
    transmute( seqnames, start, end, strand, transcript_id, gene_name, dir, distance) %>% 
    as_granges
  
  return(dd)
  
}


# helper function: make gene models


# helper functions to cluster gRNAs based on position and keep the one with the highest on-target activity per cluster
bin_grnas_v1 <- function(grnas, binwidth) {
  #'Filters grnas based on predicted activity score
  #'
  #'Selects top 4 gRNAs for each gene. The sgRNAs with maximum predicted activity are selected
  #'
  #'@param grnas Table of selected grnas
  #'@param binwidth Extend the range by hlaf of binwidth a both ends 
  
  merged_grnas <- grnas %>%
    # convert to granges object
    dplyr::mutate(start = position + 1,
                  width = 1) %>%
    as_granges() %>%
    anchor_5p() %>%
    stretch(22) %>% 
    anchor_3p() %>%
    stretch(-3) %>%
    # get middle position of gRNA
    anchor_center() %>%
    stretch(-20) %>%
    # extend by the user-provided binwidth
    stretch(binwidth) %>%
    # gRNAs that overlap after extenion belong to the same cluster
    group_by(transcript) %>%
    reduce_ranges(sgID = paste(sgID, collapse = ","),
                  predicted_score = paste(predicted_score, collapse = ",")) %>%
    as_tibble() %>% 
    separate_rows(sgID, predicted_score, sep = ",") %>% 
    # find the gRNA with the highest on-target activity per cluster
    group_by(seqnames, start, end, transcript) %>% 
    dplyr::summarise(sgID = sgID[predicted_score == max(predicted_score)]) %>% 
    ungroup()
  
  # retrieve all information about the selected gRNAs originally provided in 'grnas' 
  binned_grnas <- grnas %>%
    inner_join(merged_grnas, by = c("seqnames", "transcript", "sgID")) %>%
    dplyr::select(-c("start", "end"))
  
  return(binned_grnas)
  
}

# helper functions to cluster gRNAs based on position and keep the one with the best combination of high on-target + low off-target activity per cluster
bin_grnas_v2 <- function(grnas, binwidth, max_reduction_activity) {
  #'Filters grnas based on predicted activity score and off-target effects
  #'
  #'Selects top 4 gRNAs for each gene. The sgRNAs with maximum predicted 
  #'activity and low off-target effects are selected. Maximum activity is 
  #'created after adjusting for off-target score 
  #'
  #'@param grnas Table of selected grnas
  #'@param binwidth Extend the range by hlaf of binwidth a both ends 
  #'@param max_reduction_activity Threshold for predicted activity compromised for off-target effects
  
  merged_grnas <- grnas %>%
    # convert to granges object
    dplyr::mutate(start = position + 1,
                  width = 1) %>%
    as_granges() %>%
    anchor_5p() %>%
    stretch(22) %>% 
    anchor_3p() %>%
    stretch(-3) %>%
    # get middle position of gRNA
    anchor_center() %>%
    stretch(-20) %>%
    # extend by the user-provided binwidth
    stretch(binwidth) %>%
    # gRNAs that overlap after extenion belong to the same cluster
    group_by(transcript) %>%
    reduce_ranges(sgID = paste(sgID, collapse = ","),
                  off_target_stringency = paste(off_target_stringency, collapse = ","),
                  predicted_score = paste(predicted_score, collapse = ",")) %>% 
    as_tibble() %>% 
    separate_rows(sgID, off_target_stringency, predicted_score, sep = ",") %>% 
    dplyr::mutate(off_target_stringency = as.integer(off_target_stringency),
                  predicted_score = as.double(predicted_score)) %>% 
    # heuristic to prefer gRNAs with high on-target activity scores AND low off-target scores
    group_by(seqnames, start, end, transcript) %>% 
    dplyr::arrange(desc(predicted_score), .by_group = T) %>%
    dplyr::mutate(predicted_score_best = max(predicted_score),
                  off_target_stringency_best = off_target_stringency[predicted_score == max(predicted_score)],
                  # if the highest scoring gRNA has a bad off-target score (>2) AND some other gRNAs have better off-target scores and not much worse on-target scores (how much worse is tolerated is defined by the user-provided parameter max_reduction activity), then instead of the original best gRNA we rather take the best out of these other gRNAs
                  predicted_score_chosen = ifelse(unique(off_target_stringency_best) > 2 & 
                                                    (sum((predicted_score_best - predicted_score < max_reduction_activity) + (off_target_stringency_best > off_target_stringency) == 2) > 0),
                                                  predicted_score[which((predicted_score_best - predicted_score < max_reduction_activity) + (off_target_stringency_best > off_target_stringency) == 2)[1]],
                                                  max(predicted_score))) %>% 
    dplyr::summarise(sgID = sgID[predicted_score == predicted_score_chosen]) %>% 
    ungroup()
  
  # retrieve all information about the selected gRNAs originally provided in 'grnas' 
  binned_grnas <- grnas %>%
    inner_join(merged_grnas, by = c("seqnames", "transcript", "sgID")) %>%
    dplyr::select(-c("start", "end"))
  
  return(binned_grnas)
  
}

# Violin plot to compare predicted score of selected vs all gRNAs
predicted_score_plot <- function(gRNAs_gg6_all, gRNAs_gg6_top4_v2, title = '') {
  #'Compare predicted score between two gRNA tables
  #'
  #'Violin plot to compare predicted score of two gRNA Tables
  #'
  #'@param gRNAs_gg6_all First table of selected grnas
  #'@param gRNAs_gg6_top4_v2 Second table of selected grnas
  #'@param title String for plot title
  bind_rows(
    species = bind_rows(
      `all designed` = gRNAs_gg6_all %>% dplyr::select(sgID, predicted_score),
      selected = gRNAs_gg6_top4_v2 %>% dplyr::select(sgID, predicted_score),
      .id = "gRNAs"
    ),
    .id = "species"
  ) %>%
    dplyr::mutate(
      gRNAs = factor(gRNAs, levels = c("all designed", "selected")),
      species = factor(species, levels = c("species")),
      category = paste0(species, "_", gRNAs)
    ) %>%
    ggplot(aes(x = gRNAs, y = predicted_score, fill = category)) +
    geom_violin(draw_quantiles = 0.5) +
    facet_wrap( ~ species) +
    theme_bw(base_size = 21) +
    theme(plot.title = element_text(hjust = 0.5))+
    scale_fill_manual(values = c("#FFD0D9", "#B0144F"), guide = "none") +
    ylab("predicted activity") + 
    ggtitle(title)
}

off_target_stringency_plot <- function(gRNAs_top4_v1, gRNAs_top4_v2, title = '') {
  #'Compare the stringency metrics
  #'
  #'Compare the stringency metrics between two gRNA tables
  #'
  #'@param gRNAs_top4_v1 First table of selected grnas
  #'@param gRNAs_top4_v2 Second table of selected grnas
  #'@param title String for plot title
  
  bind_rows(species = bind_rows(v1 = gRNAs_top4_v1,
                                v2 = gRNAs_top4_v2,
                                .id = "version"),
            .id = "species") %>% 
    dplyr::count(species, version, off_target_stringency) %>% 
    dplyr::mutate(species = factor(species, levels = c('species'))) %>% 
    ggplot(aes(x = version, y = n, fill = as.factor(off_target_stringency))) +
    geom_bar(stat = "identity", position = "stack") +
    theme_bw() +
    facet_wrap(~species) +
    scale_fill_manual(values = rev(wes_palette("Darjeeling1")), name = "off_target_stringency") + 
    ggtitle(title)+
    theme(plot.title = element_text(hjust = 0.5))
}

######################################
## Folder structure for Genome
######################################

create_folder_structure <- function(input_parameters){
  # Genome folder
  dir.create(paste0(input_parameters$working_dir, '/', input_parameters$genome))
  
  # Use genome folder
  genome_folder_path <- paste0(input_parameters$working_dir, '/', input_parameters$genome)
  
  # gRNA design folders
  dir.create(sprintf('%s/design_input_files', genome_folder_path))
  dir.create(sprintf('%s/design_input_files/bt_indexed_genomes_and_promoters', genome_folder_path))
  dir.create(sprintf('%s/design_output_files', genome_folder_path))
  dir.create(sprintf('%s/intermediate_files/', genome_folder_path))
  dir.create(paste0(genome_folder_path, '/intermediate_files/', input_parameters$genome ,'_bt_folder'))
  dir.create(sprintf('%s/shell_scripts', genome_folder_path))
  dir.create(sprintf('%s/SLURM_files', genome_folder_path))
  dir.create(sprintf('%s/RDS', genome_folder_path))
  dir.create(sprintf('%s/figures', genome_folder_path))
  dir.create(sprintf('%s/gRNAs', genome_folder_path))
}

######################################
## Seqnames Check
######################################

seqnames_check <- function(gtf_file, atac, np = NA){
  
  # Get seqnames of each input
  gtf_seqnames <- seqnames(gtf_file)
  atac_seqnames <- seqnames(atac)
  
  # If np is missing, ignore
  if (!is_empty(np)){
    np_seqnames <- seqnames(np)
  
  # Common seqnames
  common_seqnames <- c(Reduce(intersect, list(gtf_seqnames, np_seqnames, atac_seqnames)))
  
  # Check if required seqnames exist in all 3
  seqnames_missing <- list('gtf' = input_parameters$seqnames_to_keep[!(input_parameters$seqnames_to_keep %in% gtf_seqnames)], 
                           'np' = input_parameters$seqnames_to_keep[!(input_parameters$seqnames_to_keep %in% np_seqnames)], 
                           'atac' = input_parameters$seqnames_to_keep[!(input_parameters$seqnames_to_keep %in% atac_seqnames)])
  
  } else{
    # Common seqnames
    common_seqnames <- c(Reduce(intersect, list(gtf_seqnames, atac_seqnames)))
    
    # Check if required seqnames exist in atac and GTF
    seqnames_missing <- list('gtf' = input_parameters$seqnames_to_keep[!(input_parameters$seqnames_to_keep %in% gtf_seqnames)], 
                             'atac' = input_parameters$seqnames_to_keep[!(input_parameters$seqnames_to_keep %in% atac_seqnames)])
  }
  
  
  return(seqnames_missing)
  
}


######################################
# Extract TSS wrapper
######################################

extract_tss <- function(gtf_file, TF_list, gene_name = 'gene_name', transcript_type = 'protein_coding', type = 'exon', exon_number = 1) {
  #'Compare the stringency metrics
  #'
  #'Compare the stringency metrics between two gRNA tables
  #'
  #'@param gtf_file GTF File input
  #'@param TF_list Gene/Transcription factors for which exons and TSS are to be extracted
  #'@param gene_name Column name of gene is GTF file
  #'@param transcript_type GTF Transcript type plyranges::filter; eg., protein-encoding
  #'@param type GTF type: for plyranges::filter; here exon
  #'@param exon_number exon_number = 1 by default
  
  TF_gtf_exons <- gtf_file %>% 
    plyranges::filter(gene_name %in% TF_list &
                        transcript_type == transcript_type &
                        type == type)
  
  # TSS of the selected TFs
  TF_gtf_tss <- TF_gtf_exons %>% 
    filter(exon_number == exon_number) %>% 
    anchor_5p() %>% 
    mutate(width = 1)
  
  # Save exons and TSS from GTS
  saveRDS(TF_gtf_exons, sprintf("%s/RDS/TF_gtf_exons_%s.rds", input_parameters$working_dir, input_parameters$genome))
  saveRDS(TF_gtf_tss, sprintf("%s/RDS/TF_gtf_tss_%s.rds", input_parameters$working_dir, input_parameters$genome))
  
  return(list('TF_gtf_exons' = TF_gtf_exons %>% as_granges(), 'TF_gtf_tss' = TF_gtf_tss %>% as_granges()))
}

######################################
## Nanopore TSS generation wrapper
######################################
#!!! Generalise function 
generate_np_TSS <- function(np, TF_gtf_exons, input_parameters){
  
  if (input_parameters$np_expr_gff_path != ""){
    
  
  # Nanopore data annotation from GTF
  TF_np_gtf <- annotate_nanopore_transcripts(np, ref = TF_gtf_exons, min_overlap = 0.6, set_seqlevel_style = input_parameters$seq_level)
  
  # !!! Change string names to take in common string
  # !!! Does not make sense for just one GTF
  TF_np <-  as_tibble(TF_np_gtf) %>%
    dplyr::mutate(annot_evidence = case_when(!is.na(gene_name) ~ "gtf",
                                             T ~ "ambiguous"))
  
  # Save Nanopore based TF GRange with ambiguous evidence
  saveRDS(TF_np, sprintf("%s/RDS/TF_np_annot_%s_compare_evidence.rds", input_parameters$working_dir, input_parameters$genome))
  
  #Plot: Annotation evidence of nanopore transcripts
  png(sprintf("%s/figures/nanopore_annotation_evidence_%s.png", input_parameters$working_dir, input_parameters$genome) , width = 700, height = 450)
  TF_np %>% 
    dplyr::transmute(gtf = as.integer(grepl("gtf", annot_evidence)),
                     ambiguous = as.integer(annot_evidence == "ambiguous")) %>% 
    data.frame() %>% 
    UpSetR::upset(sets = c("ambiguous", "gtf"), order.by = "freq", keep.order = T, mainbar.y.label = "Annotation evidence of\nnanopore transcripts", text.scale = 1.5)
  dev.off()
  
  # Resolve ambiguous gene ID-gene name duplications
  TF_np <- TF_np %>% 
    group_by(gene_id) %>% 
    dplyr::filter(n_distinct(gene_name) == 1 | annot_evidence %in% c("gtf")) %>% 
    ungroup() %>% 
    dplyr::select(seqnames, start, end, width, strand, source, score, gene_id, transcript_id, gene_name, annot_evidence) %>% 
    as_granges()
  
  # Save Nanopore based TF GRange after filtering out ambiguous evidence
  saveRDS(TF_np, sprintf("%s/RDS/TF_np_annot_%s.rds", input_parameters$working_dir, input_parameters$genome))
  
  # Pull out nanopore-based TSS
  TF_np_tss <- TF_np %>% 
    as_granges() %>% 
    anchor_5p() %>% 
    mutate(width = 1)
  
  # Save Nanopore TSS
  saveRDS(TF_np_tss, sprintf("%s/RDS/TF_np_tss_annot_%s.rds", input_parameters$working_dir, input_parameters$genome))
  
  } else {
    TF_np_tss <- data.table()
  }
  
  return(TF_np_tss)
}

######################################
## Nanopore TSS checks
######################################
nanopore_tss_checks <- function(TF_np_tss, TF_list){
  #'Print checks for Nanopore TSS
  #'
  #'Checks for genes with multiple IDs, number of IDs mapped to multiple gene names,
  #'Also gives number of transcription factors in nanopore transcript 
  #'
  #'@param TF_np_tss Annotated Nanopore transcript TSSs

  # Number of genes with multiple gene IDs
  print(sprintf("Number of genes with multiple gene IDs: %i" , TF_np_tss %>% 
                  as_tibble() %>% 
                  group_by(gene_id) %>% 
                  mutate(n_gene_name_per_id = n_distinct(gene_name)) %>% 
                  dplyr::filter(n_gene_name_per_id > 1) %>% 
                  nrow()))
  
  #print("Number of genes with multiple gene IDs: ")
  print(TF_np_tss %>% 
          as_tibble() %>% 
          group_by(gene_id) %>% 
          summarise(n_gene_name_per_id = n_distinct(gene_name)) %>% 
          dplyr::count(n_gene_name_per_id))
  
  # Number of gene IDs with multiple gene names 
  print("Number of genes IDs with multiple gene names: ")
  print(TF_np_tss %>% 
          as_tibble() %>% 
          group_by(gene_name) %>% 
          summarise(n_gene_id_per_name = n_distinct(gene_id)) %>% 
          dplyr::count(n_gene_id_per_name))
  
  # How many of the TFs have at least 1 nanopore transcript?
  print("TFs have at least 1 nanopore transcript: ")
  print(table(TF_list %in% (TF_np_tss %>% as_tibble() %>% pull(gene_name) %>% unique())))
}


######################################
## Joining TSS Sources
######################################
combine_tss_sources <- function(input_parameters, TF_tss_all_sources, TF_np_tss = NA){
  
if (!is_empty(TF_np_tss)){
  
  # TSS from GTF
  TF_tss_all_sources <- as_tibble(TF_gtf_tss) %>% 
    dplyr::rename(transcript_id.gtf = transcript_id)
  
  # Distance of np TSS with all sources
  dist_tss_np_all <- as_tibble(TF_np_tss) %>%
    left_join(TF_tss_all_sources, by = c("gene_name", "seqnames", "strand"), suffix = c("", ".other"), relationship = "many-to-many") %>%
    rowwise() %>%
    dplyr::mutate(dist_np_other = ifelse(start >= start.other & start <= end.other,
                                         0, min(c(abs(start - start.other),
                                                  abs(start - end.other))))) %>%
    ungroup()

#Plot: Distribution of min distances per TSS and choose cutoff
p1 <- dist_tss_np_all %>%
  group_by(gene_name, start) %>%
  dplyr::filter(dist_np_other == min(dist_np_other)) %>%
  ungroup() %>%
  distinct(gene_name, dist_np_other) %>%
  ggplot(aes(x = dist_np_other)) +
  geom_density() +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_vline(xintercept = 100, linetype = "dashed", color = "red4") +
  xlab("distance between each nanopore TSS\nand the closest GTF TSS") +
  ggtitle("All")

p2 <- dist_tss_np_all %>%
  group_by(gene_name, start) %>%
  dplyr::filter(dist_np_other == min(dist_np_other)) %>%
  ungroup() %>%
  distinct(gene_name, dist_np_other) %>%
  dplyr::filter(dist_np_other < 500) %>%
  ggplot(aes(x = dist_np_other)) +
  geom_density() +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_vline(xintercept = 100, linetype = "dashed", color = "red4") +
  xlab("Distance between each nanopore TSS\nand the closest GTF TSS") +
  ggtitle("Distances below 500 bp")
p1 + p2

# Save plot
ggsave(sprintf("%s/figures/min_dist_np_other_tss_%s.png", input_parameters$working_dir, input_parameters$genome), width = 9, height = 4)

}
  else{
    
    # TSS from GTF
    TF_tss_all_sources <- as_tibble(TF_gtf_tss)
    
    # TSS from GTF only
    dist_tss_np_all <- TF_tss_all_sources %>% 
      dplyr::mutate(dist_np_other = 0)
  }
    
    return(list('TF_tss_all_sources' = TF_tss_all_sources, 'dist_tss_np_all' = dist_tss_np_all))
  
}

######################################
## Merged and Filtered TSS
######################################
create_TF_TSS <- function(input_parameters, dist_tss_np_all, atac){
  
  
if (input_parameters$np_expr_gff_path != ""){
  

# Only keep nanopore TSS that have a GTF TSS closer than 100bp 
TF_TSS <- dist_tss_np_all %>% 
  dplyr::filter(dist_np_other < input_parameters$max_dist_other_TSS) %>% 
  as_granges() %>% 
  # only keep nanopore TSS that overlap with an ATAC-seq peak
  join_overlap_intersect(atac %>% plyranges::select(atac_qValue = qValue, atac_signalValue = signalValue)) %>% 
  # merge TSS that are closer than stretch_for_TSS_merge
  anchor_center() %>% 
  stretch(input_parameters$stretch_for_TSS_merge) %>% 
  # as_tibble()
  group_by(gene_name) %>% 
  # keep track of the evidence, by counting up the number supporting nanopore mRNAs as well as the number of gencode transcript with support level 1/2 or manual annotation (HAVANA)
  reduce_ranges_directed(n_np = n_distinct(transcript_id[is.na(transcript_id)==F]),
                         n_gtf = n_distinct(na.omit(transcript_id.gtf[is.na(transcript_id.gtf)==F])),
                         np_transcripts = paste(unique(transcript_id[is.na(transcript_id)==F]), collapse=","),
                         gtf_transcripts = paste(unique(transcript_id.gtf[is.na(transcript_id.gtf)==F]), collapse=","),
                         atac_signalValue = max(atac_signalValue)) %>% 
  stretch(-input_parameters$stretch_for_TSS_merge) %>% 
  as_tibble() %>% 
  rowwise() %>% 
  # ad hoc evidence score to rank TSS
  mutate(evidence = sum(n_np*3 + n_gtf + atac_signalValue/1000)) %>% 
  group_by(gene_name) %>% 
  arrange(desc(evidence)) %>% 
  mutate(rnk = 1:n()) %>% 
  rowwise() %>% 
  mutate(tss_id = sprintf("%s_%s", gene_name, rnk),
         `TSS source` = paste(c("nanopore_atac", c("gtf")[c(n_gtf > 0)]), collapse = "_")) %>% 
  ungroup() 

} else{
  TF_TSS <- dist_tss_np_all %>% 
    dplyr::filter(dist_np_other < input_parameters$max_dist_other_TSS) %>% 
    as_granges() %>% 
    # only keep nanopore TSS that overlap with an ATAC-seq peak
    join_overlap_intersect(atac %>% plyranges::select(atac_qValue = qValue, atac_signalValue = signalValue)) %>% 
    # merge TSS that are closer than stretch_for_TSS_merge
    anchor_center() %>% 
    stretch(input_parameters$stretch_for_TSS_merge) %>% 
    # as_tibble()
    group_by(gene_name) %>% 
    # keep track of the evidence, by counting up the number supporting nanopore mRNAs as well as the number of gencode transcript with support level 1/2 or manual annotation (HAVANA)
    reduce_ranges_directed(n_gtf = n_distinct(transcript_id[is.na(transcript_id)==F]),
                           gtf_transcripts = paste(unique(transcript_id[is.na(transcript_id)==F]), collapse=","),
                           atac_signalValue = max(atac_signalValue)) %>% 
    stretch(-input_parameters$stretch_for_TSS_merge) %>% 
    as_tibble() %>% 
    rowwise() %>% 
    # ad hoc evidence score to rank TSS
    mutate(evidence = sum(n_gtf + atac_signalValue/1000)) %>% 
    group_by(gene_name) %>% 
    arrange(desc(evidence)) %>% 
    mutate(rnk = 1:n()) %>% 
    rowwise() %>% 
    mutate(tss_id = sprintf("%s_%s", gene_name, rnk),
           `TSS source` = paste(c("gtf")[c(n_gtf > 0)], collapse = "_")) %>% 
    ungroup()
}

  # Label genes with multiple TSSs in a new column - not removed
  TF_TSS_filt <- TF_TSS %>%
    as_tibble() %>%
    group_by(gene_name) %>%
    dplyr::mutate(evidence_ratio = evidence / max(evidence),
                  distance = start - min(start),
                  multi_tss = ifelse(length(unique(tss_id)) > 1, TRUE, FALSE)) %>%
    ungroup() %>%
    as_granges()
  
  # Save the list of filtered TSS
  saveRDS(TF_TSS_filt %>% as_granges(), sprintf("%s/RDS/TF_tss_%s.rds", input_parameters$working_dir, input_parameters$genome))

return(TF_TSS_filt)

}


######################################
## Merged and Filtered TSS: checks
######################################
filtered_TSS_checks <- function(input_parameters, TF_TSS){
  
# Number of unique TFs in selection
sprintf('Number of unique TFs in selection: %i', length(unique(TF_TSS$gene_name)))

# Number of TSSs per TF
print('Number of TSSs per TF: ')
print(TF_TSS %>% 
        dplyr::count(gene_name, name = "n_TSSs_per_TF") %>% 
        dplyr::count(n_TSSs_per_TF, name = "n_TF"))

#Plot: TSS width for merged and filtered TSS
ggplot(TF_TSS, aes(x = width)) +
  geom_histogram()+
  theme_bw() +
  xlab(sprintf("Width of merged & filtered TSSs in %s (bp)", input_parameters$genome))
ggsave(sprintf("%s/figures/TSS_width_%s.png", input_parameters$working_dir, input_parameters$genome), width = 5, height = 4)

# Different lines of TSS evidence
png(sprintf("figures/sources_of_TSSs_%s.png", input_parameters$genome), width = 700, height = 450)
TF_TSS %>% 
  dplyr::transmute(nanopore = 1,
                   ATAC = 1, 
                   GTF = n_gtf > 0) %>% 
  data.frame() %>% 
  UpSetR::upset(sets = c("nanopore", "ATAC", "GTF"), order.by = "freq", keep.order = T, mainbar.y.label = "Number of TSSs", text.scale = 1.5)
dev.off()

}


######################################
## gRNA Design tool inputs
######################################
# Combine all TSS and create bed file of tss regions
gRNA_design_TSS_regions <- function(gtf_file, input_parameters, type_sub = "transcript"){
  
  # GTF TSSs +-500 bp
  TSS_file <- gtf_file %>% 
    filter(type == type_sub) %>% 
    anchor_5p() %>% 
    mutate(width=1) %>%
    anchor_center() %>% 
    stretch(1000)
  
  # Change seqlevels based on input
  #seqlevelsStyle(TSS_file) <- input_parameters$seq_level
  
  return(TSS_file)
}

gRNA_design_tool_inputs <- function(input_parameters, TF_TSS_filt, gtf_file, np = NA){
  
  # Write TSSs to a tssTable.txt file (input format required by the gRNA design tool)
  TF_tssTable <- TF_TSS_filt %>%
    # move all coordinates right to convert to a zero-based coordinate system
    as_granges() %>% 
    shift_left(1) %>% 
    as_tibble %>% 
    # keep columns required by design tool
    transmute(gene = gene_name,
              transcripts = tss_id, 
              position = start, 
              strand = strand, 
              chromosome = seqnames,
              `TSS source` = TSS.source,
              `primary TSS` = paste0("(",start,", ",end,")"),
              #`secondary TSS` = paste0("(",start,", ",end,")"))
              `cage peak ranges` = "[]")
  
  # Save Table
  write_delim(TF_tssTable, sprintf("%s/design_input_files/%s_tssTable.txt", input_parameters$working_dir, input_parameters$genome), delim = "\t", quote = "none")

  # Define X and MT chromosomes
  MT_chr <- if_else(input_parameters$seq_level == 'UCSC', 'chrM', 'M')
  X_chr <- if_else(input_parameters$seq_level == 'UCSC', 'chrX', 'X')
  
  # Bed files
  gtf_tss_regions <- gRNA_design_TSS_regions(gtf_file, input_parameters, type_sub = "transcript")
  
  # Check if nanopore/expression file exists
  if (!is_empty(np)){
    # Get Nanopore TSS regions
    np_tss_regions <- gRNA_design_TSS_regions(np, input_parameters, type_sub = "mRNA")
    
    # Sort out chromosome names
    seq_levels_to_keep <- c(Reduce(intersect, list(seqlevels(gtf_tss_regions), seqlevels(np_tss_regions))))
    
    # Remove MT and X chromosomes
    seq_levels_to_keep <- seq_levels_to_keep[!seq_levels_to_keep %in% c(X_chr, MT_chr)]
    
    #!!! only for two inputs right now: will need to be a loop/function for 2+ GTFs
    gtf_tss_regions <- keepSeqlevels(gtf_tss_regions, seq_levels_to_keep, pruning.mode = "coarse")
    np_tss_regions <- keepSeqlevels(np_tss_regions, seq_levels_to_keep, pruning.mode = "coarse")
    
    # combine GTF and nanopore
    all_tss_regions <- bind_ranges(gtf_tss_regions, np_tss_regions) %>% 
      reduce_ranges(gene = paste(unique(gene_name),sep="_")) %>% 
      mutate(start = ifelse(start < 0, 1, start),
             tss_id = paste0("tss_", seqnames, "_", start))
    names(all_tss_regions) <- all_tss_regions$tss_id
  }
  
  else{
    # combine GTF and nanopore
    all_tss_regions <- gtf_tss_regions %>% 
      reduce_ranges(gene = paste(unique(gene_name),sep="_")) %>% 
      mutate(start = ifelse(start < 0, 1, start),
             tss_id = paste0("tss_", seqnames, "_", start))
    names(all_tss_regions) <- all_tss_regions$tss_id
}
  
  # Write bed file
  write_bed(all_tss_regions, sprintf("%s/design_input_files/bt_indexed_genomes_and_promoters/all_tss_regions_noMT_%s.bed", input_parameters$working_dir, input_parameters$genome))
}

######################################
## Shell script for bowtie_indices
######################################

bowtie_indices_script <- function(input_parameters){
  
  # Get genome file folder
  #genome_folder <- dirname(input_parameters$genome_fa)
  
  # Add script for filtering noMT genome file  
  if (!file.exists(sprintf('%s/design_input_files/bt_indexed_genomes_and_promoters/genome_%s_noMT.fa', input_parameters$working_dir, input_parameters$genome))){
    
      #!!! Check chromosome list
    noMTfasta_op <- c('#Get genome FASTA (without the MT chromosome)',
      sprintf("samtools faidx %s %s > %s/design_input_files/bt_indexed_genomes_and_promoters/genome_%s_noMT.fa \n", input_parameters$genome_fa, paste(input_parameters$seqnames_to_keep, collapse = " "), input_parameters$working_dir, input_parameters$genome),
      
      '#Get genome FASTA (without the MT chromosome)',
      sprintf("sbatch --wait --job-name=bowtie_index_genome --wrap='bowtie-build %s/design_input_files/bt_indexed_genomes_and_promoters/genome_%s_noMT.fa %s/design_input_files/bt_indexed_genomes_and_promoters/genome_%s_noMT' \n", input_parameters$working_dir, input_parameters$genome, input_parameters$working_dir, input_parameters$genome))
  } else{
    noMTfasta_op <- '#Fasta with noMT exisits'
  }
    
  # # Add script for 2bit genome file
  # Create filepath for 2bit genome file
  twobit_path <- paste0(input_parameters$working_dir, '/design_input_files/bt_indexed_genomes_and_promoters/', input_parameters$genome, '.2bit')
  
  if ((input_parameters$genome_2bit == '') & (!file.exists(twobit_path))){
    
    twobit_op <- c('#Create Genome 2bit',
    sprintf("faToTwoBit %s %s \n", input_parameters$genome_fa, twobit_path))
    
  } else if (file.exists(twobit_path)){
    twobit_op <- sprintf('#2Bit genome file exists at %s', twobit_path)
    input_parameters$genome_2bit = twobit_path
  }
    
    else {
      twobit_op <- sprintf('#2Bit genome file exists at %s', input_parameters$genome_2bit)
    }
 
  
  # Create Shell script for bowtie_indices
  writeLines(c('#!/bin/bash', 
               paste0('#SBATCH --error=',input_parameters$working_dir,'/SLURM_files/bowtie_indices.%J.err'),
               paste0('#SBATCH --output=',input_parameters$working_dir,'/SLURM_files/bowtie_indices.%J.out'),
               '#SBATCH --cpus-per-task=10',
               '#SBATCH --mem=30G',
               
               '#Change working directory',           
               sprintf("cd %s/design_input_files/bt_indexed_genomes_and_promoters/ \n", input_parameters$working_dir),
               
               # Whether filtered fasta script is needed
               noMTfasta_op,
               
               # Whether twobit genomescript is needed
               twobit_op,
               
               '#2bit to FASTA for selected TFs',
               sprintf("/opt/bin/twoBitToFa -bed=all_tss_regions_noMT_%s.bed %s all_tss_regions_noMT_%s.fa \n", input_parameters$genome, input_parameters$genome_2bit, input_parameters$genome),
               
               '#Bowtie Indexing',
               sprintf("sbatch --wait --job-name=bowtie_index_promoters_%s --wrap='bowtie-build all_tss_regions_noMT_%s.fa all_tss_regions_noMT_%s' \n", input_parameters$genome, input_parameters$genome, input_parameters$genome)), 
             
             # Saving to shell scripts
             sprintf("%s/shell_scripts/create_bowtie_indices.sh", input_parameters$working_dir))
  system(sprintf('chmod u+r+x %s/shell_scripts/create_bowtie_indices.sh', input_parameters$working_dir))
  
}

######################################
## Shell script for gRNA Design
######################################
gRNA_design_script <- function(input_parameters){
  
# Create Shell script for gRNA Design
writeLines(c('#!/bin/bash', 
             paste0('#SBATCH --error=',input_parameters$working_dir,'/SLURM_files/gRNA_design.%J.err'),
             paste0('#SBATCH --output=',input_parameters$working_dir,'/SLURM_files/gRNA_design.%J.out'),
             '#SBATCH --cpus-per-task=10',
             '#SBATCH --mem=30G \n',
             "cd /opt/miniconda3/bin/",
             "source activate",
             "conda activate /data/home/termeg/.conda/envs/CRISPRiaDesign", 
             'export PATH="$PATH:/opt/bin" \n',
             sprintf("cd %s \n", input_parameters$working_dir),
             sprintf("sbatch --wait --job-name=design_gRNAs_workflow --wrap='
python %s --tss %s/design_input_files/%s_tssTable.txt \ --model %s \ --genome %s \ --atac %s \ --btp  %s/design_input_files/bt_indexed_genomes_and_promoters/all_tss_regions_noMT_%s \ --btg  %s/design_input_files/bt_indexed_genomes_and_promoters/genome_%s_noMT \ --btf  intermediate_files/%s_bt_folder \ -o design_output_files/TF_gRNAs_%s -n %i'", input_parameters$gRNA_function_path, input_parameters$working_dir, input_parameters$genome, input_parameters$gRNA_model_path, input_parameters$genome_fa, input_parameters$atac_bw_path, input_parameters$working_dir, input_parameters$genome, input_parameters$working_dir, input_parameters$genome, input_parameters$genome, input_parameters$genome, input_parameters$grna_num),
"conda deactivate"), 
"./shell_scripts/design_grnas.sh", sep ='\n') 
system(sprintf('chmod u+r+x %s/shell_scripts/design_grnas.sh', input_parameters$working_dir))
}

######################################
## gRNA Design results - wrapper
######################################
grna_design_results <- function(input_parameters)
{
  # top_n_gRNAs
  gRNAs <- data.table::fread(sprintf("%s/design_output_files/TF_gRNAs_%s_top%s.csv", input_parameters$working_dir, input_parameters$genome, input_parameters$grna_num)) %>% 
    group_by(gRNA_sequence) %>%
    dplyr::filter(length(unique(gene)) == 1) %>% #Remove gRNAs that target multiple genes
    ungroup()
  
  ## take the gRNA with the highest on-target activity per bin
  gRNAs_top4_v1 <- gRNAs %>% 
    dplyr::rename(seqnames = chromosome) %>% 
    dplyr::mutate(predicted_score = as.double(format(predicted_score, digits = 6))) %>% 
    bin_grnas_v1(input_parameters$binwidth) %>% 
    group_by(transcript) %>% 
    slice_max(order_by = predicted_score, n = 4) %>% 
    ungroup()
  
  ## take the gRNA with the best combination of high on-target and low off-target activity per bin
  gRNAs_top4_v2 <- gRNAs %>% 
    dplyr::rename(seqnames = chromosome) %>% 
    dplyr::mutate(predicted_score = as.double(format(predicted_score, digits = 6))) %>% 
    bin_grnas_v2(input_parameters$binwidth, input_parameters$max_reduction_activity)  %>% 
    group_by(transcript) %>% 
    slice_max(order_by = predicted_score, n = 4) %>% 
    ungroup() %>% 
    mutate(binwidth = input_parameters$binwidth)
  
  # Reduce binwidth parameter if not all transcripts have 4 gRNAs
  criteria <- gRNAs_top4_v2  %>% group_by(transcript) %>% summarise(num = length(sgID))  %>% filter(num < 4)
  if (nrow(criteria) > 0){
    # Get trannscripts for which the criteria will be relaxed
   relaxed_criteria_transcripts <- criteria %>% pull(transcript)
   rel_cat_gRNAs <- gRNAs %>% filter(transcript %in% relaxed_criteria_transcripts)
   
   # Drop the transcripts for final result
   gRNAs_top4_v2 <- gRNAs_top4_v2 %>% 
     filter(!(transcript %in% relaxed_criteria_transcripts))
   
   # Relax binwidth
   relaxed_binwidth = input_parameters$binwidth - 3
   
   # Re-run with relaxed parameter:
   while (relaxed_binwidth > 0){
     
   # Get top gRNAs for subset transcripts
   relaxed_gRNAs_top4_v2 <- rel_cat_gRNAs %>% 
     dplyr::rename(seqnames = chromosome) %>% 
     dplyr::mutate(predicted_score = as.double(format(predicted_score, digits = 6))) %>%
     bin_grnas_v2(relaxed_binwidth, input_parameters$max_reduction_activity)  %>% 
     group_by(transcript) %>% 
     slice_max(order_by = predicted_score, n = 4) %>% 
     ungroup()
   
   # Filter for transcripts that did not satisfy the criteria
   relaxed_criteria_transcripts <- relaxed_gRNAs_top4_v2  %>% group_by(transcript) %>% summarise(num = length(sgID))  %>% filter(num < 4) %>% pull(transcript)
   
   # Join back the transcripts that did satisfy the criteria to final result
   satisfied_gRNAs_top_v2 <- relaxed_gRNAs_top4_v2 %>% 
     filter(!(transcript %in% relaxed_criteria_transcripts)) %>% 
     mutate(binwidth = relaxed_binwidth) # Note the value of relaxed binwidth
   
   gRNAs_top4_v2 <- rbind(gRNAs_top4_v2, satisfied_gRNAs_top_v2) #rbindto final results
   
   # Remove the selected gRNAs from the next iteration
   rel_cat_gRNAs <- gRNAs %>% filter(transcript %in% relaxed_criteria_transcripts)
   
   # Further reduce binwidth
   relaxed_binwidth <- relaxed_binwidth - 3
   }
   
   # If any transcripts remain, add them with binwidth = 0
   relaxed_gRNAs_top4_v2$binwidth <- 0
   gRNAs_top4_v2 <- rbind(gRNAs_top4_v2, relaxed_gRNAs_top4_v2)
   
  }
  
  # Save to RDS 
  saveRDS(gRNAs_top4_v2, sprintf("%s/RDS/gRNAs_%s_top4.rds", input_parameters$working_dir, input_parameters$genome))
  
  # load all gRNAs
  gRNAs_all <- data.table::fread(sprintf("%s/design_output_files/TF_gRNAs_%s_unfiltered.csv", input_parameters$working_dir, input_parameters$genome), skip = 2, col.names = c("sgID", "gene", "transcript", "gRNA_sequence", "genomic_sequence", "predicted_score", "39_nearTSS", "31_nearTSS", "21_genome", "31_2_nearTSS", "31_3_nearTSS", "seqnames", "position", "strand"), select = c(1:11, 14, 18, 19)) %>% 
    dplyr::mutate(off_target_stringency = case_when(`31_nearTSS` & `21_genome` ~ 0,
                                                    `31_nearTSS` ~ 1,
                                                    `21_genome` ~ 2,
                                                    `31_2_nearTSS` ~ 3,
                                                    `31_3_nearTSS` ~ 4),
                  is_in_final_lib = sgID %in% gRNAs_top4_v2$sgID) %>% 
    dplyr::select(-c( "39_nearTSS", "31_nearTSS", "21_genome", "31_2_nearTSS", "31_3_nearTSS")) %>% 
    drop_na(off_target_stringency) 
  saveRDS(gRNAs_all, sprintf("%s/RDS/gRNAs_%s_all.rds", input_parameters$working_dir, input_parameters$genome))
  
  gRNA_all_gr <- gRNAs_all %>% 
    dplyr::mutate(start = position + 1,
                  width = 1) %>% 
    as_granges() %>% 
    anchor_5p() %>%
    stretch(22) %>% 
    anchor_3p() %>%
    stretch(-3)
  #saveRDS(gRNA_all_gr, "shiny_app/data_files/gorGor6/gRNAs/gRNAs.rds")
  saveRDS(gRNA_all_gr, sprintf("%s/gRNAs/gRNAs.rds", input_parameters$working_dir))
  
  # Plots
  pred_plot <- predicted_score_plot(gRNAs_all, gRNAs_top4_v2, title = input_parameters$genome)
  ggsave(filename = sprintf("%s/figures/predicted_score_plot_%s.png", input_parameters$working_dir, input_parameters$genome), plot = pred_plot, width = 9, height = 4)
  
  off_target_plot <- off_target_stringency_plot(gRNAs_top4_v1, gRNAs_top4_v2, title = input_parameters$genome)
  ggsave(filename = sprintf("%s/figures/off_target_stringency_plot_%s.png", input_parameters$working_dir, input_parameters$genome), plot = off_target_plot, width = 9, height = 4)
  
  # pull out off-target scores
  off_target_scores <- gRNAs %>%
    dplyr::count(transcript, off_target_stringency)

  # # check genes with the worst off-target scores
  # off_target_scores %>%
  #   pivot_wider(names_from = "off_target_stringency", names_prefix = "off_target_", values_from = "n", values_fill = 0) %>%
  #   arrange(desc(off_target_4)) %>%
  #   head()

  # plot off-target scores
  p1 <- off_target_scores %>%
    group_by(transcript) %>%
    dplyr::mutate(off_target1234 = sum(n[off_target_stringency %in% 1:4])) %>%
    ungroup() %>%
    arrange(desc(off_target1234)) %>%
    dplyr::mutate(transcript = factor(transcript, levels = transcript[off_target_stringency == 0]),
                  off_target_stringency = as.factor(off_target_stringency)) %>%
    ggplot(aes(x = transcript, y = n, fill = off_target_stringency)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = rev(wes_palette("Darjeeling1"))) +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    ggtitle(input_parameters$genome) +
    ylab("number of gRNAs")


  p1 + plot_layout(guides = "collect")
  ggsave(sprintf("%s/figures/off_target_scores_per_gene_plot_%s.png", input_parameters$working_dir, input_parameters$genome), width = 14, height = 5)
}

######################################
## App Input data preparation wrapper
######################################

create_app_data_files <- function(input_parameters){
  # Set working directory
  workdir <- sprintf("%s/data_files", input_parameters$shiny_working_dir)
  if (!dir.exists(workdir)){
    dir.create(workdir)
  }
  setwd(workdir)
  
  # Change working directory to genome folder
  #input_parameters$working_dir <- paste0(input_parameters$working_dir, '/', input_parameters$genome)
  
  # Source helper functions
  source(input_parameters$shiny_helper_path)
  
  # Create data list
  data4shiny <- list('organism' = list(genome = input_parameters$genome,
                                       gtf = list(Liftoff = input_parameters$gtf_file_path),
                                       nanopore_reads =  input_parameters$np_expr_gff_path,
                                       nanopore_transcripts = sprintf("%s/RDS/TF_np_annot_%s.rds", input_parameters$working_dir, input_parameters$genome),
                                       atac = list('organism' = input_parameters$atac_bw_path),
                                       tss = list("final" = sprintf("%s/RDS/TF_tss_%s.rds", input_parameters$working_dir, input_parameters$genome)),
                                       grnas = sprintf("%s/gRNAs/gRNAs.rds", input_parameters$working_dir), gene_models = c()))
  
  # Change name to genome name
  names(data4shiny) <- input_parameters$genome
  names(data4shiny[[input_parameters$genome]]$atac) <- input_parameters$genome
  
  # Create species folder
  sapply(names(data4shiny), function(i) {system(paste0("mkdir ", workdir, '/', data4shiny[[i]][["genome"]]))})
  
  # Change TF_list to include TF from Gorilla and Orang
  TF_list <- readRDS(input_parameters$TF_list_path)
  
  # Data setup for app
  # Creating Gene models
  set_up_data_for_shiny <- function(species, data4shiny, TF_list, workdir) {
    
    # GTF
    for (gtf_type in names(data4shiny[[species]]$gtf)) {
      
      # Make gene_models folder
      system(paste("mkdir",  glue::glue("{workdir}/{data4shiny[[species]]$genome}/gene_models/")))
      
      # Create gene models for all GTF files
      if (gtf_type %in% c("GENCODE", "Liftoff")) {
        
        gene_models <- makeGviz_geneModel_from_gencode(gtf = read_gff(data4shiny[[species]]$gtf[[gtf_type]]) %>% filter(gene_name %in% TF_list))
        saveRDS(gene_models, paste0(workdir,'/', data4shiny[[species]]$genome, "/gene_models/", gtf_type, ".rds"))
        
      }
      
      if (gtf_type == "ENSEMBL") {
        
        gene_models <- makeGviz_geneModel_from_ensembl(gtf = read_gff(data4shiny[[species]]$gtf[[gtf_type]]) %>% filter(gene_name %in% TF_list))
        saveRDS(gene_models, paste0(workdir,'/', data4shiny[[species]]$genome, "/gene_models/", gtf_type, ".rds"))
        
      }
      
    }
    
    # Add nanopore expression data to app if it exists
    if (input_parameters$np_expr_gff_path != ''){
      
      # Make nanopore_reads folder
      system(paste("mkdir",  glue::glue("{workdir}/{data4shiny[[species]]$genome}/nanopore_reads/")))
      
      # nanopore reads
      system(glue::glue("ln -s {data4shiny[[species]]$nanopore_reads} {workdir}/{data4shiny[[species]]$genome}/nanopore_reads/nanopore.bam"))
      system(glue::glue( "ln -s {data4shiny[[species]]$nanopore_reads}.bai {workdir}/{data4shiny[[species]]$genome}/nanopore_reads/nanopore.bam.bai" ))
      
      # nanopore transcripts
      gene_models <- makeGviz_geneModel_from_np(gtf = readRDS(data4shiny[[species]]$nanopore_transcripts) %>% filter(gene_name %in% TF_list))
      saveRDS(gene_models, paste0(workdir, '/', data4shiny[[species]]$genome, "/gene_models/Nanopore.rds"))
    }
    
    # ATAC-seq data
    atac_dir <- glue::glue("{workdir}/{data4shiny[[species]]$genome}/atac/")
    system(paste("mkdir", atac_dir))
    
    i = names(data4shiny[[species]]$atac)
    j = data4shiny[[species]]$atac[[i]]  
    system(glue::glue("mkdir {atac_dir}{i}"))
    
    n=0
    nfn <- glue::glue("{atac_dir}{i}/{i}_{n}.bw")
    system(glue::glue("ln -s {j} {nfn}"))
    
    # TSS
    for (tss_type in names(data4shiny[[species]]$tss)) {
      
      # Make tss folder
      system(paste("mkdir",  glue::glue("{workdir}/{data4shiny[[species]]$genome}/tss/")))
      
      system(glue::glue("ln -s {data4shiny[[species]]$tss[[tss_type]]} {workdir}/{data4shiny[[species]]$genome}/tss/{tss_type}.rds"))
      
    }
    
    # gRNAs
    
    # Create gRNAs folder
    system(paste("mkdir", glue::glue("{workdir}/{data4shiny[[species]]$genome}/gRNAs/")))
    system(glue::glue("ln -s {data4shiny[[species]]$grnas} {workdir}/{data4shiny[[species]]$genome}/gRNAs/gRNAs.rds"))
    
  }
  
  # Run the data creation function for required species
  set_up_data_for_shiny(input_parameters$genome, data4shiny, TF_list, workdir)
}

######################################
## END
######################################
