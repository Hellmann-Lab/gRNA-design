library(tidyverse)
library(topGO)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(biomaRt)
library(yaml)

# Set working directory
setwd('/data/share/htp/perturb-seq/gRNA_design_workflow/')

# Read user paramters
input_parameters <- read_yaml('./shiny_app/input.yaml') 

# UDFs
pull_expr_percent <- function(cntrl_cnts, TF_list, percent_expr) {
  #'Expression (%) data table
  #'
  #'Returns a data table with expression percent for each TF in TF_list 
  #'
  #'@param cntrl_cnts Expression counts for all TFs
  #'@param TF_list List of Transcription factors
  #'@param percent_expr Expression threshold
  expr_human_cntrl <- data.frame(gene_name = rownames(cntrl_cnts),
                                 percent_expr = rowSums(cntrl_cnts > 0) / ncol(cntrl_cnts) * 100) %>% 
    right_join(data.frame(gene_name = TF_list)) %>% 
    dplyr::mutate(percent_expr = replace_na(percent_expr, 0))
}


## TF expression --------------------------------------------------------

# TF list
TF_list <- readRDS(input_parameters$TF_list_path)

## Paralogs -------------------------------------------------------------

# get hsap paralog and macfas homolog info for these TFs from BioMart database
# gorGor6 not available in ensemble?
mart <- useEnsembl(biomart = "ensembl", 
                   dataset = "hsapiens_gene_ensembl")
orthoPara <- getBM(attributes = c("external_gene_name",
                                  "ensembl_gene_id",
                                  "hsapiens_paralog_associated_gene_name",
                                  "hsapiens_paralog_ensembl_gene",
                                  "hsapiens_paralog_perc_id",
                                  "hsapiens_paralog_subtype",
                                  "mmulatta_homolog_associated_gene_name",
                                  "mmulatta_homolog_ensembl_gene",
                                  "mmulatta_homolog_perc_id",
                                  "mmulatta_homolog_orthology_type",
                                  "ggorilla_homolog_associated_gene_name",
                                  "ggorilla_homolog_ensembl_gene",
                                  "ggorilla_homolog_perc_id",
                                  "ggorilla_homolog_orthology_type",
                                  "pabelii_homolog_associated_gene_name",
                                  "pabelii_homolog_ensembl_gene",
                                  "pabelii_homolog_perc_id",
                                  "pabelii_homolog_orthology_type"
                                  ),
                   filters    = "hgnc_symbol",
                   values     = TF_list, 
                   mart       = mart)

# keep only human-specific paralogs and keep only real paralogs that are not just alternative transcripts
hsap_gtf <- plyranges::read_gff("/data/share/htp/perturb-seq/genome_data/GRCh38_withdCas9/genes/genes.gtf")
table(orthoPara$hsapiens_paralog_subtype)

# 
paralogs <- orthoPara %>% 
  dplyr::filter(hsapiens_paralog_subtype %in% c("Primates", "Hominidae", "Homininae", "Hominoidea", "Homo sapiens")) %>% 
  rowwise() %>% 
  dplyr::mutate(exon_overlap_fraction = (plyranges::find_overlaps_within(hsap_gtf %>%  
                                                                           filter(type == "exon" & str_split(gene_id, "\\.", simplify = T)[,1] == ensembl_gene_id), hsap_gtf %>% 
                                                                           filter(type == "exon" & str_split(gene_id, "\\.", simplify = T)[,1] == hsapiens_paralog_ensembl_gene))  %>% 
                                           GenomicRanges::reduce(ignore.strand=T) %>% 
                                           width() %>% 
                                           sum()) / 
                  (min(hsap_gtf %>% 
                         filter(type == "exon" & str_split(gene_id, "\\.", simplify = T)[,1] %in% ensembl_gene_id) %>% 
                         GenomicRanges::reduce(ignore.strand=T) %>% 
                         width() %>% 
                         sum(),
                       hsap_gtf %>% 
                         filter(type == "exon" & str_split(gene_id, "\\.", simplify = T)[,1] %in% hsapiens_paralog_ensembl_gene) %>% 
                         GenomicRanges::reduce(ignore.strand=T) %>% 
                         width() %>% 
                         sum())))
paralogs <- paralogs %>% 
  dplyr::transmute(type = ifelse(hsapiens_paralog_subtype == "Primates", "primate specific", "human specific"),
                   paralog_gene_name = hsapiens_paralog_associated_gene_name,
                   paralog_ensembl_id = hsapiens_paralog_ensembl_gene,
                   percent_identity = format(hsapiens_paralog_perc_id, digits = 3, nsmall = 1),
                   percent_overlap = format(exon_overlap_fraction*100, digits = 3, nsmall = 1)) %>% 
  split(paralogs$external_gene_name)
saveRDS(paralogs, "TF_info/paralogs.rds")


## phastCons scores of longest CCDS -----------------------------------------------------

# gtf
hg38_gtf <- plyranges::read_gff("/data/share/htp/perturb-seq/genome_data/GRCh38_withdCas9/genes/genes.gtf")

# longest CCDS
TF_longest_CCDS <- hg38_gtf %>% 
  as_tibble() %>% 
  dplyr::filter(gene_name %in% TF_list & transcript_type == "protein_coding" & type == "CDS" & tag == "CCDS") %>% 
  group_by(gene_name, gene_id, transcript_id, ccdsid) %>% 
  dplyr::mutate(length_CDS = sum(width)) %>% 
  ungroup() %>% 
  distinct(gene_name, gene_id, ccdsid, length_CDS, seqnames, start, end, width, strand, source, transcript_support_level) %>% 
  group_by(gene_name) %>% 
  dplyr::filter(length_CDS == max(length_CDS)) %>% 
  dplyr::filter(length(unique(ccdsid)) == 1 | transcript_support_level == min(transcript_support_level, na.rm = T)) %>% 
  ungroup() %>% 
  distinct(gene_name, gene_id, ccdsid, seqnames, start, end, width, strand)
length(unique(TF_longest_CCDS$gene_name))
length(unique(TF_longest_CCDS$ccdsid))

# helper function (fixed)
analysePhastCons_v2 <- function( bigWigFile, gr, probcut=0.9){
  phyloP  <- rtracklayer::import(bigWigFile, 
                                 which= gr,
                                 as="NumericList")
  sumP<- sapply(phyloP, function(x){ 
    n<-length(x)
    c(n, sum(x), sum(x>probcut) )
  }) %>% t() %>% data.frame()
  names(sumP)<- c("bp","sumP","neg_n")
  
  sumP <- as_tibble(phyloP@metadata$ranges) %>% dplyr::select(-strand) %>% bind_cols(sumP) %>% distinct()
  
  gr %>% as_tibble() %>% inner_join(sumP, by = c("seqnames", "start", "end", "width"))
}

# calculate phastCons scores (with fixed function)
TF_longest_CCDS_phastCons <- analysePhastCons_v2(bigWigFile = "/data/share/htp/perturb-seq/TF_selection_gRNA_design_exp3/phastCons100way/hg38.phastCons100way.bw", gr = TF_longest_CCDS %>% as_granges(), 
                                              probcut = 0.9)


# summarise phastCons scores across all bps
TF_phastCons <- TF_longest_CCDS_phastCons %>%
  group_by(gene_name, gene_id)  %>% 
  dplyr::summarise(meanCons = sum(sumP) / sum(bp),
                   fracCons = sum(neg_n) / sum(bp),
                   bp = sum(bp)) %>% 
  ungroup()
saveRDS(TF_phastCons, "TF_info/phastCons.rds")


## GO terms -------------------------------------------------------------

go2gene <- topGO::annFUN.org("BP", mapping = "org.Hs.eg.db", ID = "symbol")
gene2go <- topGO::inverseList(go2gene)
tf2go <- gene2go[names(gene2go) %in% TF_list]
tf2go <- lapply(tf2go, 
                
                function(ids) {
                  
                  data.frame(GO_ID = ids,
                             GO_term = AnnotationDbi::Term(ids),
                             row.names = NULL) %>% 
                    dplyr::filter(!(GO_ID %in% c("GO:0000122", "GO:0006355", "GO:0006357", "GO:0010468", "GO:0045944", "GO:0045892", "GO:0045893", "GO:0010629")))
                  
                  
                })
saveRDS(tf2go, "TF_info/tf2go.rds")


## TF families

# the most common interpro terms grouped into TF families (loosely based on Vaquerizas et al., 2009, their assignment is in genomes/interpro2family_Vaquerizas.txt)
interpro2family <- readRDS("/data/share/htp/perturb-seq/TF_selection_gRNA_design_exp3/RDS/interpro2family.rds")

# interpro info from bioMart
mart <- useEnsembl(biomart = "ensembl", 
                   dataset = "hsapiens_gene_ensembl")
tf_families <- getBM(attributes=c("external_gene_name", "interpro_description"),
                     filters    = "hgnc_symbol",
                     values     = TF_list, 
                     mart       = mart) %>% 
  dplyr::rename(gene_name = external_gene_name) %>% 
  dplyr::filter(interpro_description != "")

# add TF family labels (label = other if the interpro term is rare)
tf_families <- tf_families %>% 
  left_join(interpro2family, by = c("interpro_description" = "interpro_term")) %>% 
  # drop_na() %>% 
  # dplyr::select(-interpro_description) %>% 
  # distinct() %>% 
  # right_join(data.frame(gene_name = tf.symb.list)) %>% 
  group_by(TF_family) %>% 
  dplyr::mutate(n_gene = length(unique(gene_name)),
                TF_family = ifelse(n_gene < 5 | is.na(TF_family), "other/not annotated", TF_family)) %>% 
  dplyr::select(-n_gene) %>% 
  group_by(gene_name) %>% 
  arrange(TF_family) %>% 
  ungroup() %>% 
  distinct()
saveRDS(tf_families %>% dplyr::select(-gene_name) %>% split(tf_families$gene_name), "TF_info/TF_families.rds")





