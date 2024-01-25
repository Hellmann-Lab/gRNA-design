# Import libraries
library(yaml)

# Set working directory
setwd('/data/share/htp/perturb-seq/gRNA_design_workflow/shiny_app/data_files/')
input_parameters <- read_yaml('/data/share/htp/perturb-seq/gRNA_design_workflow/shiny_app/input_chimp.yaml') 
workdir <- '/data/share/htp/perturb-seq/gRNA_design_workflow/shiny_app/data_files'

# Change working directory to genome folder
#input_parameters$working_dir <- paste0(input_parameters$working_dir, '/', input_parameters$genome)

# Source helper functions
source('/data/share/htp/perturb-seq/gRNA_design_workflow/shiny_app/helper_functions.R')

# Create data list
data4shiny <- list('organism' = list(genome = input_parameters$genome,
                                     gtf = list(Liftoff = input_parameters$gtf_file_path),
                                     nanopore_reads =  input_parameters$np_expr_gff_path,
                                     nanopore_transcripts = sprintf("%s/RDS/TF_np_annot_%s.rds", input_parameters$working_dir, input_parameters$genome),
                                     atac = list('organism' = input_parameters$atac_bw_path),
                                     tss = list("final" = sprintf("%s/RDS/TF_tss_%s.rds", input_parameters$working_dir, input_parameters$genome)),
                                     grnas = sprintf("%s/gRNAs/gRNAs.rds", input_parameters$working_dir), gene_models = c()))
                    
# Create species folder
sapply(names(data4shiny), function(i) {system(paste0("mkdir ", workdir, '/', data4shiny[[i]][["genome"]]))})

# Create sub folders
#sapply(names(data4shiny[['organism']]), function(i) {system(paste0("mkdir ", workdir, data4shiny[['organism']][["genome"]], '/' ,i))})

# 

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
      saveRDS(gene_models, paste0(workdir, data4shiny[[species]]$genome, "/gene_models/Nanopore.rds"))
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
set_up_data_for_shiny("organism", data4shiny, TF_list, workdir)
