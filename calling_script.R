#!/usr/bin/env Rscript

# Import libraries
library(yaml)
library(optparse)

######################################
## Input parameters
######################################

# Setting list of options
option_list = list(make_option(c("--sourcecode_path"),
                               type = "character",
                               help = "Path to sourcecode funcions",
                               metavar = "character",
                               default = '/data/share/htp/perturb-seq/gRNA_design_workflow/sourcecode.R'),
                   make_option(c("--paramter_path"),
                               type = "character",
                               help = "Path to Input paramters YAML file",
                               metavar = "character",
                               default = '/data/share/htp/perturb-seq/gRNA_design_workflow/shiny_app/input_chimp.yaml')
                   )

# options parser
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

# Import sourcecode
source(opt$sourcecode_path)

# Read user paramters
input_parameters <- read_yaml(opt$paramter_path)

######################################
## Folder structure for Genome
######################################

# Check if directory already exists
dir_exists_flag <- dir.exists(sprintf('%s/outputs/%s', input_parameters$working_dir, input_parameters$genome))

if ((input_parameters$rewrite_flag == FALSE) & (dir_exists_flag))
{ print('Attempt to overwrite output with rewrite flag FALSE!!! \n Script will stop execution.')
  stop()
  
} else if((input_parameters$rewrite_flag == TRUE) & (dir_exists_flag)){
  print('Overwriting output since rewrite flag TRUE!!! \n Script will continue execution.')
}

# Call folder structure function
create_folder_structure(input_parameters)

# Change working directory to genome folder
input_parameters$working_dir <- paste0(input_parameters$working_dir, '/outputs/', input_parameters$genome)

# Set working directory to genome folder
setwd(input_parameters$working_dir)

######################################
## Load Inputs: GTF and Nanopore
######################################

# Load pre-selected Genes
TF_list <- readRDS(input_parameters$TF_list_path)

# GTF file
gtf_file <- plyranges::read_gff(input_parameters$gtf_file_path)

# Load the Nanopore data
if (input_parameters$np_expr_gff_path != ""){
  
  np <- plyranges::read_gff2(input_parameters$np_expr_gff_path)
  
  # Align seqnames of GTF and nanopore
  seqlevelsStyle(gtf_file) <- input_parameters$seq_level
  seqlevelsStyle(np) <- input_parameters$seq_level
  
} else{
  np <- data.table()
  
  # Set gtf_file style
  seqlevelsStyle(gtf_file) <- input_parameters$seq_level
  
}


######################################
## ATAC-seq Data
######################################

# Load ATAC-seq data
atac <- read_narrowpeaks(input_parameters$atac_seq_narrow_peaks_path)

# Set ATAC-seq data seqstyle
seqlevelsStyle(atac) <- input_parameters$seq_level

######################################
## Seqnames Check
######################################

# Perform check
seqnames_missing <- seqnames_check(gtf_file, atac, np)

# If missing seqnames found, Halt script
if (sum(sapply(seqnames_missing, length)) > 0){
  
  cat(sprintf("Missing seqnames found in %s: %s \n", names(seqnames_missing), (seqnames_missing)))
  
  print("Halting script execustion. Check and align input seqnames!!")
  
  stop()
}

######################################
## GTF TSS
######################################

# Extract Exons and TSS of gene list from loaded GTF file
#!!! check if column names need to be parameterized/standardized
exons_tss <- extract_tss(gtf_file, TF_list)
TF_gtf_exons = exons_tss$TF_gtf_exons
TF_gtf_tss = exons_tss$TF_gtf_tss
remove(exons_tss)

#!!! Script to create file structure for outputs -done

######################################
## Nanopore expression data
######################################

# Get Nanopore TSS
TF_np_tss <- generate_np_TSS(np, TF_gtf_exons, input_parameters)

######################################
## Joining TSS from all sources 
######################################

# Distance of np TSS with all sources
combined_sources <- combine_tss_sources(input_parameters, TF_gtf_tss, TF_np_tss = TF_np_tss)
dist_tss_np_all = combined_sources$dist_tss_np_all
TF_tss_all_sources = combined_sources$TF_tss_all_sources
remove(combined_sources)

# Merged and Filtered TSS
#!! deal with multiple TSSs - done; keeping
TF_TSS_filt <- create_TF_TSS(input_parameters, dist_tss_np_all, atac)

######################################
## gRNA Design Inputs
######################################

# Save inputs
gRNA_design_tool_inputs(input_parameters, TF_TSS_filt, gtf_file, np = np)


######################################
## Shell script for bowtie_indices
######################################

bowtie_indices_script(input_parameters)

######################################
## Shell script for gRNA Design
######################################

gRNA_design_script(input_parameters)
  
######################################
## Run Shell scripts
######################################
if (input_parameters$run_scripts == T){
  
  # Run Shell script - Bowtie indices
  system(sprintf('%s/shell_scripts/create_bowtie_indices.sh', input_parameters$working_dir))
  
  # Run Shell script - gRNA Design
  system(sprintf('%s/shell_scripts/design_grnas.sh', input_parameters$working_dir))
  
  # Process gRNA results
  #!!! Check if changing number of grnas to design changes the name
  grna_design_results(input_parameters)
  
}
#!!! Check if TF_list filter is done for only gtf file - done in gtf

######################################
## Create log file
######################################
sink(sprintf("%s/%s", input_parameters$working_dir, input_parameters$logfile_path))
sprintf("Log file: %s \n", input_parameters$genome)
sprintf("GTF Exons saved to: RDS/TF_exons_%s.rds", input_parameters$genome)
sprintf("GTF TSSs saved to: RDS/TF_tss_%s.rds", input_parameters$genome)
sprintf("Annotated nanopore saved to: RDS/TF_np_annot_%s_compare_evidence.rds", input_parameters$genome)
if (input_parameters$np_expr_gff_path != ''){
  nanopore_tss_checks(TF_np_tss, TF_list)}
sprintf("Nanopore TSSs saved to: RDS/TF_np_tss_annot_%s.rds", input_parameters$genome)
sprintf("TF_np_TSS distance plot saved to: figures/min_dist_np_other_tss_%s.png", input_parameters$genome)
sprintf("TRUE: Nanopore TSS that fall < %i of GTF TSS", input_parameters$max_dist_other_TSS)
print(table(dist_tss_np_all$dist_np_other < input_parameters$max_dist_other_TSS))
sprintf("Merged adn filtered TSS saved to: RDS/TF_tss_%s.rds", input_parameters$genome)
filtered_TSS_checks(input_parameters, TF_TSS_filt)
sprintf("Sources of evidence plot saved to: figures/sources_of_TSSs_%s.png", input_parameters$genome)
sprintf("gRNA tool inputs saved to: %s/design_input_files", input_parameters$working_dir)
sprintf("Shell script for bowtie indices saved to: %s/shell_scripts/create_bowtie_indices.sh", input_parameters$working_dir)
sprintf("Shell script for gRNA design saved to: %s/shell_scripts/design_grnas.sh", input_parameters$working_dir)
if (input_parameters$run_scripts == T){
  sprintf("gRNA design results saved to: %s/design_output_files/", input_parameters$working_dir)
  sprintf("gRNA design plots saved to: %s/figures/", input_parameters$working_dir)
}
sink()


######################################
## Run shiny app
######################################

# Launch shiny app
if (input_parameters$run_shiny_app) {
  
  #Create app inputs
  create_app_data_files(input_parameters)
  
  # Get Selected TSS list for each TF
  TF_list <- readRDS(sprintf("%s/RDS/TF_tss_%s.rds", input_parameters$working_dir, input_parameters$genome))%>% 
    as_tibble() %>% 
    pull(gene_name) %>% 
    sort(decreasing = FALSE) %>% 
    unique()
  
  # Source app
  source(input_parameters$shiny_helper_path)
  setwd(input_parameters$shiny_working_dir)
  source(sprintf("%s/app.R", input_parameters$shiny_working_dir))
  
  ### run the application ###
  shinyApp(ui = ui, server = server)
}
