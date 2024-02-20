# gRNA Design Pipeline

**Table of Contents**

  - [Description](#description)
  - [Requirements](#requirements)
  - [Inputs](#inputs)
  - [Outputs](#outputs)
  - [Usage](#usage)
  - [R Session](#r_session)
  - [Contributions](#contributions)
 
## Description
The pipeline designs species-specific gRNAs for a given list of genes.
Designing gene-specific guide RNAs, RNA sequence that recognizes and directs the dCAS-KRAB
repressor complex, is a challenge in highthroughput experiments especially for non-human
species. Most CRISPRi studies chose pre-designed gRNAs published in Horlbeck et
al. (2016) (copyright CC BY-NC 4.0). Here, we generalize their CRISPRia design tool to be able to design gRNAs for other species.

The process has the following steps:

  1. Identifying Transcription Start sites: 
    - Based on the GTF and gene list provided. The TSSs are identified and used downstream. If nanopore expression data is provided, it is used as the primary evidence with the GTF providing annotation evidence (use TSSs with GTF annotation in vicinity).
    - Multiple TSSs for each gene are marked and kept.
  
  2. Designing guide RNA:
    - gRNA design is done using Horlbeck et al. (2016) (copyright CC BY-NC 4.0)
    
  3. Filtering guide RNA:
    - Selection of top gRNAs (4 per gene) based on predicted activity and off-target effects

## Requirements
The pipeline incorporates R and python scripts. The scripts utilize following versions:
R version 4.2.3 
Python 2.7.18

Use the following code to install requirements of the pipeline from requirements.txt

```bash
Rscript install_packages.R
```

## Inputs

The following files are required as inputs:

  1. Gene list: List of genes for which guide RNA design is to be run. The list is an important user input and has been mandated in the current version to avoid running the design script for the entire gene list in the genome.

  2. GTF: The `GTF` file is required to get a list of annotated TSSs along with the characteristics including gene position, length etc.
  
  3. ATAC-Seq: The ATAC-Seq data in `narrowpeaks` format is required to filter out active TSSs from the list of annotated TSS in `GTF`
  Note: ATAC-Seq data in bigwig file format is also required for plotting in the app.
  
  4. Genome: Genome `FASTA` file is required for bowtie indexing during gRNA design
  
  5. Expression data: Nanopore expression file in `GFF` format is optional if TSSs should be selected based on Nanopore expression data. The expression data is also added to the app if provided.
  
  6. YAML file: All inputs including filepaths and design parameters are provided through the YAML file.

<h2 id="#outputs">Outputs</h2>

### Shell scripts
The R script produces 2 shell scripts:
  - create_bowtie_indices.sh: Filtering genome FASTA to remove mitochondria genome and X chromosome, Bowtie indexing of genome, creation of 2bit file and bowtie indexing of tss regions
  - design_grnas.sh: Running gRNA design script based on Horlbeck, 2016.
  
### Log file
The log file contains the summary of script execution:
  - Paths of saved files and figures
  - Sanity checks on Nanopore data (if it exists)
  - Saniy checks on filtered and merged TSSs
  
### Designed gRNAs (if parameter `run_scripts == TRUE`)
The designed gRNAs are saved as:
  - unfiltered table
  - top n selected table for all genes
  - top 4 gRNAs for all genes based on predicted activity and off_target_stringency score
  
### Shiny app
The shiny app inputs are created and saved. The app is also launched if `run_shiny_app == TRUE`.
The app contains the multiple options and visualizations for targeted gRNA selection.

The user options are as follows:
  - Genome: Genome for which data is to be visualized
  - Gene: Particular gene to visualize
  - TSS: TSS to visualize in case of multiple TSSs. Select `all` to keep all
  - Gviz region: To view the entire gene (`entire gene body`) or only the TSS region (`TSS only`)
  - Padding for genome: Basepairs around thegene to include in the Gviz plot
  
Specifics in the app:
  1. Gviz plot: Visualize the genomic region of TSS. The Gviz plot has multiple tracks
    - GTF track: Gene annoataion
    - ATAC-Seq tracks: ATAC peaks for all samples
    - gRNA track: Target positions of deasigned gRNA
  2. gRNA plot: Scatterplot between gRNA target position and predicted activity
  3. Selection table: Dynamic selection table to accumulate selected gRNAs


## Usage

### Repository structure
The repository has following files and folders:
  - `calling_script.R`: R script for gRNA design
  - `calling_script_bash.sh`: Bash wrapper script for R script
  - `sourcecode.R`: Source functions for R script
  - functions: Folder containing Horlbeck gRNA design python scripts
  - shiny_app: Shiny app script; App inputs are created here
  - outputs: gRNA design script outputs are created here
  - SLURM_files: For use of bash script; SLURM files are saved here

### Filling the YAML file

All inputs to the script are provided through the YAML file. A template YAML file has been provided with the parameters and their explanations below:

```
# Flag to rewrite folder; if true rewrites output; The folder directory is preserved
rewrite_flag : TRUE

# Log file path
logfile_path : "logfile.txt"

# Working direcctory
working_dir: "/path/to/working/dir"

# TRUE: Run shell scripts for Horlbeck gRNA design 
# FALSE: Stop after shell script generation
run_scripts: TRUE

# GTF
gtf_file_path : "/path/to/gtf_file"

# TF list : required user input - Load pre-selected TFs
TF_list_path : "/path/to/TF_list"

# Nanopore GFF - load Gorilla iPSC Nanopore data
np_expr_gff_path : "/path/to/nanopore_gff"

# ATAC seq - load Gorilla iPSC ATAC-seq data
atac_seq_narrow_peaks_path : "/path/to/atac_seq_narrow_peaks"

# Genome - Organism genome
genome : 'gorGor6'

# Seqnames level of genome - Imporant to synchronize chromosome names before gRNA design
seq_level: 'UCSC'

# TSS thresholds
max_dist_other_TSS : 100 # Minimum distance between GTF and nanopore TSS required to accept nanopore TSS
stretch_for_TSS_merge : 100 # how close TSSs should be merged?

# gRNA Design
grna_num : 40 # Number of gRNAs to design
binwidth : 10 # Minimum distance between target sites of two gRNAs for same TSS
max_reduction_activity : 0.15 # Concession given to predicted activity to select gRNAs with lower off-target effects
gRNA_function_path : '/path/to/gRNA_function_python_file'
gRNA_model_path : '/path/to/gRNA_function_pickle_file'
genome_fa : '/path/to/genome_FASTA'
genome_2bit : '/path/to/genome_2bit' # If left empty, the code creates a 2bit genome file in the directory of genome FASTA file
atac_bw_path: '/path/to/atac_bigwig'

# shiny App
run_shiny_app : TRUE # If shiny app should be launched after gRNA design
shiny_working_dir : "/path/to/shiny_app_script"
shiny_helper_path : "/path/to/shiny_app_helper_functions"
```
The template is also saved in the repository as input.yaml

### Run via Bash wrapper script

The script can be run via the `calling_script_bash.sh`. The script has 4 variables:
  - sourcecode_path - File path of sourcecode.R file
  - paramter_path - File path of YAML parameter file
  - jobname - Name of the job name in SLURM. Example: 'grna_pipeline'
  - output_path - Output for SLURM file of submitted job
  
Follow these steps to run gRNA design:
  1. Fill and save YAML file with desired parameters
  2. Update the calling script bash script and save
  3. Check if the script has run permission
  4. Call the script in terminal and run.
  
The outputs are generated and saved within outputs under the folder named after the input data genome.

### Run via R calling script

The script can also be run directly via the `calling_script.R`. The inputs of the script are:
  - sourcecode_path - File path of sourcecode.R file
  - paramter_path - File path of YAML parameter file
  
They can be changed in the `option_list` directly. The run steps are as follows:
  1. Fill and save YAML file with desired parameters
  2. Run script.


## R Session Info

```
R version 4.2.3 (2023-03-15)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Devuan GNU/Linux 4 (chimaera)

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.13.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] grid      stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] cowplot_1.1.1          plotly_4.10.2          DT_0.28                shinyBS_0.61.1        
 [5] shiny_1.7.5            glue_1.6.2             docstring_1.0.0        GenomicFeatures_1.50.4
 [9] AnnotationDbi_1.60.2   Biobase_2.58.0         reticulate_1.31        gridExtra_2.3         
[13] ggplotify_0.1.1        lmerTest_3.1-3         lme4_1.1-34            Matrix_1.6-0          
[17] ggtext_0.1.2           data.table_1.14.8      rtracklayer_1.58.0     Gviz_1.42.1           
[21] patchwork_1.1.3        wesanderson_0.3.6      universalmotif_1.16.0  biomaRt_2.54.1        
[25] kableExtra_1.3.4       knitr_1.43             plyranges_1.18.0       yaml_2.3.7            
[29] Rsamtools_2.14.0       Biostrings_2.66.0      XVector_0.38.0         GenomicRanges_1.50.2  
[33] GenomeInfoDb_1.34.9    IRanges_2.32.0         ATACseqQC_1.22.0       S4Vectors_0.38.1      
[37] BiocGenerics_0.44.0    lubridate_1.9.2        forcats_1.0.0          stringr_1.5.0         
[41] dplyr_1.1.2            purrr_1.0.2            readr_2.1.4            tidyr_1.3.0           
[45] tibble_3.2.1           ggplot2_3.4.3          tidyverse_2.0.0        optparse_1.7.3        

loaded via a namespace (and not attached):
  [1] rappdirs_0.3.3                R.methodsS3_1.8.2             ragg_1.2.5                   
  [4] bit64_4.0.5                   DelayedArray_0.24.0           R.utils_2.12.2               
  [7] rpart_4.1.19                  KEGGREST_1.38.0               TFBSTools_1.36.0             
 [10] RCurl_1.98-1.12               AnnotationFilter_1.22.0       generics_0.1.3               
 [13] lambda.r_1.2.4                RSQLite_2.3.1                 bit_4.0.5                    
 [16] tzdb_0.4.0                    webshot_0.5.5                 xml2_1.3.5                   
 [19] httpuv_1.6.11                 SummarizedExperiment_1.28.0   DirichletMultinomial_1.40.0  
 [22] xfun_0.40                     jquerylib_0.1.4               hms_1.1.3                    
 [25] evaluate_0.21                 promises_1.2.1                fansi_1.0.4                  
 [28] restfulr_0.0.15               progress_1.2.2                caTools_1.18.2               
 [31] dbplyr_2.3.3                  DBI_1.1.3                     htmlwidgets_1.6.2            
 [34] futile.logger_1.4.3           ellipsis_0.3.2                crosstalk_1.2.0              
 [37] backports_1.4.1               GenomicScores_2.10.0          annotate_1.76.0              
 [40] deldir_1.0-9                  MatrixGenerics_1.10.0         vctrs_0.6.3                  
 [43] ensembldb_2.22.0              cachem_1.0.8                  withr_2.5.0                  
 [46] BSgenome_1.66.3               vroom_1.6.3                   checkmate_2.2.0              
 [49] GenomicAlignments_1.34.1      prettyunits_1.1.1             getopt_1.20.3                
 [52] svglite_2.1.1                 cluster_2.1.4                 lazyeval_0.2.2               
 [55] seqLogo_1.64.0                crayon_1.5.2                  labeling_0.4.2               
 [58] edgeR_3.42.4                  pkgconfig_2.0.3               nlme_3.1-162                 
 [61] ProtGenerics_1.30.0           nnet_7.3-19                   rlang_1.1.1                  
 [64] lifecycle_1.0.3               filelock_1.0.2                BiocFileCache_2.6.1          
 [67] AnnotationHub_3.6.0           dichromat_2.0-0.1             VennDiagram_1.7.3            
 [70] randomForest_4.7-1.1          matrixStats_1.0.0             graph_1.76.0                 
 [73] Rhdf5lib_1.20.0               boot_1.3-28.1                 base64enc_0.1-3              
 [76] png_0.1-8                     viridisLite_0.4.2             rjson_0.2.21                 
 [79] bitops_1.0-7                  R.oo_1.25.0                   KernSmooth_2.23-22           
 [82] rhdf5filters_1.10.1           blob_1.2.4                    regioneR_1.30.0              
 [85] gridGraphics_0.5-1            jpeg_0.1-10                   CNEr_1.34.0                  
 [88] scales_1.2.1                  memoise_2.0.1                 magrittr_2.0.3               
 [91] plyr_1.8.8                    zlibbioc_1.44.0               compiler_4.2.3               
 [94] BiocIO_1.8.0                  RColorBrewer_1.1-3            cli_3.6.1                    
 [97] ade4_1.7-22                   htmlTable_2.4.1               formatR_1.14                 
[100] Formula_1.2-5                 MASS_7.3-60                   tidyselect_1.2.0             
[103] stringi_1.7.12                textshaping_0.3.6             highr_0.10                   
[106] locfit_1.5-9.8                ChIPpeakAnno_3.32.0           latticeExtra_0.6-30          
[109] sass_0.4.7                    VariantAnnotation_1.44.1      polynom_1.4-1                
[112] tools_4.2.3                   timechange_0.2.0              parallel_4.2.3               
[115] rstudioapi_0.15.0             TFMPvalue_0.0.9               foreign_0.8-84               
[118] farver_2.1.1                  digest_0.6.33                 BiocManager_1.30.21.1        
[121] pracma_2.4.2                  Rcpp_1.0.11                   gridtext_0.1.5               
[124] BiocVersion_3.16.0            later_1.3.1                   httr_1.4.7                   
[127] motifStack_1.42.0             biovizBase_1.46.0             colorspace_2.1-0             
[130] rvest_1.0.3                   XML_3.99-0.14                 splines_4.2.3                
[133] yulab.utils_0.0.6             RBGL_1.74.0                   multtest_2.54.0              
[136] systemfonts_1.0.4             preseqR_4.0.0                 xtable_1.8-4                 
[139] jsonlite_1.8.7                nloptr_2.0.3                  futile.options_1.0.1         
[142] poweRlaw_0.70.6               UpSetR_1.4.0                  R6_2.5.1                     
[145] Hmisc_5.1-0                   pillar_1.9.0                  htmltools_0.5.6              
[148] mime_0.12                     fastmap_1.1.1                 minqa_1.2.5                  
[151] BiocParallel_1.34.2           interactiveDisplayBase_1.36.0 codetools_0.2-19             
[154] utf8_1.2.3                    bslib_0.5.1                   lattice_0.21-8               
[157] numDeriv_2016.8-1.1           curl_5.0.2                    gtools_3.9.4                 
[160] GO.db_3.16.0                  interp_1.1-4                  roxygen2_7.2.3               
[163] survival_3.5-5                limma_3.56.2                  rmarkdown_2.24               
[166] InteractionSet_1.26.1         munsell_0.5.0                 rhdf5_2.42.1                 
[169] GenomeInfoDbData_1.2.9        HDF5Array_1.26.0              reshape2_1.4.4               
[172] gtable_0.3.4   
```

## Contributions
The repository is written by Anita Térmeg and Daksh Pratap Singh Pamar