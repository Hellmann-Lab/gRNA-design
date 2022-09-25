# Cross-species gRNA-design

## Introduction

This is a detailed documentation of the ***Cross-species gRNA design*** shiny app created for the Perturb-seq project of the Hellmmann-Enard Lab. As part of this project, we plan to perturb a selection of transcription factors (TFs) using single-cell CRISPRi screens in primate iPS cells, infer gene regulatory networks (GRNs) based on the outcome of the perturbations, then quantitatively compare these GRNs aross species. 

During the experimental design, we selected TFs based on previous data that might be interesting to study and designed single-guide RNAs (gRNAs) to target them in the human and cynomolgus macaque genomes. The app aims to present all information about the TF selection and gRNA design in a structured and interactive manner.

The main steps of the experimental design are summarised on *Figure 1*.

<p>
  <img align="center"
  src="gRNA_design_pipeline.svg"
  alt="pipeline"></p>
  <p align="center"><em><strong>Figure 1.</strong> Main steps of the TF selection and gRNA design</em></p>
  
We considered only those genes as potential targets that fulfill certain basic criteria:

 - located on an autosome
 - annotated as a transcription factor with at least one associated motif (based on the [JASPAR 2022 vertebrate core and unvalidated collections](#2) and the [IMAGE database](#3))
 - has TSSs with sufficient evidence in both the human and the cynomolgus macaque genomes
 
We found 1109 TFs that passed these criteria, these are the ones included in the app. We then hand-picked 76 TFs based on various characteristics (expression level, robustness and cross-species conservation of regulons in co-expression networks, protein sequence conservation, TF family annotations and functional importance in iPSCs) to experimentally perturb. 

We identified the TSSs of the 76 genes in the human and macaque genomes, then designed gRNAs using the model published in [*Horlbeck, 2016*](#1) to target these genomic loci. We compiled species-specific libraries by selecting gRNAs with the highest and most comparable predicted activity scores across species. The gRNA libraries contain 4 gRNAs for each of the TSS and each of the species, as well as non-targeting gRNAs as negative controls. These libraries will be transduced into human and cynomolgus macaque iPS cell lines that inducibly express dCAs9-KRAB to achieve the TF perturbations. Later we would like to expand the spectrum of the species to the gorilla and orang-utan as well. 
  
## Data

All data files and scripts required to run the app are available at: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7111658.svg)](https://doi.org/10.5281/zenodo.7111658)

## Input

### Expression and module robustness cutoffs

Minimum expression and module robustness required. The potential target genes in the drop-down menu are filtered based on these cutoffs.

- **Expression measure:** The % of cells expressing a given TF in the human/cynomolgus iPS cells based on unpublished scRNA-seq data.
- **Module robustness measure:** The robustness/preservation of the TF regulons based on the human and cynomolgus early neural differentiation lineage in our unpublished scRNA-seq data. We considered two aspects of robustness, originally described in [Langfelder et al, 2008](#4):  

  - Density: How well-connected are the genes of the regulon overally?
  - Connectivity: How consistent are the interaction patterns inferred within the module if we compare them across biological replicates from the same species?  
- The histograms show the expression and module robustness distributions for all TFs that passed our basic criteria and appear in the count matrix of our expression data.
- The darkgrey colouring on the histograms indicates TFs that pass the cutoff in the given species.
- To be included for the further analysis, TFs have to pass the cutoff in both species.

### Target gene

- Normally, the drop-down menu contains all TFs that passed the basic criteria, appear in the count matrix of our expression data and are above the expression and module robustness cutoffs.

- If both cutoffs are set to the minimum (expression cutoff = 0 and module robustness cuotff = -3), the drop-down menu contains all TFs that passed the basic criteria, regardless whether we have expression data about them or not.

- An interactive plot of cross-species network conservation is displayed to aid the choice of TFs. Mouse over data points to find out which TFs have the most conserved/diverged regulons.

### Genome

For each gene, gRNAs were designed in both the human (hg38) and the cynomolgus macaque (macFas6) genomes. Choose a genome to decide which design should be displayed.

### Include Horlbeck design?

If set to "Yes", the gRNAs from *Horlbeck et al., 2016* (lifted over from hg19 to hg38) are also displayed in the hg38 genome alongside our own design.

### Padding for genomic region

The number of bps to display on both sides of the promoter region on the Gviz plot.

## Output - TF characteristics
 
### Expression in iPSCs

The same histogram as before, with the expression level of the chosen TF marked in red.
  
### Module robustness

The same histogram as before, with the robustness measure of the chosen TF marked in red.
  
### Cross-species conservation:

Conservation of the TF regulons between human and cynomolgus based on the early neural differentiation lineage in our unpublished scRNA-seq data. We considered two aspects of conservation, originally described in [Langfelder et al, 2008](#4):

  - Density: Does a module that is well-connected in one species also stay well-connected in the other species?
  - Connectivity: How similar are the interaction pattern between the regulons of a TF in the two different species?
    
### TFBS motifs

TFBS motif logos based on the information content matrices for each annotated motif of the TF. The full-length information contents are displayed below the plot.
  
### Protein sequence conservation

Mean phastCons score of the TF CDS based on multiple alignments of 29 primate genome sequences to the human genome ([phastCons30way](#4)). The histogram shows the distribution of the meanCons scores of all genes in the drop-down menu, the score of the chosen TF is marked in red.
  
### Functional annotation

Annotated GO terms for the chosen TF. Terms describing TF activity or DNA binding and terms that are not end nodes in the GO hierarchy were excluded.
  
### Primate- and human-specific paralogs

Primate- and human-specific paralogs from the BioMart database. 
  
## Output - gRNA selection

### Gviz plot

A plot displaying the genomic region around the TSSs of the chosen gene. The width of the region can be adjusted via the input parameter *Padding for genomic region*.

Tracks:

- **GENCODE annotation**: Human (GRCh38.p13) or cynomolgus macaque (Macaca_fascicularis_6.0.105) reference annotation for the gene of interest. If the chosen gene has an annotated symbol, the transcripts are filtered by the symbol, otherwise the transcripts are filtered by the genomic region and strand.

- **Long-read RNA-seq data**: Nanopore reads and coverage from human and cynomolgus iPS cells.

- **Designed gRNAs**: The targeted genomic positions and predicted activity scores of the designed gRNAs.

- **ATAC-seq data**: ATAC-seq peaks from human and cynomolgus iPS cells (two individuals for each of the species).

### Scores & positions for the gRNAs

An interactive plot displaying the targeted genomic positions and predicted activity scores of the gRNAs designed for each of the TSS regions. 

- gRNAs selected in the final CRISPRi library are marked by green, the ones not selected are coloured red. 

- If the *Include Horlbeck design?* input parameter is set to "Yes", the gRNAs from [*Horlbeck, 2016*](#1) are also displayed (coloured black). Please note that the scores from the Horlbeck design and the scores from our own design are not directly comparable, because they are calculated based on slightly different models. For our purposes, the published model was re-trained using ATAC-seq peaks from human and cynomolgus iPS cells as accessibility data.

- gRNAs can be selected by clicking or brushing data points with the mouse. Detailed information about the selected gRNAs will appear in the table below.

### Selected gRNAs

A table summarising the available information about the gRNAs selected by the user.

Table entries:

- **sgID**: Unique identifier of the gRNA in the chosen genome (format: geneName_strand_genomicPosition.lengthWithPAM-tssId).

- **source**: The design where the gRNA was found.

  - horlbeck: the gRNAs in Horlbeck et al., 2016 lifted over from hg19 to hg38
  - hg38 design: our own design for hg38
  - mf6 design: our own design for macFas6
  
- **tss_id**: The transcriptional start site (TSS) region that the gRNA targets. We identified the TSSs by integrating evidence from the GENCODE annotation, long read RNA-seq and bulk ATAC-seq data as well as reciprocal best BLAT hits of the human TSS in case of the non-human primates. The TSS were ranked based on available evidence, with a lower number meaning more evidence, and named according to the ranks (format: geneName_rank). For the gRNA design, we merged TSSs of a gene that are closer than 1kb into TSS regions. We refer to these regions by the concatenated names of the merged TSSs (format: geneName_rank1,geneName_rank2).

- **predicted_activity**: Activity score calculated based on Horlbeck et al., 2016. 

  Model:

  - in case of source "horlbeck": the original elastic net in Horlbeck et al., 2016
  - in case of source "hg38 design" and "mf6 design": a modified version of this elastic net trained with only one type of accessibility data (human iPSC ATAC-seq data as a replacement of the original DNase data)
  
  The scores are within the range [0, 1] with a higher score meaning a higher activity/more efficient knock-down. The strongest predictor of the activity is the position relative to the TSS, including both distance from the TSS and avoidance of canonical nucleosome-occupied regions.

- **off_target_stringency**: The level of the most stringent off-target filter the gRNA passed, with a lower level meaning a more stringent filter/less off-target effects.

  Filters applied:

  - 31_nearTSS: passed if there is no off-target site within the TSS regions with a mismatch threshold of 31
  - 21_genome: passed if there is no off-target site within the whole genome with a mismatch threshold of 21
  - 31_2_nearTSS: passed if there is <=1 off-target site within the TSS regions with a mismatch threshold of 31
  - 31_3_nearTSS: passed if there are <=2 off-target sites within the TSS regions with a mismatch threshold of 31
  
  Levels:

  - 0: the gRNA passed the filters 31_nearTSS & 21_genome
  - 1: the gRNA passed the filter 31_nearTSS but not 21_genome
  - 2: the gRNA passed the filter 21_genome but not 31_nearTSS
  - 3: the gRNA passed the filter 31_2_nearTSS but not 31_nearTSS or 21_genome
  - 4: the gRNA passed the filter 31_3_nearTSS but not 31_nearTSS, 21_genome or 31_2_nearTSS
  
  If a gRNA does not pass any of these filters, it is discarded from the design output.

- **gRNA sequence**: The 20-nt long sequence of the gRNA. In case of source "horlbeck", this means the original hg19 sequence which occassionally differs from the sequence at the position where the gRNA is displayed in hg38 genome space.

- **match_in_other_genome**: Category based on the matching between the gRNAs designed for hg38 and the gRNAs designed for macFas6.

  - zero_mismatch: the gRNA has a perfectly matching counterpart in the other genome
  - one_mismatch: the gRNA has a counterpart with 1 mismatch in the other genome
  - no_cyno_gRNA/no_human_gRNA: the gRNA does not have a match with ≤1 mismatch in the other genome

- **is_in_lib**: "Yes" if the gRNA is in the final CRISPRi library, "No" if it is not. 

New gRNAs can be added to table by selecting data points on the *Scores & positions for the gRNAs* plot. The entries can be cleared using the button *Clear selected*. The table can also be searched, filtered and exported.

## References
<a id="1">[1]</a> 
Horlbeck, M. A., Gilbert, L. A., Villalta J. E. *et al*. (2016). **Compact and highly active next-generation libraries for CRISPR-mediated gene repression and activation** *eLife* **5**:e19760.

<a id="2">[2]</a> 
Castro-Mondragon J. A., Riudavets-Puig R., Rauluseviciute I. *et al*. (2022) **JASPAR 2022: the 9th release of the open-access database of transcription factor binding profiles** *Nucleic Acids Research* **50**:D165–D173.

<a id="3">[3]</a> 
Madsen J. G. S., Rauch A., Van Hauwaert E. L. *et al*. (2018) **Integrated analysis of motif activity and gene expression changes of transcription factors** *Genome Res.* **28**:243-255.

<a id="4">[4]</a> 
Langfelder P., Luo R., Oldham M. C., Horvath S. (2008) **Is my network module preserved and reproducible?** *PLoS Comput Biol* **7**:e1001057


<a id="5">[5]</a> 
Siepel, A., Bejerano, G., Pedersen, J. S. *et al.* (2005). **Evolutionarily conserved elements in vertebrate, insect, worm, and yeast genomes.** *Genome Research* **15**:1034–1050.

## `R` session info

``` r
R version 4.1.3 (2022-03-10)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Devuan GNU/Linux 3 (beowulf)

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/openblas/libblas.so.3
LAPACK: /usr/lib/x86_64-linux-gnu/libopenblasp-r0.3.5.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8       
 [4] LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
[10] LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] grid      stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] cowplot_1.1.1         universalmotif_1.12.4 plotly_4.10.0         Gviz_1.38.4          
 [5] DT_0.24               plyranges_1.14.0      GenomicRanges_1.46.1  GenomeInfoDb_1.30.1  
 [9] IRanges_2.28.0        S4Vectors_0.32.4      BiocGenerics_0.40.0   forcats_0.5.2        
[13] stringr_1.4.1         dplyr_1.0.10          purrr_0.3.4           readr_2.1.2          
[17] tidyr_1.2.0           tibble_3.1.8          ggplot2_3.3.6         tidyverse_1.3.2      
[21] shinyBS_0.61.1        shiny_1.7.2          

loaded via a namespace (and not attached):
  [1] readxl_1.4.1                backports_1.4.1             Hmisc_4.7-1                
  [4] systemfonts_1.0.4           BiocFileCache_2.2.1         lazyeval_0.2.2             
  [7] splines_4.1.3               crosstalk_1.2.0             BiocParallel_1.28.3        
 [10] digest_0.6.29               ensembldb_2.18.4            htmltools_0.5.3            
 [13] fansi_1.0.3                 magrittr_2.0.3              checkmate_2.1.0            
 [16] memoise_2.0.1               BSgenome_1.62.0             googlesheets4_1.0.1        
 [19] cluster_2.1.4               tzdb_0.3.0                  Biostrings_2.62.0          
 [22] modelr_0.1.9                matrixStats_0.62.0          prettyunits_1.1.1          
 [25] jpeg_0.1-9                  colorspace_2.0-3            blob_1.2.3                 
 [28] rvest_1.0.3                 rappdirs_0.3.3              textshaping_0.3.6          
 [31] haven_2.5.1                 xfun_0.32                   crayon_1.5.1               
 [34] RCurl_1.98-1.8              jsonlite_1.8.0              VariantAnnotation_1.40.0   
 [37] survival_3.4-0              glue_1.6.2                  gtable_0.3.1               
 [40] gargle_1.2.0                zlibbioc_1.40.0             XVector_0.34.0             
 [43] DelayedArray_0.20.0         scales_1.2.1                DBI_1.1.3                  
 [46] Rcpp_1.0.9                  viridisLite_0.4.1           xtable_1.8-4               
 [49] progress_1.2.2              htmlTable_2.4.1             foreign_0.8-82             
 [52] bit_4.0.4                   Formula_1.2-4               htmlwidgets_1.5.4          
 [55] httr_1.4.4                  RColorBrewer_1.1-3          ellipsis_0.3.2             
 [58] farver_2.1.1                pkgconfig_2.0.3             XML_3.99-0.10              
 [61] sass_0.4.2                  nnet_7.3-17                 dbplyr_2.2.1               
 [64] deldir_1.0-6                utf8_1.2.2                  labeling_0.4.2             
 [67] tidyselect_1.1.2            rlang_1.0.3                 later_1.3.0                
 [70] AnnotationDbi_1.56.2        munsell_0.5.0               cellranger_1.1.0           
 [73] tools_4.1.3                 cachem_1.0.6                cli_3.3.0                  
 [76] generics_0.1.3              RSQLite_2.2.15              broom_1.0.1                
 [79] fastmap_1.1.0               ragg_1.2.2                  yaml_2.3.5                 
 [82] knitr_1.40                  bit64_4.0.5                 fs_1.5.2                   
 [85] AnnotationFilter_1.18.0     KEGGREST_1.34.0             mime_0.12                  
 [88] xml2_1.3.3                  biomaRt_2.50.3              compiler_4.1.3             
 [91] rstudioapi_0.14             filelock_1.0.2              curl_4.3.2                 
 [94] png_0.1-7                   reprex_2.0.2                bslib_0.4.0                
 [97] stringi_1.7.8               GenomicFeatures_1.46.5      lattice_0.20-45            
[100] ProtGenerics_1.26.0         Matrix_1.4-1                vctrs_0.4.1                
[103] pillar_1.8.1                lifecycle_1.0.1             jquerylib_0.1.4            
[106] data.table_1.14.2           bitops_1.0-7                httpuv_1.6.5               
[109] rtracklayer_1.54.0          R6_2.5.1                    BiocIO_1.4.0               
[112] latticeExtra_0.6-30         promises_1.2.0.1            gridExtra_2.3              
[115] dichromat_2.0-0.1           MASS_7.3-58                 assertthat_0.2.1           
[118] SummarizedExperiment_1.24.0 rjson_0.2.21                withr_2.5.0                
[121] GenomicAlignments_1.30.0    Rsamtools_2.10.0            GenomeInfoDbData_1.2.7     
[124] parallel_4.1.3              hms_1.1.2                   rpart_4.1.16               
[127] MatrixGenerics_1.6.0        googledrive_2.0.0           biovizBase_1.42.0          
[130] Biobase_2.54.0              lubridate_1.8.0             base64enc_0.1-3            
[133] interp_1.1-3                restfulr_0.0.15    
```
