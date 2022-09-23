# Cross-species gRNA-design

## Background

This repository contains all data and code required to run the ***Cross-species gRNA design app*** created for the Perturb-seq project of the Hellmmann-Enard Lab. As part of this project, we plan to perturb a selection of transcription factors (TFs) using single-cell CRISPRi screens in primate iPS cells, infer gene regulatory networks (GRNs) based on the outcome of the perturbations, then quantitatively compare these GRNs aross species. 

As a first step, we selected TFs based on previous data that might be interesting to perturb, identified their transcriptional start sites (TSSs) and designed single-guide RNAs (gRNAs) using the model published in [*Horlbeck, 2016*](#1) to target these genomic loci. We then compiled species-specific libraries by selecting gRNAs with the highest and most comparable predicted activity scores across species (*Figure 1*). These libraries will be transduced into human and cynomolgus macaque iPS cell lines that inducibly express dCAs9-KRAB to achieve TF repression. Later we would like to expand the spectrum of the species to the gorilla and orang-utan as well.

<p>
  <img align="center"
  src="gRNA_design_pipeline.svg"
  alt="pipeline"></p>
  <p align="center"><em><strong>Figure 1.</strong> Main steps of the TF selection and gRNA design</em></p>
  
Basic criteria for a target gene:

 - autosomal
 - annotated as a transcription factor with at least one associated motif
 - has TSSs with sufficient evidence in both the human and the cynomolgus macaque genomes (hg38 and macFas6, respctively)
 
We found 1109 TFs that passed these criteria, these are the ones included in the app. We then selected 76 TFs -- corresponding to 80 TSSs -- based on various characteristics (expression levels, robustness and cross-species conservation of regulons in co-expression networks, protein sequence conservation, TF family annotations and functional importance in iPSCs). The gRNA libraries contain 4 gRNAs for each of these 80 TSS and each of the 2 species, as well as non-targeting gRNAs as negative controls. The app aims to present all characteristics we collected for the 1109 candidate TFs, the designed gRNAs and our final selections in an easily searchable and visual way.
  
## Data



## User guide for the app

### Input parameters

#### Expression and module robustness cutoffs

- **expression measure:** The % of cells expressing a given TF in the human/cyno iPS cells based on unpublished scRNA-seq data.

- **module robustness measure:** The robustness/preservation of the TF regulons based on the human and cyno early neural differentiation lineage in our unpublished scRNA-seq data. We considered two aspects of robustness:

  1. Density: How well-connected are the genes of the regulon overally?
  2. Connectivity: How consistent are the interaction patterns inferred within the module if we compare them across biological replicates from the same species?
  
- The histograms show the expression and module robustness distributions for all TFs that passed our basic criteria and have available expression data.
  
- The darkgrey colouring on the histograms indicates TFs that pass the cutoff in the given species.

- To be included for the further analysis, TFs have to pass the cutoff in both species.

#### Target gene

- Normally, the drop-down menu contains all TFs that passed the basic criteria, have available expression data and are above the expression and module robustness cutoffs.
  
- If both cutoffs are set to the minimum (expression cutoff = 0 and module robustness cuotff = -3), the drop-down menu contains all TFs that passed the basic criteria, regardless of whether we have expression data about them or not.
  
- An interactive plot of cross-species regulatory network conservation is displayed to aid the choice of TFs. Mouse over data points to find out which are the TFs with the most conserved/diverged regulons.

#### Genome

For each gene, gRNAs were designed in both the human (hg38) and the cyno (macFas6) genome. Choose a genome to decide which design should be displayed.

#### Include Horlbeck design?

If set to yes, the gRNAs from *Horlbeck et al., 2016* (lifted over from hg19 to hg38) are also displayed in the hg38 genome alongside our own design.

#### Padding for genomic region

The number of bps to display on both sides of the promoter region for the chosen TF on the Gviz plot.

### Output

#### TF characteristics
 
  - **Expression in iPSCs:** The same histogram as before, with the expression of the chosen TF marked in red.
  
  - **Module robustness:** The same histogram as before, with the robustness measure of the chosen TF marked in red.
  
  - **Cross-species conservation:** Preservation of overall connectedness (density) and connectivity patterns (connectivity) between human and cyno based on the early neural differentiation lineage in the NPC differentiation and EB data. Detailed steps:

    1. integrate the NPC differentiation data and the neural lineage of the EB data
    2. keep only iPS cells and the early stage of differentiation
    3. infer clonewise networks
    4. calculate average network and based on this, get the pruned regulons of the TFs
    5. calculate the Z-density for all human/cyno clones and the Z-connectivity for all human-human/cyno-cyno clone pairs
    6. test density-based conservation: perform two-sided t-tests to determine if there is a significant difference between the Z-denisties of human clones and the Z-densities of cyno clones
    7. test connectivity-based conservation: perform one-sided t-tests to determine if  Z-connectivities of within-species clone pairs are significantly higher than the Z-connectivities of across-species clone pairs
    8. in both cases, the p-value of the t-test is a measure of conservation, with a lower value meaning stronger divergence
    
  - **TFBS motifs:** TFBS motif logos based on the information content matrices for each annotated motif of the TF. The full-length information contents are displayed below the plot.
  
  - **Protein sequence conservation:** Mean phastCons score of the TF CDS based on multiple alignments of 29 primate genome sequences to the human genome (phastCons30way). The histogram shows the distribution of the meanCons scores of all genes in the drop-down menu, the score of the chosen TF is marked in red.
  
  - **Functional annotation:** Annotated GO terms for the chosen TF. Terms describing TF activity or DNA binding and terms that are not end nodes in the GO hierarchy were excluded.
  
  - **Paralogs:** Primate- and human-specific paralog from the BioMart database. 
  
#### gRNA selection

## References
<a id="1">[1]</a> 
Horlbeck, M. A., Gilbert, L. A., Villalta J. E. *et al*. (2016). **Compact and highly active next-generation libraries for CRISPR-mediated gene repression and activation** *eLife* **5**:e19760.



