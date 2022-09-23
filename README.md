# Cross-species gRNA-design

## Background

This repository contains all data and code required to run the ***Cross-species gRNA design app*** created for the Perturb-seq project of the Hellmmann-Enard Lab. As part of this project, we plan to perturb a selection of transcription factors (TFs) using single-cell CRISPRi screens in primate iPS cells, infer gene regulatory networks (GRNs) based on the outcome of the perturbations, then quantitatively compare these GRNs aross species. As a first step, we selected TFs based on previous data that might be interesting to perturb, identified their transcriptional start sites (TSSs) and designed single-guide RNAs (gRNAs) using the model published in [*Horlbeck, 2016*](#1) to target these genomic loci. We then compiled species-specific libraries by selecting gRNAs with the highest and most comparable predicted activity scores. The steps of the experimental design are summarised on *Figure 1*.

<p>
  <img align="center"
  src="gRNA_design_pipeline.pdf"
  alt="flowchart"></p>
  <p align="center"><em><strong>Figure 1.</strong> Main steps of the TF selection and gRNA design</em></p>

## User guide for the app

### Input parameters:

#### Step 1: choose expression and module robustness criteria

- **expression measure:** The % of cells expressing a given TF in the human/cyno iPS cells of the NPC differentiation data.

- **module robustness measure:** Overall connectedness (density) and consistency of inferred regulatory links between biological replicates of the same species (connectivity) based on the human and cyno early neural differentiation lineage in the NPC differentiation and EB data. Detailed steps:

  1. integrate the NPC differentiation data and the neural lineage of the EB data
  2. keep only iPS cells and the early stage of differentiation
  3. infer clonewise networks
  4. calculate average network and based on this, get the pruned regulons of the TFs
  5. calculate the Z-density for all human/cyno clones and the Z-connectivity for all human-human/cyno-cyno clone pairs
  6. calculate the mean across all clones or clone pairs (1 Z-density and 1 Z-connectivity measure per species per module)
  7. fit a linear regression line over the Z-connectivity - module size data points of all modules, separately in the two species and calculate residuals
  8. calculate the Z-scores of the Z-densities and the Z-scores of the Z-connectivity residuals, together for both species
  9. calculate the mean of the density and connectivity Z-scores in human/cyno for the module of interest
  
- The histograms show the expression and module robustness distributions of all TFs that...

  - ...have at least 1 annotated motif
  - ...have at least 1 well-supported TSS 
  - ...are not X- or Y-linked
  - ...are in the NPC differentiation data count matrix
  - ...have a regulon size >= 20
  
- The darkgrey colouring on the histograms indicates TFs that pass the cutoff in the given species.

- To be included for the further analysis, TFs have to pass the cutoff in both species.

#### Step 2: choose the the TF to be perturbed

- Normally, the drop-down menu contains all TFs that...

  - ...have at least 1 annotated motif
  - ...have at least 1 well-supported TSS 
  - ...are not X- or Y-linked
  - ...are in the NPC differentiation data count matrix
  - ...have a regulon size >= 20
  - ...pass the expression cutoff in both species
  - ...pass the module robustness cutoff in both species
  
- If both cutoffs are set to the minimum (expression cutoff = 0 and module robustness cuotff = -3), the drop-down menu contains all TFs that...

  - ...have at least 1 annotated motif
  - ...have at least 1 well-supported TSS 
  - ...are not X- or Y-linked
  
- An interactive plot of cross-species regulatory network conservation is displayed to aid the choice of TFs. Mouse over data points to find out which are the TFs with the most conserved/diverged regulons.

#### Step3: choose the genome

- For each gene, gRNAs were designed in both the human (hg38) and the cyno (macFas6) genome. Choose a genome to decide which design should be displayed.

#### Step4: choose additional parameters

- **Include Horlbeck design?**: If set to yes, the gRNAs from *Horlbeck et al., 2016* (lifted over from hg19 to hg38) are also displayed in the hg38 genome and appear as matched gRNAs in the macFas6 genome.

- **Padding for genomic region**: The number of bps to display on both sides of the TSS regions for the chosen TF.

- **Binning parameter for the merging of close-by gRNAs:** The number of bps the gRNAs will be extended/reduced to before finding overlaps between them and keeping only the highest scoring gRNA from each cohort. The higher the binning parameter, the less gRNAs will be kept and the larger the distance will be between them.

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



