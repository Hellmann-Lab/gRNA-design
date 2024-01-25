#!/bin/bash

# Parameters and paths
sourcecode_path=/data/share/htp/perturb-seq/gRNA_design_workflow/sourcecode.R
paramter_path=/data/share/htp/perturb-seq/gRNA_design_workflow/shiny_app/input_chimp.yaml
jobname=chimp_grna_pipeline
output_path=/data/share/htp/perturb-seq/gRNA_design_workflow/SLURM_files

# Submit job
sbatch --mem=30G --cpus-per-task=10 --job-name=${jobname} --error=${output_path}/chimp_grna_pipeline.%J.err --output=${output_path}/chimp_grna_pipeline.%J.out --wrap="Rscript /data/share/htp/perturb-seq/gRNA_design_workflow/calling_script.R --sourcecode_path ${sourcecode_path} --paramter_path ${paramter_path}"