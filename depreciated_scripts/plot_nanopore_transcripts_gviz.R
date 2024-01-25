np <- readRDS("/data/share/htp/perturb-seq/TF_selection_gRNA_design_exp3/RDS/TF_np_annot_hg38.rds")

np_gene_models <- np %>% 
  as_tibble() %>% 
  transmute(chromosome = seqnames, 
            start, end, width, strand,
            transcript = transcript_id,
            symbol = gene_name)

np_track <- GeneRegionTrack(np_gene_models %>% filter(symbol == "POU5F1"), 
                chromosome = "chr6",
                genome = "hg38",
                # display gene name for each transcript
                showId = TRUE, 
                geneSymbol = TRUE, 
                # # if the symbol is not annotated, use the ENSEMBL ID as gene name
                # symbol = gn,
                name= "Nanopore",
                fontsize = 10, cex.title = 0.8,
                fill = "darkblue", col = "darkblue", col.axis = "black", col.title = "black")

plotTracks(np_track,
           collapseTranscripts = F, shape = "arrow", 
           from = 31164337 - 2000, 
           to = 31180731 + 2000,
           title.width = 1.1,
           col.grid='grey' ,
           fontsize=11,
           cex.main = 1)
