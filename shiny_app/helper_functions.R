# Import libraries
library(docstring)
library(tidyverse)
library(plyranges)
library(glue)


#### TF CHARACTERISTICS TAB ####


### Plot expression histogram with the chosen expression cutoff ###
pull_TFs_by_gen <- function(expr, gen, percent_expr){
  #'Returns a filtered TF list 
  #'
  #'Uses expression threshold and selected genome to filter down TF list.  
  #'
  #'@param expr Expression table for all genomes
  #'@param gen Genome
  #'@param percent_expr Expression threshold
  
  species <- case_when(gen == 'hg38' ~ 'human',
                       gen == 'macFas6' ~ 'cyno',
                       gen == 'gorGor6' ~ 'gorilla',
                       gen == 'ponAbe3' ~ 'orang')
  
  # Create column for expression column
  expr_col <- paste0('percent_expr_', species)
  
  return(sort(expr %>% 
                filter(!! sym(expr_col) >= percent_expr) %>% 
                dplyr::pull(gene_name))) 
  
}

plot_expr_hist <- function(expr, percent_expr, is_expr, species_name) {
  #'Returns a histogram as ggplot object
  #'
  #'Creates a histogram for % of cells with number of TFs expressed 
  #'
  #'@param expr Expression table for all genomes
  #'@param percent_expr Expression threshold
  #'@param is_expr Name of logical column to include gene in plot; genome specific
  #'@param species_name Name of species 
  
  exhist <- expr %>%
    dplyr::mutate(species = species_name) %>%
    ggplot(aes(x = !! sym(percent_expr), fill = !! sym(is_expr))) +
    geom_histogram(col = "white",
                   binwidth = 2.5,
                   center = 1.25) +
    scale_fill_manual(values = c("grey30", "grey80"), guide = "none")  +
    theme_bw(base_size = 14) +
    ylab("Number of TFs expressed") +
    facet_grid( ~ species) +
    theme(axis.title.x = element_blank())
  return(exhist)
}

plot_expr_init <- function(expr, cutoff) {
  #'Returns a grid of histograms
  #'
  #'Creates a histogram for % of cells with number of TFs expressed for all species 
  #'
  #'@param expr Expression table for all genomes
  #'@param cutoff Expression threshold
  
  # if the expression cutoff is set to minimum, color the entire plot dark grey - percent_expr_human > cutoff (see below) would not cover everything
  if (cutoff == 0) {
    
    expr <- expr %>%
      dplyr::mutate(is_expr_human = T,
                    is_expr_cyno = T,
                    is_expr_gorilla = T,
                    is_expr_orang = T)
    
    # otherwise colour bins above the cutoff dark grey, bins below cutoff light grey
  } else{
    
    expr <- expr %>% 
      dplyr::mutate(is_expr_human = factor(percent_expr_human > cutoff, levels = c(T, F)), # > instead of >= for aesthetic reasons, otherwise bins might get 2 colours
                    is_expr_cyno = factor(percent_expr_cyno > cutoff, levels = c(T, F)),
                    is_expr_gorilla = factor(percent_expr_gorilla > cutoff, levels = c(T, F)),
                    is_expr_orang = factor(percent_expr_orang > cutoff, levels = c(T, F)))
    
  }
  
  # histogram for human
  exhist_h <- plot_expr_hist(expr, percent_expr = 'percent_expr_human', 
                             is_expr = 'is_expr_human', species_name = "human")
  
  # histogram for cynomolgus
  exhist_c <- plot_expr_hist(expr, percent_expr = 'percent_expr_cyno', 
                             is_expr = 'is_expr_cyno', species_name = "cynomolgus")
  
  # histogram for gorilla
  exhist_g <- plot_expr_hist(expr, percent_expr = 'percent_expr_gorilla', 
                             is_expr = 'is_expr_gorilla', species_name = "gorilla")
  
  # histogram for orang
  exhist_o <- plot_expr_hist(expr, percent_expr = 'percent_expr_orang', 
                             is_expr = 'is_expr_orang', species_name = "orang")
  
  title <- ggdraw() +
    draw_label("Percent of cells",
               size = 15,
               x = 0.5,
               hjust = 0) +
    theme(plot.margin = margin(0, 40, 10, 0),
          plot.background = element_rect(fill = "white", color = "transparent"))
  
  plot_grid(plot_grid(exhist_h,
                      exhist_c,
                      exhist_g,
                      exhist_o,
                      rel_widths = c(1.05, 1, 1, 1)),
            title,
            ncol = 1,
            rel_heights = c(1, 0.05))
  
}


### Plot expression histogram with the chosen expression cutoff and the gene of interest marked ###

plot_expr <- function(gn, expr, cutoff) {
  #'Returns a grid of histograms
  #'
  #'Creates a histogram for % of cells with number of TFs expressed for all species 
  #'
  #'@param gn Gene name
  #'@param expr Expression table for all genomes
  #'@param cutoff Expression threshold
  
  # if the expression cutoff is set to minimum, color the entire plot dark grey - percent_expr_human > cutoff (see below) would not cover everything
  if (cutoff == 0) {
    
    expr <- expr %>%
      dplyr::mutate(is_expr_human = T,
                    is_expr_cyno = T,
                    is_expr_gorilla = T,
                    is_expr_orang = T)
    
    # otherwise colour bins above the cutoff dark grey, bins below cutoff light grey
  } else{
    
    expr <- expr %>%
      dplyr::mutate(is_expr_human = factor(percent_expr_human > cutoff, levels = c(T, F)), # > instead of >= for aesthetic reasons
                    is_expr_cyno = factor(percent_expr_cyno > cutoff, levels = c(T, F)),
                    is_expr_gorilla = factor(percent_expr_gorilla > cutoff, levels = c(T, F)),
                    is_expr_orang = factor(percent_expr_orang > cutoff, levels = c(T, F)))
    
  }
  
  # histogram for human
  exhist_h <- plot_expr_hist(expr, percent_expr = 'percent_expr_human', 
                             is_expr = 'is_expr_human', species_name = "human")
  
  # histogram for cynomolgus
  exhist_c <- plot_expr_hist(expr, percent_expr = 'percent_expr_cyno', 
                             is_expr = 'is_expr_cyno', species_name = "cynomolgus")
  
  # histogram for gorilla
  exhist_g <- plot_expr_hist(expr, percent_expr = 'percent_expr_gorilla', 
                             is_expr = 'is_expr_gorilla', species_name = "gorilla")
  
  # histogram for orang
  exhist_o <- plot_expr_hist(expr, percent_expr = 'percent_expr_orang', 
                             is_expr = 'is_expr_orang', species_name = "orang")
  
  title <- ggdraw() +
    draw_label("Percent of cells",
               size = 15,
               x = 0.5,
               hjust = 0) +
    theme(plot.margin = margin(0, 40, 10, 0),
          plot.background = element_rect(fill = "white", color = "transparent"))
  
  plot_grid(plot_grid(exhist_h +
                        # add dashed line and gene label
                        geom_segment(x = expr$percent_expr_human[expr$gene_name == gn], xend = expr$percent_expr_human[expr$gene_name == gn], y = max(ggplot_build(exhist_h)$data[[1]]$count)*1.02, yend = 0, linetype = "dashed", size = 0.4, color = "red") +
                        geom_label(data = expr %>% dplyr::filter(gene_name == gn), aes(label = gene_name, y = max(ggplot_build(exhist_h)$data[[1]]$count)*1.05), size = 2.5, color = "red", fill = "transparent"),
                      exhist_c +
                        # add dashed line and gene label
                        geom_segment(x = expr$percent_expr_cyno[expr$gene_name == gn], xend = expr$percent_expr_cyno[expr$gene_name == gn], y = max(ggplot_build(exhist_c)$data[[1]]$count)*1.02, yend = 0, linetype = "dashed", size = 0.4, color = "red") +
                        geom_label(data = expr %>% dplyr::filter(gene_name == gn), aes(label = gene_name, y = max(ggplot_build(exhist_c)$data[[1]]$count)*1.05), size = 2.5, color = "red", fill = "transparent"),
                      exhist_g +
                        # add dashed line and gene label
                        geom_segment(x = expr$percent_expr_gorilla[expr$gene_name == gn], xend = expr$percent_expr_gorilla[expr$gene_name == gn], y = max(ggplot_build(exhist_g)$data[[1]]$count)*1.02, yend = 0, linetype = "dashed", size = 0.4, color = "red") +
                        geom_label(data = expr %>% dplyr::filter(gene_name == gn), aes(label = gene_name, y = max(ggplot_build(exhist_g)$data[[1]]$count)*1.05), size = 2.5, color = "red", fill = "transparent"),
                      
                      exhist_o +
                        # add dashed line and gene label
                        geom_segment(x = expr$percent_expr_orang[expr$gene_name == gn], xend = expr$percent_expr_orang[expr$gene_name == gn], y = max(ggplot_build(exhist_o)$data[[1]]$count)*1.02, yend = 0, linetype = "dashed", size = 0.4, color = "red") +
                        geom_label(data = expr %>% dplyr::filter(gene_name == gn), aes(label = gene_name, y = max(ggplot_build(exhist_o)$data[[1]]$count)*1.05), size = 2.5, color = "red", fill = "transparent"),
                      
                      rel_widths = c(1.05, 1)),
            title,
            ncol = 1,
            rel_heights = c(1, 0.05))
  
}


### Get full-length IC scores for all motifs of the chosen gene ###

get_IC <- function(gn, motif_IC) {
  #'Returns a plot with TFBS motif
  #'
  #'Creates a Motif plot for given gene 
  #'
  #'@param gn Gene name
  #'@param motif_IC Motif information table
  
  
  motif_IC %>% 
    dplyr::filter(SYMBOL == gn) %>% 
    dplyr::mutate(text = paste0(motif_id, ": IC = ", format(IC, digits = 3))) %>% 
    pull(text) %>% 
    paste(collapse = "\n")
  
}


### Plot phastCons histogram with the chosen gene marked ###

plot_phastCons <- function(gn, gn_list, phastCons) {
  #'Returns a phastCons histogram plot
  #'
  #'Creates a phastCons histogram with the chosen gene marked
  #'
  #'@param gn Gene name
  #'@param gn_list Filtered TF list
  #'@param phastCons phastCons information table
  
  p <- phastCons %>% 
    dplyr::filter(gene_name %in% gn_list) %>% 
    ggplot(aes(x=meanCons)) +
    geom_histogram(col="white", binwidth = 0.025, center = 0.0125, fill = "grey30") +
    theme_bw(base_size = 14) +
    ylab("Number of TFs")
  
  # add dashed line and gene label
  p + geom_segment(x = phastCons$meanCons[phastCons$gene_name == gn], 
                   xend = phastCons$meanCons[phastCons$gene_name == gn], 
                   y = max(ggplot_build(p)$data[[1]]$count)*1.02, 
                   yend = 0, 
                   linetype = "dashed", size = 0.4, color = "red") +
    geom_label(data = phastCons %>% dplyr::filter(gene_name == gn), aes(label = gene_name, y = max(ggplot_build(p)$data[[1]]$count)*1.05), size = 2.5, color = "red", fill = "transparent")
  
}

#### gRNA-DESIGN TAB ####

### Gviz plot ###

set_seqlevels <- function(genome_ranges, seqnames_style){
  #'Returns a genome range with suitable seqlevelsStyle
  #'
  #'Changes the seqlevelsStyle of the genome ranges to match given input. Can be 
  #'improved to automatically adjust to the most common style. 
  #'
  #'@param genome_ranges List of genome ranges available for the species
  #'@param seqnames_style Name of seqlevelStyle to be used for Gviz
  
  if (seqnames_style == 'NCBI'){
    return(genome_ranges %>% 
             dplyr::mutate(chromosome = gsub("chr", "", chromosome)))
  }
  
  else if (seqnames_style == 'UCSC'){
    genome_ranges <- genome_ranges %>% 
      dplyr::mutate(chromosome = gsub("chr", "", chromosome)) %>% dplyr::mutate(chromosome = paste("chr", chromosome, sep =""))
    return(genome_ranges)
  }
}

make_gRNA_gviz <- function(gn, tss_sel, gen, offset, gene_models, tss, grnas, atacBW, nanoporeBAM, gviz_focus){
  #'Returns a Gviz plot
  #'
  #'Creates a Gviz plot for given TF TSS. Includes ATAC-Seq and available gene annotations
  #'
  #'@param gn TF/gene name
  #'@param gen Genome
  #'@param offset Padding value for Gviz around TSS
  #'@param gene_models List of genome models available for the species
  #'@param tss TSS as genome range
  #'@param grnas Table of grnas designed for the particular TF-TSS
  #'@param atacBW ATAC-Seq peak data
  #'@param nanoporeBAM Nanopore expression data
  #'@param gviz_focus Logical trigger to plot Gviz from app
  
  # Ensure the chromosome names of all tracks have the same format
  if (gen %in% c('hg38', 'macFas6')){
    seqnames_style = 'NCBI'
  }
  
  else {
    seqnames_style = 'UCSC'
  }
  
  # Gene models
  for (num_gm in 1:length(gene_models)){
    
    gene_models[[num_gm]] <- set_seqlevels(gene_models[[num_gm]], seqnames_style)
  }
  
  # split data by is_in_final_lib
  grnas_filt <- grnas %>% 
    plyranges::filter(gene == gn) %>% 
    dplyr::mutate(is_in_final_lib = factor(is_in_final_lib, levels = c(F, T)))
  
  # Filter only if one is selected
  if (tss_sel != 'all'){
    grnas_filt <- grnas_filt %>% 
      plyranges::filter(transcript == tss_sel)
  }
  seqlevelsStyle(grnas_filt) <- seqnames_style
  gRNA_data_split <- split(grnas_filt, as_tibble(grnas_filt)$is_in_final_lib) 
  
  # pull out the TSS regions for the selected gene
  tss_filt <- tss[["final"]] %>% 
    plyranges::filter(gene_name == gn) %>% 
    arrange(rnk)
  
  # Filter only if one is selected
  if (tss_sel != 'all'){
    tss_filt <- tss_filt %>% 
      plyranges::filter(tss_id == tss_sel)
  }
  
  
  # filter gene models for the gene of interest
  gene_models_filt <- lapply(gene_models, function(gene_model) {
    
    gene_model %>% 
      dplyr::filter(symbol == gn)
  })
  
  names(gene_models_filt) <- names(gene_models)
  
  # get chromosome
  if (nrow(as_tibble(gene_models_filt[["Nanopore"]])) > 0) {
    
    chr <- as.character(unique(as_tibble(gene_models_filt[["Nanopore"]])$chromosome))
    
  } else {
    
    chr <- as.character(unique(as_tibble(gene_models_filt[["Liftoff"]])$chromosome))
    
  }
  
  # get region boundaries
  
  if (gviz_focus == "tss" & nrow(as_tibble(tss_filt)) > 0) {
    
    min_pos <- min(as_tibble(tss_filt)$start)
    
  } else {
    
    min_pos <- lapply(c(list(tss_filt), gene_models_filt), function(gr) {
      
      as_tibble(gr)$start
      
    }) %>% unlist() %>% min(na.rm = T)
    
  }
  
  min <- min_pos - offset
  
  if (gviz_focus == "tss" & nrow(as_tibble(tss_filt)) > 0) {
    
    max_pos <- max(as_tibble(tss_filt)$end)
    
  } else {
    
    max_pos <- lapply(c(list(tss_filt), gene_models_filt), function(gr) {
      
      as_tibble(gr)$end
      
    }) %>% unlist() %>% max(na.rm = T)
    
  }
  
  max <- max_pos + offset
  
  # set correct chromosome format
  options(ucscChromosomeNames = F)
  
  # plot axis with genomic coordinates
  axis  <- GenomeAxisTrack(genome = gen)
  
  # plot GENCODE annotation
  gene_model_tracks <- lapply(names(gene_models), function(name) {
    
    # gene_model_attrs <-  gene_models[[name]] %>%
    #                         filter(chromosome == chr) %>% 
    #                         select(symbol)
    # gene_model_attrs <- gene_model_attrs %>% 
    #                         mutate(fill.col = ifelse(symbol == gn, 
    #                                                  'darkblue', 'darkgray'))
    
    GeneRegionTrack(gene_models[[name]], 
                    chromosome = chr,
                    genome = gen,
                    # display gene name for each transcript
                    showId = TRUE, 
                    geneSymbol = TRUE, 
                    # # if the symbol is not annotated, use the ENSEMBL ID as gene name
                    # symbol = gn,
                    name= paste0(name, "\ntranscripts"),
                    fontsize = 10, cex.title = 0.8,
                    fill = 'darkgray', 
                    col = 'darkgray', 
                    col.axis = "black", 
                    col.title = "black")
    
    
  })
  names(gene_model_tracks) <- names(gene_models)

  # Add features to tracks
  for (name in names(gene_model_tracks)){
    feature(gene_model_tracks[[name]]) <- ifelse(symbol(gene_model_tracks[[name]]) == gn, "highlight", "rest")
  }
  
  
  # gene_tracks <- lapply(names(gene_models_filt), function(name) {
  #   
  # highlight_gene <- gene_models_filt[[name]] %>% 
  #   filter(symbol == gn)
  # 
  # highlighted_gene_track <- GeneRegionTrack(highlight_gene, 
  #                                           chromosome = chr,                           
  #                                           genome = gen,
  #                                           # display gene name for each transcript
  #                                           showId = TRUE, 
  #                                           geneSymbol = TRUE, 
  #                                           # # if the symbol is not annotated, use the ENSEMBL ID as gene name
  #                                           # symbol = gn,
  #                                           name= paste0(name, "\ntranscripts"),
  #                                           fontsize = 10, cex.title = 0.8,
  #                                           fill = "darkblue", 
  #                                           col = "darkblue", 
  #                                           col.axis = "black", 
  #                                           col.title = "black")
  # 
  # })
  # names(gene_tracks) <- names(gene_models_filt)
  
  # plot ATAC-seq coverage 
  atac_colors = c("organism"="tan2")
  names(atac_colors) <- input_parameters$genome
  
  atac_tracks <- lapply(names(atacBW), function(name) {
    
    spec <- strsplit(name,"\\_")[[1]][1]
    
    DataTrack( range = atacBW[[name]], 
               type = 'h', 
               chromosome = chr,
               ## select colour according to the species
               col = unname(atac_colors[spec]) ,
               name = paste0("ATAC-seq\n(", name, ")"),
               window = -1, windowSize = 100, genome = gen,
               col.title="black", cex.title=0.8,
               col.axis="black") 
    
  })
    
  
  # If nanopore data does not exist then pass
  if (!is.na(nanoporeBAM)){
    # plot nanopore reads and coverage
    nanopore_track <- AlignmentsTrack(nanoporeBAM,
                                      chromosome = chr,
                                      genome = gen,
                                      name = "Nanopore\nreads", 
                                      height=0.2, coverageHeight = 0.1, minCoverageHeight = 0,
                                      window = -1, windowSize = 100,
                                      cex.title=0.8,
                                      col.axis="black", col.title="black")
    
    list_tracks <- c(gene_model_tracks[names(gene_model_tracks)], nanopore_track, atac_tracks)
  } else {
    list_tracks <- c(gene_model_tracks[names(gene_model_tracks)], atac_tracks)
    nanopore_track <- NULL
  }
  
  
  # highlight TSS on LIFTOFF/GENCODE/ENSEMBL + nanopore + ATAC tracks
  if (nrow(as_tibble(tss_filt)) > 0) {
    
    tss_track <- HighlightTrack(trackList = list_tracks,
                                start = start(tss_filt) - 10, 
                                end   = end(tss_filt) + 10,
                                chromosome = chr,
                                fill = c("chartreuse4", "chartreuse2", "darkolivegreen1", "darkseagreen1")[1:nrow(as_tibble(tss_filt))],
                                col  = c("chartreuse4", "chartreuse2", "darkolivegreen1", "darkseagreen1")[1:nrow(as_tibble(tss_filt))])
    
    
  } else {
    
    tss_track <- list_tracks
  }
  
  if (nrow(as_tibble(grnas_filt)) > 0) {
    
    # plot the 3 groups (selected in final library, not selected in final library, horlbeck) separately with different colours
    gRNA_scores <- lapply(names(gRNA_data_split), function(n) { 
      
      DataTrack(gRNA_data_split[[n]] %>% 
                  plyranges::select(predicted_score) %>% 
                  plyranges::mutate(strand = "+"),
                chromosome = chr, 
                genome = gen,
                col = c("grey70", "red2")[which(names(gRNA_data_split) == n)],
                name="Designed \ngRNAs",
                col.axis="black", col.title="black", cex.title=0.8 , 
                ylim =c(min(as_tibble(grnas_filt)$predicted_score)*0.9, max(as_tibble(grnas_filt)$predicted_score)*1.1))  
      
    } )
    
    ## overlay the 3 tracks
    grna_track <- OverlayTrack( gRNA_scores ,size=4, alpha=0.5, name="gRNA\n scores",
                                window = -1, windowSize = 100,
                                legend=T )
    
  } else {
    
    grna_track <- NULL
    
  }
  
  # combine into a single plot
  # Create track list
  track.list <- c(axis, tss_track, grna_track)
  
  # get sizes for each track
  if (!is.null(nanopore_track)){
    track.sizes <- c(0.7, rep(1.2, length(gene_model_tracks)), 3, rep(1.2, length(atac_tracks) + !is.null(grna_track)))
  } else {
    track.sizes <- c(0.7, rep(1.2, length(gene_model_tracks)), rep(1.2, length(atac_tracks) + !is.null(grna_track)))
  }
 
  plotTracks(track.list, 
              collapseTranscripts = F, shape = "arrow", 
              from = min, 
              to = max,
              title.width = 1.1,
              col.grid='grey' ,
              sizes = track.sizes,
              fontsize=11,
              main = paste0("chromosome ", chr),
              cex.main = 1,
              highlight = "darkblue") # Added highlight as feature for gene selected
  
  
}


### Create table summarising the information about the designed gRNAs ###

format_gRNA_table <- function(gn, tss_sel, genome, grnas){
  #'Returns a table with gRNA information
  #'
  #'Creates a gRNA information table for given TF-TSS. The gRNAs are included in 
  #'the table based on user input in the app.
  #'
  #'@param gn TF/gene name
  #'@param genome Genome
  #'@param grnas Table of grnas designed for the particular TF-TSS
  
  if (tss_sel != 'all'){
    
    grnas %>%
      plyranges::filter((gene == gn) & (transcript == tss_sel)) %>% 
      as_tibble() %>% 
      dplyr::mutate(predicted_score = round(predicted_score, digits = 2),
                    is_in_final_lib = factor(is_in_final_lib, levels = c(F, T))) %>% 
      dplyr::select(sgID, transcript, predicted_activity = predicted_score, off_target_stringency, gRNA_sequence, is_in_final_lib, start)
  }
  
  else {
    grnas %>%
    plyranges::filter(gene == gn) %>% 
    as_tibble() %>% 
    dplyr::mutate(predicted_score = round(predicted_score, digits = 2),
                  is_in_final_lib = factor(is_in_final_lib, levels = c(F, T))) %>% 
    dplyr::select(sgID, transcript, predicted_activity = predicted_score, off_target_stringency, gRNA_sequence, is_in_final_lib, start)}
  
  
  
}


### Plot the scores and positons of the gRNAs ###

gRNA_selection_plot<- function(tab, genome){
  #'Returns a scatter plot of selected gRNA
  #'
  #'Scatter plot of selected gRNA filtered in the pipeline. Each gene has at 
  #'most 4 selected gRNAs.
  #'
  #'@param tab Output of format_gRNA_table function from above
  #'@param genome Genome
  
  if (nrow(tab) > 0) {
    
    colors <- c("grey70", "red2")
    names(colors) <- c(F, T)
    
    p <- tab %>% 
      ggplot(aes(x=start,y=predicted_activity,
                 col=is_in_final_lib)) +
      geom_point(size=2.5) +
      xlab("Position") + ylab("Score")+
      scale_color_manual(values = colors[names(colors) %in% unique(tab$is_in_final_lib)]) +
      facet_grid(.~transcript,scales = "free")+
      theme_bw(base_size = 16)+
      labs(color = "Selected in final library?") +
      theme(axis.text.x=element_blank(),
            plot.margin = margin(5.5,
                                 case_when(length(unique(tab$transcript)) == 1 ~ 300,
                                           length(unique(tab$transcript)) == 2 ~ 150,
                                           T ~ 5.5),
                                 5.5,5.5))
    p
    
    return(p)
    
  } else {
    
    return(NULL)
    
  }
  
  
  
}

### DATA PREPARATION ###

makeGviz_geneModel_from_ucsc <- function(gtf){
  #'Gene model from UCSC annotations
  #'
  #'Returns a gene model for input gff/gtf file of UCSC annotation type. Only 
  #'protein-coding CDS are selected.
  #'
  #'@param gtf Annotation file in gtf format
  gtf %>%
    as_tibble() %>% 
    filter(type %in% c("exon","3UTR","CDS","5UTR")) %>% 
    transmute(chromosome =seqnames, 
              start, end, width, strand, 
              type = as.character(type),
              gene = gene_id,
              exon = as.character(exon_id),
              transcript = transcript_id,
              symbol = gene_name) %>% 
    group_by(exon) %>% 
    dplyr::mutate( feature = case_when("CDS"  %in% type ~ "protein_coding",
                                       "5UTR" %in% type ~ "utr",
                                       "3UTR" %in% type ~ "utr",
                                       T ~ type)) %>% 
    dplyr::select( -type) %>% ungroup
}

makeGviz_geneModel_from_np <- function(gtf){
  #'Gene model from Nanopore(pinfish) annotations
  #'
  #'Returns a gene model for input gff/gtf file of Nanopore(pinfish) annotation 
  #'type.
  #'
  #'@param gtf Annotation file in gtf format
  gtf %>% 
    as_tibble() %>% 
    transmute(chromosome = seqnames, 
              start, end, width, strand,
              transcript = transcript_id,
              symbol = gene_name)
}

makeGviz_geneModel_from_gencode<- function(gtf){
  #'Get rank from TSS_ID
  #'
  #'Splits the existing tss_id to retrieve rank as number (Internal Function)
  #'
  #'@param tss_id TSS ID
  tmp<- gtf %>% as_tibble %>% 
    filter(type %in% c("exon","UTR","CDS")) %>% 
    transmute(chromosome =seqnames, 
              start, end, width, strand, type,
              gene = gene_id,
              exon = exon_id,
              transcript = transcript_id,
              symbol = gene_name, gene_type) %>% 
    group_by(exon) %>% 
    dplyr::mutate( feature = case_when("CDS" %in% type ~ "protein_coding",
                                       "UTR" %in% type ~ "utr",
                                       T ~ gene_type)) %>% 
    dplyr::select( -gene_type, -type) %>% ungroup
  
}

makeGviz_geneModel_from_ensembl<- function(gtf){
  #'Get rank from TSS_ID
  #'
  #'Splits the existing tss_id to retrieve rank as number (Internal Function)
  #'
  #'@param tss_id TSS ID
  tmp<-gtf %>% as_tibble %>% 
    filter(type %in% c("exon","type")) %>% 
    transmute(chromosome =seqnames, 
              start, end, width, strand, 
              type = as.character(type),
              gene = gene_id,
              exon = as.character(exon_id),
              transcript = transcript_id,
              symbol = gene_name) %>% 
    group_by(exon) %>% 
    dplyr::mutate( feature = case_when("CDS" %in% type ~ "protein_coding",
                                       grepl("utr", type) ~ "utr",
                                       T ~ type)) %>% 
    dplyr::select( -type) %>% ungroup
  
}