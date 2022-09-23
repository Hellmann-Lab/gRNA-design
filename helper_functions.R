
#### TF CHARACTERISTICS TAB ####


### Plot expression histogram with the chosen expression cutoff ###

plot_expr_init <- function(expr, cutoff) {
  
  # if the expression cutoff is set to minimum, color the entire plot dark grey - percent_expr_human > cutoff (see below) would not cover everything
  if (cutoff == 0) {
    
    expr <- expr %>%
      dplyr::mutate(is_expr_human = T,
                    is_expr_cyno = T)
    
    # otherwise colour bins above the cutoff dark grey, bins below cutoff light grey
  } else{
    
    expr <- expr %>% 
      dplyr::mutate(is_expr_human = factor(percent_expr_human > cutoff, levels = c(T, F)), # > instead of >= for aesthetic reasons, otherwise bins might get 2 colours
                    is_expr_cyno = factor(percent_expr_cyno > cutoff, levels = c(T, F)))
    
  }
  
  # histogram for human
  exhist_h <- expr %>% 
    dplyr::mutate(species = "human") %>% 
    ggplot(aes(x=percent_expr_human, fill = is_expr_human)) +
    geom_histogram(col="white", binwidth = 2.5, center = 1.25) +
    scale_fill_manual(values=c("grey30", "grey80"), guide = "none")  +
    theme_bw(base_size = 14) + 
    ylab("Number of TFs expressed") +
    facet_grid(~species) +
    theme(axis.title.x = element_blank())
  
  # histogram for cynomolgus
  exhist_c <- expr %>% 
    dplyr::mutate(species = "cynomolgus") %>% 
    ggplot(aes(x=percent_expr_cyno, fill = is_expr_cyno)) +
    geom_histogram(col="white", binwidth = 2.5, center = 1.25) +
    scale_fill_manual(values=c("grey30", "grey80"), guide = "none")  +
    theme_bw(base_size = 14) + 
    ylab("Number of TFs expressed") +
    facet_grid(~species) +
    theme(axis.title = element_blank())
  
  title <- ggdraw() +
    draw_label("Percent of cells",
               size = 15,
               x = 0.5,
               hjust = 0) +
    theme(plot.margin = margin(0, 40, 10, 0),
          plot.background = element_rect(fill = "white", color = "transparent"))
  
  plot_grid(plot_grid(exhist_h,
                      exhist_c,
                      rel_widths = c(1.05, 1)),
            title,
            ncol = 1,
            rel_heights = c(1, 0.05))
  
}


### Plot expression histogram with the chosen expression cutoff and the gene of interest marked ###

plot_expr <- function(gn, expr, cutoff) {
  
  # if the expression cutoff is set to minimum, color the entire plot dark grey - percent_expr_human > cutoff (see below) would not cover everything
  if (cutoff == 0) {
    
    expr <- expr %>%
      dplyr::mutate(is_expr_human = T,
                    is_expr_cyno = T)
    
    # otherwise colour bins above the cutoff dark grey, bins below cutoff light grey
  } else{
    
    expr <- expr %>%
      dplyr::mutate(is_expr_human = factor(percent_expr_human > cutoff, levels = c(T, F)), # > instead of >= for aesthetic reasons
                    is_expr_cyno = factor(percent_expr_cyno > cutoff, levels = c(T, F)))
    
  }
  
  # histogram for human
  exhist_h <- expr %>%
    dplyr::mutate(species = "human") %>%
    ggplot(aes(x=percent_expr_human, fill = is_expr_human)) +
    geom_histogram(col="white", binwidth = 2.5, center = 1.25) +
    scale_fill_manual(values=c("grey30", "grey80"), guide = "none")  +
    theme_bw(base_size = 14) +
    ylab("Number of TFs expressed") +
    facet_grid(~species) +
    theme(axis.title.x = element_blank())
  
  # histogram for cynomolgus
  exhist_c <- expr %>%
    dplyr::mutate(species = "cynomolgus") %>%
    ggplot(aes(x=percent_expr_cyno, fill = is_expr_cyno)) +
    geom_histogram(col="white", binwidth = 2.5, center = 1.25) +
    scale_fill_manual(values=c("grey30", "grey80"), guide = "none")  +
    theme_bw(base_size = 14) +
    ylab("Number of TFs expressed") +
    facet_grid(~species) +
    theme(axis.title = element_blank())
  
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
                      rel_widths = c(1.05, 1)),
            title,
            ncol = 1,
            rel_heights = c(1, 0.05))
  
}


### Plot robustness histogram with the chosen expression cutoff ###

plot_pres_init <- function(pres, cutoff) {
  
  # colour bins above the cutoff dark grey, bins below cutoff light grey
  pres <- pres %>%
    dplyr::mutate(is_pres_human = factor(pres_human >= cutoff, levels = c(T, F)),
                  is_pres_cyno = factor(pres_cyno >= cutoff, levels = c(T, F)))
  
  # histogram for human
  exhist_h <- pres %>%
    dplyr::mutate(species = "human") %>%
    ggplot(aes(x=pres_human, fill = is_pres_human)) +
    geom_histogram(col="white", binwidth = 0.25, center = 0.125) +
    scale_fill_manual(values=c("grey30","grey80"), guide = "none")  +
    theme_bw(base_size = 14) +
    ylab("Number of modules") +
    facet_grid(~species) +
    theme(axis.title.x = element_blank())
  
  # histogram for cynomolgus
  exhist_c <- pres %>%
    dplyr::mutate(species = "cynomolgus") %>%
    ggplot(aes(x=pres_cyno, fill = is_pres_cyno)) +
    geom_histogram(col="white", binwidth = 0.25, center = 0.125) +
    scale_fill_manual(values=c("grey30", "grey80"), guide = "none")  +
    theme_bw(base_size = 14) +
    ylab("Number of modules") +
    facet_grid(~species) +
    theme(axis.title = element_blank())
  
  title <- ggdraw() +
    draw_label("Robustness measure",
               size = 15,
               x = 0.5,
               hjust = 0) +
    theme(plot.margin = margin(0, 40, 10, 0),
          plot.background = element_rect(fill = "white", color = "transparent"))
  
  plot_grid(plot_grid(exhist_h,
                      exhist_c,
                      rel_widths = c(1.05, 1)),
            title,
            ncol = 1,
            rel_heights = c(1, 0.05))
  
}


### Plot robustness histogram with the chosen expression cutoff and the gene of interest marked ###

plot_pres <- function(gn, pres, cutoff) {
  
  # colour bins above the cutoff dark grey, bins below cutoff light grey
  pres <- pres %>% 
    dplyr::mutate(is_pres_human = factor(pres_human >= cutoff, levels = c(T, F)),
                  is_pres_cyno = factor(pres_cyno >= cutoff, levels = c(T, F)))
  
  # histogram for human
  exhist_h <- pres %>% 
    dplyr::mutate(species = "human") %>% 
    ggplot(aes(x=pres_human, fill = is_pres_human)) +
    geom_histogram(col="white", binwidth = 0.25, center = 0.125) +
    scale_fill_manual(values=c("grey30","grey80"), guide = "none")  +
    theme_bw(base_size = 14) + 
    ylab("Number of modules") +
    facet_grid(~species) +
    theme(axis.title.x = element_blank())
  
  # histogram for cynomolgus
  exhist_c <- pres %>% 
    dplyr::mutate(species = "cynomolgus") %>% 
    ggplot(aes(x=pres_cyno, fill = is_pres_cyno)) +
    geom_histogram(col="white", binwidth = 0.25, center = 0.125) +
    scale_fill_manual(values=c("grey30","grey80"), guide = "none")  +
    theme_bw(base_size = 14) + 
    ylab("Number of modules") +
    facet_grid(~species) +
    theme(axis.title = element_blank())
  
  title <- ggdraw() +
    draw_label("Robustness measure",
               size = 15,
               x = 0.5,
               hjust = 0) +
    theme(plot.margin = margin(0, 40, 10, 0),
          plot.background = element_rect(fill = "white", color = "transparent"))
  
  plot_grid(plot_grid(exhist_h +
                        # add dashed line and gene label
                        geom_segment(x = pres$pres_human[pres$hub == gn], xend = pres$pres_human[pres$hub == gn], y = max(ggplot_build(exhist_h)$data[[1]]$count)*1.02, yend = 0, linetype = "dashed", size = 0.4, color = "red") +
                        geom_label(data = pres %>% dplyr::filter(hub == gn), aes(label = hub, y = max(ggplot_build(exhist_h)$data[[1]]$count)*1.05), size = 2.5, color = "red", fill = "transparent"),
                      exhist_c +
                        # add dashed line and gene label
                        geom_segment(x = pres$pres_cyno[pres$hub == gn], xend = pres$pres_cyno[pres$hub == gn], y = max(ggplot_build(exhist_c)$data[[1]]$count)*1.02, yend = 0, linetype = "dashed", size = 0.4, color = "red") +
                        geom_label(data = pres %>% dplyr::filter(hub == gn), aes(label = hub, y = max(ggplot_build(exhist_c)$data[[1]]$count)*1.05), size = 2.5, color = "red", fill = "transparent"),
                      rel_widths = c(1.05, 1)),
            title,
            ncol = 1,
            rel_heights = c(1, 0.05))
  
}


### Plot interactive cross-species conservation scatterplot ###

plot_cons_init <- function(gn_list, cons) {
  
  p <- cons %>% 
    dplyr::filter(hub %in% gn_list) %>%
    ggplot(aes(x = size, y = -log10(p.value), text = hub)) + # set hub as the variable to show when hovering over data points
    geom_point(shape = 21, alpha = 0.5, size = 0.2, fill = "black", color = "black") +
    theme_bw(base_size = 12) +
    facet_wrap(~pres_aspect, scales = "free") +
    xlab("Regulon size") +
    theme(legend.position='none',
          axis.title.y = element_blank()) 
  
  # make the plot interactive with the help of plotly's hover tooltip option
  ggplotly(p, tooltip = "text") %>% 
    layout(yaxis = list(title = list(text = "-log<sub>10</sub><i>p</i>-value",
                                     font = list(family = "Arial",
                                                 size = 16,
                                                 color = "black"))))
  
}


### Plot interactive cross-species conservation scatterplot with the gene of interest marked ###

plot_cons <- function(gn, gn_list, cons) {
  
  p <- cons %>% 
    dplyr::filter(hub %in% gn_list) %>%
    dplyr::mutate(is_selected = hub == gn) %>% 
    ggplot(aes(x = size, y = -log10(p.value), text = hub)) + # set hub as the variable to show when hovering over data points
    geom_point(aes(fill = is_selected, color = is_selected, size = is_selected, alpha = is_selected), shape = 21) +
    scale_alpha_manual(values = c(0.5, 1), guide = "none") +
    scale_size_manual(values = c(0.2, 1), guide = "none") +
    scale_fill_manual(values = c("black", "red"), guide = "none") +
    scale_color_manual(values = c("black", "red"), guide = "none") +
    theme_bw(base_size = 12) +
    facet_wrap(~pres_aspect, scales = "free") +
    xlab("Regulon size") +
    theme(legend.position='none',
          axis.title.y = element_blank()) +
    # add gene label
    geom_text(data = cons %>% filter(hub == gn),
              aes(label = hub), nudge_y = cons %>% 
                dplyr::filter(hub %in% gn_list) %>%
                arrange(pres_aspect) %>% 
                group_by(pres_aspect) %>% 
                dplyr::summarise(nudge = (max(-log10(p.value)) - min(-log10(p.value)))*0.04) %>% 
                pull(nudge), 
              color = "red", size = 3, fontface = "bold")
  
  # make the plot interactive with the help of plotly's hover tooltip option
  ggplotly(p, tooltip = "text") %>% 
    layout(yaxis = list(title = list(text = "-log<sub>10</sub><i>p</i>-value",
                                     font = list(family = "Arial",
                                                 size = 16,
                                                 color = "black"))))
  
}


### Get full-length IC scores for all motifs of the chosen gene ###

get_IC <- function(gn, motif_IC) {
  
  motif_IC %>% 
    dplyr::filter(SYMBOL == gn) %>% 
    dplyr::mutate(text = paste0(motif_id, ": IC = ", format(IC, digits = 3))) %>% 
    pull(text) %>% 
    paste(collapse = "\n")
  
}


### Plot phastCons histogram with the chosen gene marked ###

plot_phastCons <- function(gn, gn_list, phastCons) {
  
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


### Pull out the chromosome where the gene is located ###

get_chromosome <- function(gn, genome, tss) {
  
  tss %>% 
    as_tibble() %>% 
    filter(gene_name == gn) %>% 
    pull(seqnames) %>% 
    unique()
  
}


### Get the full genomic range to plot (containing all TSS regions of the chosen gene) ###

format_gene_range<- function(tss, gn, offset=2e3){
  
  tss %>% 
    ungroup %>%  
    filter(gene_name == gn) %>% 
    as_tibble %>% 
    group_by(seqnames, strand, gene_name) %>% 
    summarise(start = min(start), end = max(end)) %>% 
    as_granges %>% 
    stretch(offset)
  
}


### Pre-format data for the Gviz plot and apply the plotting function ###

make_gRNA_gviz<- function(gn, genome, gene_models, tss, grnas, offset=1e4, horlbeck, lib){
  
  # pull out designed gRNAs for the selected gene and add variable is_in_lib
  
  ## if the user selected the Horlbeck design to be included, 3 options for is_in_lib: yes, no, horlbeck
  if (horlbeck) {
    
    gRNA_data <- grnas %>% 
      plyranges::filter( gene == gn)  %>% 
      plyranges::mutate(is_in_lib = factor(ifelse(source == "horlbeck",
                                                  "horlbeck",
                                                  ifelse(sgID %in% (lib %>% 
                                                                      dplyr::filter(source == genome) %>% 
                                                                      pull(sgID)),
                                                         "yes",
                                                         "no")),
                                           levels = c("yes", "no", "horlbeck")))
    
    ## if the user did not select the Horlbeck design to be included, 2 options for is_in_lib: yes, no
  } else {
    
    gRNA_data <- grnas %>% 
      plyranges::select(-sgID.match, source.match, predicted_score.match, off_target_stringency.match) %>% 
      plyranges::filter( gene == gn) %>% 
      plyranges::mutate(is_in_lib = factor(ifelse(sgID %in% (lib %>% 
                                                               dplyr::filter(source == genome) %>%
                                                               pull(sgID)),
                                                  "yes",
                                                  "no"),
                                           levels = c("yes", "no")))
    
  }
  
  # split data by is_in_lib
  gRNA_data_split <- split(gRNA_data, as_tibble(gRNA_data)$is_in_lib)           
  
  # pull out the TSS regions for the selected gene
  tss <- tss %>% 
    plyranges::filter(gene_name == gn)
  
  # pull out Gviz gene models for the selected gene
  
  ## if the symbol is annotated, simply filter by the symbol
  if (gn %in% as_tibble(gene_models)$symbol) {
    
    gene.models = gene_models %>% 
      plyranges::filter(symbol==gn) %>% 
      plyranges::mutate(chromosome = gsub("chr", "", chromosome))
    
    ## if the symbol is not annotated, get all transcripts that are located in the region of interest + on the strand of interest
  } else {
    
    gene.models = gene_models %>% 
      plyranges::filter(start >= min(as_tibble(tss)$start) - offset & start <= max(as_tibble(tss)$end) + offset &
                          end >= min(as_tibble(tss)$start) - offset & end <= max(as_tibble(tss)$end) + offset &
                          is.na(symbol) &
                          strand == unique(as_tibble(tss)$strand) &
                          chromosome == as.character(unique(as_tibble(tss)$seqnames))) %>% 
      plyranges::mutate(chromosome = gsub("chr", "", chromosome))
    
  }
  
  # use the plotting function to create the Gviz plot
  gRNA_plot_function(PRange = format_gene_range(tss=tss, gn=gn, offset=offset) ,
                     highlight.gr = tss,
                     gRNA_data = gRNA_data_split,
                     gen = genome,
                     gene.models = gene.models, 
                     nanoporeBAM = glue::glue("data_files/{genome}/nanopore.bam"),
                     atacBW = list.files(glue::glue("data_files/{genome}/atac"), recursive=T, full.names = T))
  
}


### Gviz plot ###

gRNA_plot_function <- function(PRange,
                               highlight.gr, 
                               gene.models,
                               nanoporeBAM,
                               atacBW,
                               gRNA_data,
                               gen) {

  # get chromosome
  chr <- as.character(unique(seqnames(PRange)))
  
  # set correct chromosome format
  options(ucscChromosomeNames = FALSE)
  
  # plot axis with genomic coordinates
  axis   <- GenomeAxisTrack(genome = gen)
  

  
  # plot GENCODE annotation
  gencode.track <- GeneRegionTrack(gene.models, 
                                   chromosome = chr,
                                   genome = gen,
                                   # display gene name for each transcript
                                   showId = TRUE, 
                                   geneSymbol = TRUE, 
                                   # if the symbol is not annotated, use the ENSEMBL ID as gene name
                                   symbol = gene.models %>% 
                                              as_tibble() %>% 
                                              rowwise() %>% 
                                              dplyr::mutate(symbol = na.omit(c(symbol, gene))[1]) %>% 
                                              pull(symbol),
                                   name= "GENCODE \nannotation",
                                   fontsize = 10, cex.title = 0.8,
                                   fill = "darkblue", col = "darkblue", col.axis = "black", col.title = "black")
    
  # plot nanopore reads and coverage
  nanoporeTrack <- AlignmentsTrack(nanoporeBAM,
                                   chromosome = chr,
                                   genome = gen,
                                   name = "Nanopore \ndata", 
                                   height=0.2, coverageHeight = 0.1, minCoverageHeight = 0,
                                   window = -1, windowSize = 100,
                                   cex.title=0.8,
                                   col.axis="black", col.title="black")

  # plot gRNA positions
  
  ## plot the 3 groups (selected in final library, not selected in final library, horlbeck) separately with different colours
  gRNA_scores <- lapply(names(gRNA_data), function(n) { 
    
    DataTrack(gRNA_data[[n]] %>% 
                plyranges::select(predicted_score) %>% 
                plyranges::mutate(strand = "+"),
              chromosome = chr, 
              genome = gen,
              col = c("green4", "red4", "black")[which(names(gRNA_data) == n)],
              name="Designed \ngRNAs",
              col.axis="black", col.title="black", cex.title=0.8 , 
              ylim = c(min(as_tibble(gRNA_data)$predicted_score)*0.9, max(as_tibble(gRNA_data)$predicted_score)*1.1))  
    
    } )
  
  ## overlay the 3 tracks
  gRNATrack <- OverlayTrack( gRNA_scores ,size=4, alpha=0.5, name="gRNA\n scores",
                             window = -1, windowSize = 100,
                             legend=T )
  # plot ATAC-seq coverage 
  
  ## colours
  annot.colors = c( "Cyno"="tan2", "Human"="tomato3")
  
  ## plot each replicate (2 human and 2 cynomolgus) separately
  atac.tracks<- lapply(atacBW , function(i) {
    
    ## get species name
    spec <- strsplit(i,"/")[[1]][4]

    DataTrack( range = i, 
               type = 'h', 
               chromosome = chr,
               ## select colour according to the species
               col = unname(annot.colors[spec]) ,
               name = paste0("ATAC-seq \n", spec),
               window = -1, windowSize = 100, genome = gen,
               col.title="black", cex.title=0.8,
               col.axis="black") 
    
    })
  
  
  # add grey boxes to mark the TSS regions
  ht <- HighlightTrack( trackList = c(gRNATrack, atac.tracks) ,
                        start = start(highlight.gr) - 500, 
                        end   = end(highlight.gr) + 500,
                        chromosome = chr, 
                        fill = "grey",
                        col  = "grey")

  # combine into a single plot
  track.list<- c(axis, gencode.track, nanoporeTrack, ht)
  plotTracks( track.list, 
              collapseTranscripts = F, shape = "arrow", 
              from = start(PRange), 
              to = end(PRange),
              title.width = 1.1,
              col.grid='grey' ,
              sizes = c(0.75, 1.5, 2, 
                        c(1,rep(1,length(atacBW)))),
              fontsize=11)

}


### Create table summarising the information about the designed gRNAs ###

format_gRNA_table <- function(gn, genome, horlbeck, tss, grnas, lib){
  
  promoters <- tss %>% 
    plyranges::filter(gene_name == gn) %>% 
    plyranges::select(tss_id) %>%
    stretch(1000) %>% 
    reduce_ranges(tss_id= paste(sort(tss_id),collapse = ",")) %>% 
    sort() %>% 
    plyranges::mutate(tss_id = factor(tss_id, levels = tss_id))
  
  
  if (horlbeck) {
    
    gRNA_data <- grnas %>% 
      plyranges::filter( gene == gn) %>% 
      plyranges::mutate(is_in_lib = factor(ifelse(source == "horlbeck",
                                              "horlbeck",
                                              ifelse(sgID %in% (lib %>% 
                                                                  dplyr::filter(source == genome) %>% 
                                                                  pull(sgID)),
                                                     "yes",
                                                     "no")),
                                       levels = c("yes", "no", "horlbeck")))
    
  } else {
    
    gRNA_data <- grnas %>% 
      plyranges::select(-sgID.match, -source.match, -predicted_score.match, -off_target_stringency.match) %>% 
      plyranges::filter( gene == gn) %>% 
      plyranges::mutate(is_in_lib = factor(ifelse(sgID %in% (lib %>% 
                                                           dplyr::filter(source == genome) %>%
                                                           pull(sgID)),
                                              "yes",
                                              "no"),
                                       levels = c("yes", "no")))
    
  }
  
   gRNA_data %>% 
     plyranges::join_overlap_inner(promoters) %>% 
     as_tibble() %>% 
     dplyr::mutate(category = factor(category),
                   predicted_score = round(predicted_score, digits = 2)) %>% 
     dplyr::select(sgID, source, tss_id, predicted_activity = predicted_score, off_target_stringency, gRNA_sequence, match_in_other_genome = category, is_in_lib, start)
               
}


### Plot the scores and positons of the gRNAs ###

gRNA_selection_plot<- function(tab, genome, horlbeck){
  
  colors <- c("green4", "red4", "black")
  names(colors) <- c("yes", "no", "horlbeck")
    
    if (genome == "macFas6" | horlbeck == F) {
      
      p <- tab %>% 
        ggplot(aes(x=start,y=predicted_activity,
                   col=is_in_lib)) +
        geom_point(size=2.5) +
        xlab("Position") + ylab("Score")+
        scale_color_manual(values = colors[names(colors) %in% unique(tab$is_in_lib)]) +
        facet_grid(.~tss_id,scales = "free")+
        theme_bw(base_size = 16)+
        labs(color = "Selected in final library?")
        theme(axis.text.x=element_blank(),
              plot.margin = margin(5.5,
                                   case_when(length(unique(tab$tss_id)) == 1 ~ 300,
                                             length(unique(tab$tss_id)) == 2 ~ 150,
                                             T ~ 5.5),
                                   5.5,5.5))
      
    } else {
      
      p <- tab %>% 
        ggplot(aes(x=start,y=predicted_activity, 
                   col=is_in_lib)) +
        geom_point(size=2,alpha=0.6) + 
        xlab("Position") + ylab("Score")+
        scale_color_manual(values= colors[names(colors) %in% unique(tab$is_in_lib)]) +
        facet_grid(.~tss_id,scales = "free")+
        theme_bw(base_size = 14)+
        labs(color = "Selected in final library?")
        theme(axis.text.x=element_blank(),
              plot.margin = margin(5.5,
                                   case_when(length(unique(tab$tss_id)) == 1 ~ 300,
                                             length(unique(tab$tss_id)) == 2 ~ 150,
                                             T ~ 5.5),
                                   5.5,5.5))
      
    }
  
  return(p)
  
}
