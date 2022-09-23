#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinyBS)
library(tidyverse)
library(GenomicRanges)
library(plyranges)
library(DT)
library(Gviz)
library(plotly)
library(universalmotif)
library(ggplot2)
library(cowplot)

### Load functions ###

source("helper_functions.R")

### Load input data ###

# list of autosomal TFs with well-supported TSSs in both species and at least 1 annotated motif
TF_list <- readRDS("TF_info/TF_list.rds")

# expression data in human and cynomolgus iPSCs
expr <- readRDS("TF_info/expr.rds") %>% 
  dplyr::filter(gene_name %in% TF_list)

# robustness/preservation of regulons based on human and cynomolgus iPSCs and early neuronal lineage
pres <- readRDS("TF_info/robustness.rds") %>% 
  dplyr::filter(hub %in% TF_list)

# cross-species conservation of regulons based on human and cynomolgus iPSCs and early neuronal lineage
cons <- readRDS("TF_info/conservation.rds") %>% 
  dplyr::filter(hub %in% TF_list)

# TFBS motif IC matrices for logos
motifs <- readRDS("TF_info/TFBS_motiflist_4logo.RDS")

# TFBS motif full-length ICs
motif_IC_sum <- readRDS("TF_info/motif_IC.RDS")

# TF protein sequence conservation (mean phastCons scores)
phastCons <- readRDS("TF_info/phastCons.rds")

# Tf family annotations based on the InterPro database
TF_families <- readRDS("TF_info/TF_families.rds")

# GO term annotations
tf2go <- readRDS("TF_info/tf2go.rds")

# primate- and human-specific paralogs based on the ENSEMBL database
paralogs <- readRDS("TF_info/paralogs.rds")

# genomes
genomes <- list.files("data_files")
all_gene_models <- lapply(genomes, function(genome) { 
  readRDS( glue::glue("data_files/{genome}/gviz_gene_models.RDS")) })
names(all_gene_models)<- genomes

# TSS genomic ranges objects
all_tss <- lapply(genomes, function(genome) { 
  readRDS( glue::glue("data_files/{genome}/tf_ipsc_tss.RDS")) })
names(all_tss)<- genomes

# designed gRNAs
all_grnas <- lapply(genomes, function(genome) {
  
  filenames <- c("tf_ipsc_grnas", "tf_ipsc_grnas_with_horlbeck")
  grnas_wwo_horlbeck <- lapply(filenames, function(filename) {
    readRDS( glue::glue("data_files/{genome}/{filename}.RDS")) 
    })
  names(grnas_wwo_horlbeck) <- c(F, T)
  return(grnas_wwo_horlbeck)
  
})
 
names(all_grnas)<- genomes


# gRNA library
lib <- readRDS("TF_info/gRNA_library.rds")

### Set initial expression and module robustness cutoffs ###

# set intial cutoffs
percent_expr_init <- 15
pres_init <- -1.5

# filter genes based on cutoffs
TF_list_filt_init <- sort(intersect(expr %>% 
                                      dplyr::filter(percent_expr_human >= percent_expr_init &
                                                      percent_expr_cyno >= percent_expr_init) %>%
                                      dplyr::pull(gene_name),
                                    pres %>% 
                                      dplyr::filter(pres_human >= pres_init &
                                                      pres_cyno >= pres_init) %>%
                                      dplyr::pull(hub)))


### Set up user interface ###

ui <- fluidPage(
  
  shiny::tags$style("ul { padding-left: 3em;}"),
  
  titlePanel("Cross-species gRNA-design"),
  br(),
  
  sidebarLayout(
    
    # side panel to set input parameters
    sidebarPanel(
      
      width = 2,
      
      # expression cutoff
      sliderInput(
        "percent_expr",
        label = "Expression cutoff (% of cells)",
        min = 0,
        max = 100,
        value = percent_expr_init,
        step = 5
      ),
      
      # module robustness cutoff
      sliderInput(
        "pres",
        label = "Module robustness cutoff (AU)",
        min = -3,
        max = 5.5,
        value = pres_init,
        step = 0.5
      ),
      
      # TF (drop-down menu with only those TFs that passed the cutoffs above)
      selectizeInput(
        "gene",
        label = paste0("Gene (n = ", length(TF_list_filt_init), " above cutoffs)"),
        choices = c("", TF_list_filt_init)
      ),
      
      # genome
      radioButtons(
        "genome",
        label = "Genome",
        choiceNames = list.files("data_files"),
        choiceValues = list.files("data_files")
      ),
      
      # Should the gRNAs from the original Horlbeck paper be shown alongside with our own design?
      radioButtons(
        "horlbeck",
        label = "Include Horlbeck design?",
        choiceNames = c("no", "yes"),
        choiceValues = c(F, T)
      ),
      
      # padding for genomic region on the Gviz plot
      sliderInput(
        "offset",
        label = "Padding for genomic region",
        min = 0,
        max = 10000,
        value = 3000,
        step = 500
      ),
      
      # Are you ready?
      actionButton("go", "Go"),
      p(),
      
      # hyperlink to GitHub READ.me
      p(
        shiny::tags$a(
          href = "https://github.com/Hellmann-Lab/gRNA-design",
          "GitHub page and detailed documentation",
          target = "_blank"
        )
      )
      
    ),
    
    # main panel with the results
    mainPanel(tabsetPanel(
      type = "tabs",
      
      # available information about the chosen TF
      tabPanel("TF characteristics",
        
        # expression in iPSCs
        fluidRow(column(
          width = 9,
          h4("Expression in iPSCs"),
          plotOutput("expr_plot", height = 400)
        )),
        
        # robustness of regulon
        fluidRow(column(
          width = 9,
          h4("Module robustness"),
          plotOutput("pres_plot", height = 400)
        )),
        
        # cross-species conservation of regulon
        fluidRow(column(
          width = 9,
          h4("Cross-species conservation"),
          plotlyOutput("cons_plot", height = 450)
        )),
        
        # TFBS motif logo
        fluidRow(column(
          width = 7,
          h4("TFBS motifs"),
          plotOutput("motifs"),
          verbatimTextOutput("IC")
        )),
        
        # protein sequence conservation (mean phastCons score)
        fluidRow(column(
          width = 6,
          h4("Protein sequence conservation"),
          plotOutput("phastCons", height = 450)
        )),
        
        # TF family annotation
        fluidRow(column(
          width = 7,
          h4("TF family"),
          dataTableOutput("TF_families")
        )),
        
        # GO term annotation
        fluidRow(column(
          width = 7,
          h4("Functional annotation"),
          dataTableOutput("GO_terms")
        )),
        
        # Primate- and human-specific paralogs
        fluidRow(column(
          width = 7,
          h4("Primate- and human-specific paralogs"),
          dataTableOutput("paralogs")
        ))
        
      ),
      
      # gRNA design
      tabPanel("gRNA design",
        
        # Gviz plot       
        fluidRow(column(width = 12,
                        br(),
                        h4("Genomic region"))),
        
        fluidRow(column(
          width = 12,
          align = "center",
          shiny::span(textOutput("chromosome"), style = "font-size: 14px")
        )),
        
        fluidRow(column(
          width = 12,
          plotOutput("gviz_plot", height = 550)
        )),
        
        # Scores and positions of the gRNAs
        fluidRow(column(
          width = 12,
          br(),
          h4("Scores & positions of the gRNAs"),
          shiny::tags$br(style = "line-height: 10px"),
          plotOutput(
            "gRNA_selection",
            height = 450,
            # Equivalent to: click = clickOpts(id = "plot_click")
            click = "plot1_click",
            brush = brushOpts(id = "plot1_brush")
          )
        )),
        
        # table containing the selected gRNAs
        fluidRow(column(
          width = 5,
          br(),
          h4("Selected gRNAs"),
          dataTableOutput("brush_info"),
          br(),
          actionButton("clear", "Clear selected")
        )),
        
        # explanations of the table entries
        fluidRow(
          column(
            width = 12,
            br(),
            h4("Guide to the gRNA table entries"),
            br(),
            p(strong("sgID:")),
            p("Unique identifier of the gRNA in the chosen genome (format: geneName_strand_genomicPosition.lengthWithPAM-tssId)."),
            br(),
            p(strong("source:")),
            p("The design where the gRNA was found."),
            shiny::tags$ul(shiny::tags$li("horlbeck: the gRNAs in", em("Horlbeck et al., 2016"), "lifted over from hg19 to hg38"),
                           shiny::tags$li("hg38 design: our own design for hg38"),
                           shiny::tags$li("mf6 design: our own design for macFas6")),
            br(),
            p(strong("tss_id:")),
            p("The transcriptional start site (TSS) region that the gRNA targets. We identified the TSSs by integrating evidence from the GENCODE annotation, long read RNA-seq and bulk ATAC-seq data as well as reciprocal best BLAT hits of the human TSS in case of the non-human primates. The TSS were ranked based on available evidence, with a lower number meaning more evidence, and named according to the ranks (format: geneName_rank). For the gRNA design, we merged TSSs of a gene that are closer than 1kb into TSS regions. We refer to these regions by the concatenated names of the merged TSSs (format: geneName_rank1,geneName_rank2)."),
            br(),
            p(strong("predicted_activity:")),
            p("Activity score calculated based on", em("Horlbeck et al., 2016.")),
            p("Model:"),
            shiny::tags$ul(shiny::tags$li("in case of source 'horlbeck': the original elastic net in", em("Horlbeck et al., 2016")),
                           shiny::tags$li("in case of source 'hg38 design' and 'mf6 design': a modified version of this elastic net trained with only one type of accessibility data (human iPSC ATAC-seq data as a replacement of the original DNase data)")),
            p("The scores are within the range [0, 1] with a higher score meaning a higher activity/more efficient knock-down. The strongest predictor of the activity is the position relative to the TSS, including both distance from the TSS and avoidance of canonical nucleosome-occupied regions."),
            br(),
            p(strong("off_target_stringency:")),
            p("The level of the most stringent off-target filter the gRNA passed, with a lower level meaning a more stringent filter/less off-target effects."),
            p("Filters applied:"),
            shiny::tags$ul(shiny::tags$li("31_nearTSS: passed if there is no off-target site within the TSS regions with a mismatch threshold of 31"),
                           shiny::tags$li("21_genome: passed if there is no off-target site within the whole genome with a mismatch threshold of 21"),
                           shiny::tags$li("31_2_nearTSS: passed if there is <=1 off-target site within the TSS regions with a mismatch threshold of 31"),
                           shiny::tags$li("31_3_nearTSS: passed if there are <=2 off-target sites within the TSS regions with a mismatch threshold of 31")),
            p("Levels:"),
            shiny::tags$ul(shiny::tags$li("0: the gRNA passed the filters 31_nearTSS & 21_genome"),
                           shiny::tags$li("1: the gRNA passed the filter 31_nearTSS but not 21_genome"),
                           shiny::tags$li("2: the gRNA passed the filter 21_genome but not 31_nearTSS"),
                           shiny::tags$li("3: the gRNA passed the filter 31_2_nearTSS but not 31_nearTSS or 21_genome"),
                           shiny::tags$li("4: the gRNA passed the filter 31_3_nearTSS but not 31_nearTSS, 21_genome or 31_2_nearTSS")),
            p("If a gRNA does not pass any of these filters, it is discarded from the design output."),
            br(),
            p(strong("gRNA sequence:")),
            p("The 20-nt long sequence of the gRNA. In case of source 'horlbeck', this means the original hg19 sequence which occassionally differs from the sequence at the position where the gRNA is displayed in hg38 genome space."),
            br(),
            p(strong("match_in_other_genome:")),
            p("Category based on the matching between the gRNAs designed for hg38 (horlbeck + own design) and the gRNAs designed for mf6."),
            shiny::tags$ul(shiny::tags$li("zero_mismatch: the gRNA has a perfectly matching counterpart in the other library"),
                           shiny::tags$li("one_mismatch: the gRNA has a counterpart with 1 mismatch in the other library"),
                           shiny::tags$li("no_cyno_gRNA/no_human_gRNA: the gRNA does not have a match with ≤1 mismatch in the other library")),
            br(),
            p(strong("is_in_lib:")),
            p("'Yes' if the gRNA is in the final CRISPRi library, 'No' if it is not. We designed two species-specific libraries, one for human and one for cynomolgus, each of them contains 4 gRNAs per TSS region. Both effectiveness (high predicted activity scores) and comparability (similar predicted activity scores in the two species) was taken into account when compiling the libraries.")
        ))
      )
    ))
  )
)


### Define server logic ###

server <- function(input, output) {
  
  # expression histogram (with the chosen cutoff marked)
  output$expr_plot <- renderPlot({
    
    plot_expr_init(expr = expr,
                   cutoff = input$percent_expr)
    
  })
  
  # module robustness histogram (with the chosen cutoff marked)
  output$pres_plot <- renderPlot({
    
    plot_pres_init(pres = pres,
                   cutoff = input$pres)
    
  })
  
  # interactive cross-species conservation scatterplot (showing only the genes that pass the cutoffs above)
  output$cons_plot <- renderPlotly({
    
    plot_cons_init(gn_list = intersect(expr %>% 
                                         dplyr::filter(percent_expr_human >= input$percent_expr & percent_expr_cyno >= input$percent_expr) %>%
                                         dplyr::pull(gene_name),
                                       pres %>% 
                                         dplyr::filter(pres_human >= input$pres & pres_cyno >= input$pres) %>%
                                         dplyr::pull(hub)),
                   cons = cons)
  })
  
  # adjust the drop-down menu according to the selected cutoffs
  observe({
    
    # if the cutoffs are set to the minimum, all genes in the full TF list are included
    if (input$percent_expr == 0 & input$pres == -3) {
      
      TF_list_filt <- TF_list 
    
    # otherwise only those TFs are included that pass the cutoffs (we need to have expression data about them to check whether they pass the cutoffs -> all TFs that don't appear in the count matrix are also removed)  
    } else{
      
      TF_list_filt <- intersect(expr %>% 
                                    dplyr::filter(percent_expr_human >= input$percent_expr & percent_expr_cyno >= input$percent_expr) %>%
                                    dplyr::pull(gene_name),
                                  pres %>% 
                                    dplyr::filter(pres_human >= input$pres & pres_cyno >= input$pres) %>%
                                    dplyr::pull(hub))
      
    }
    
    updateSelectizeInput(inputId = "gene", 
                         label = paste0("Gene (n = ", length(TF_list_filt), " above cutoffs)"),
                         choices = TF_list_filt,
                         selected = input$gene)
  })
  
  # store the user's last action (will be used to remove outdated results as soon as the user changes one of the parameters)
  values <- reactiveValues()
  values$lastAction <- NULL
  observe({if (input$percent_expr != 0 | input$pres != 0 | input$gene != 0 | input$genome != 0 | input$offset != 0) {values$lastAction <- "data"}})
  observe({if (input$go != 0) {values$lastAction <- "plot"}})
  
  # upon pressing "Go"...
  observeEvent( input$go, {  
    
    # ...show expression histogram (with the selected TF marked) 
    output$expr_plot <- renderPlot({
      
      # remove the plot if the user changes some of the parameters afterwards or if there is no expression data available about the gene
      if (is.null(values$lastAction)) return(NULL)
      if (values$lastAction=="data") return(NULL) 
      if (input$gene == "") return(NULL)
      if (!(input$gene %in% expr$gene_name)) return(NULL)
      
      # show modal dialog while waiting
      showModal(modalDialog("Creating expression plot", footer = NULL))
      
      p1 <- plot_expr(gn = input$gene, 
                      expr = expr,
                      cutoff = input$percent_expr)
      
      removeModal()
      
      p1
      
    })
    
    # ...show module robustness histogram (with the selected TF marked) 
    output$pres_plot <- renderPlot({
      
      # remove the plot if the user changes some of the parameters afterwards or if there is no expression data available about the gene
      if (is.null(values$lastAction)) return(NULL)
      if (values$lastAction=="data") return(NULL) 
      if (input$gene == "") return(NULL)
      if (!(input$gene %in% expr$gene_name)) return(NULL)
      
      # show modal dialog while waiting
      showModal(modalDialog("Creating robustness plot", footer = NULL))
      
      p2 <- plot_pres( gn = input$gene, 
                       pres = pres,
                       cutoff = input$pres)
      
      removeModal()
      
      p2
      
    })
    
    
    # ...show the cross-species conservation plot (with the selected TF marked) 
    output$cons_plot <- renderPlotly({
      
      # remove the plot if the user changes some of the parameters afterwards or if there is no expression data available about the gene
      if (is.null(values$lastAction)) return(NULL)
      if (values$lastAction=="data") return(NULL) 
      if (input$gene == "") return(NULL)
      if (!(input$gene %in% expr$gene_name)) return(NULL)
      
      plot_cons( gn = input$gene, 
                 gn_list = intersect(expr %>% 
                                       dplyr::filter(percent_expr_human >= input$percent_expr & percent_expr_cyno >= input$percent_expr) %>%
                                       dplyr::pull(gene_name),
                                     pres %>% 
                                       dplyr::filter(pres_human >= input$pres & pres_cyno >= input$pres) %>%
                                       dplyr::pull(hub)),
                 cons = cons)
    })
    
    # ...show the TFBS motif logo
    output$motifs <- renderPlot({
      
      # remove the plot if the user changes some of the parameters afterwards
      if (is.null(values$lastAction)) return(NULL)
      if (values$lastAction=="data") return(NULL) 
      if (input$gene == "") return(NULL)
      
      universalmotif::view_motifs( motifs[[input$gene]])
      
    })
    
    # ...show the full-length IC of the motifs
    output$IC <- renderText({
      
      get_IC(input$gene, motif_IC_sum)
      
    })
    
    # show the protein sequence conservation histogram (with the selected TF marked)
    output$phastCons <- renderPlot({
      
      # remove the plot if the user changes some of the parameters afterwards
      if (is.null(values$lastAction)) return(NULL)
      if (values$lastAction=="data") return(NULL) 
      if (input$gene == "") return(NULL)
      
      # include only genes in the histogram that pass the expression and module robustness cutoffs
      if (input$percent_expr == 0 & input$pres == -3) {
        
        TF_list_filt <- TF_list 
        
      } else{
        
        TF_list_filt <- intersect(expr %>% 
                                      dplyr::filter(percent_expr_human >= input$percent_expr & 
                                                      percent_expr_cyno >= input$percent_expr) %>%
                                      dplyr::pull(gene_name),
                                    pres %>% 
                                      dplyr::filter(pres_human >= input$pres & 
                                                      pres_cyno >= input$pres) %>%
                                      dplyr::pull(hub))
        
      }
      
      plot_phastCons(input$gene, 
                     gn_list = TF_list_filt,
                     phastCons = phastCons %>% dplyr::filter(gene_name %in% TF_list_filt))
      
    })
    
    # ...show table containing the TF family annotations
    output$TF_families <- DT::renderDataTable({
      
      # remove the table if the user changes some of the parameters afterwards
      if (is.null(values$lastAction)) return(NULL)
      if (values$lastAction=="data") return(NULL) 
      if (input$gene == "") return(NULL)
      
      DT::datatable(TF_families[[input$gene]], rownames = F)
      
    })
    
    # ...show table containing the GO term annotations
    output$GO_terms <- DT::renderDataTable({
      
      # remove the table if the user changes some of the parameters afterwards
      if (is.null(values$lastAction)) return(NULL)
      if (values$lastAction=="data") return(NULL) 
      if (input$gene == "") return(NULL)
      
      DT::datatable(tf2go[[input$gene]], rownames = F)
      
    })
    
    # ...show table containing the primate- and human-specific paralogs
    output$paralogs <- DT::renderDataTable({
      
      # remove the table if the user changes some of the parameters afterwards
      if (is.null(values$lastAction)) return(NULL)
      if (values$lastAction=="data") return(NULL) 
      if (input$gene == "") return(NULL)
      
      DT::datatable(paralogs[[input$gene]], rownames = F)
      
    })
    
    # ...show the chromosome where the selected gene is located
    output$chromosome <- renderText({
      
      # remove if the user changes some of the parameters afterwards
      if (is.null(values$lastAction)) return(NULL)
      if (values$lastAction=="data") return(NULL) 
      if (input$gene == "") return(NULL)
      
      paste("chromosome", 
            get_chromosome( gn = input$gene, 
                            genome = input$genome,
                            tss = all_tss[[input$genome]]))
      
    })
    
    #...show Gviz plot of the genomic region around the TSS of the selected gene (GENCODE annotation, nanopore data, positions of the designed gRNAs, human and cynomolgus ATAC-seq data)
    output$gviz_plot <- renderPlot({
      
      # remove the plot if the user changes some of the parameters afterwards
      if (is.null(values$lastAction)) return(NULL)
      if (values$lastAction=="data") return(NULL) 
      if (input$gene == "") return(NULL)
      
      # show modal dialog while waiting
      showModal(modalDialog("Creating Gviz plot", footer=NULL))
      
      p3 <- make_gRNA_gviz( gn = input$gene, 
                            genome =input$genome ,
                            offset=input$offset,
                            gene_models = all_gene_models[[input$genome]],
                            tss = all_tss[[input$genome]],
                            grnas = all_grnas[[input$genome]][[input$horlbeck]],
                            horlbeck = input$horlbeck,
                            lib = lib)
      
      removeModal()
      
      p3
      
    })
    
    # show detailed plot about the positions and predicted activity scores of the gRNAs
    output$gRNA_selection <- renderPlot({
      
      # remove the plot if the user changes some of the parameters afterwards
      if (is.null(values$lastAction)) return(NULL)
      if (values$lastAction=="data") return(NULL) 
      if (input$gene == "") return(NULL)
      
      showModal(modalDialog("Creating gRNA selection plot", footer=NULL))
      
      p4 <- gRNA_selection_plot( format_gRNA_table(gn=input$gene, 
                                                   genome=input$genome,
                                                   horlbeck = input$horlbeck,
                                                   tss = all_tss[[input$genome]],
                                                   grnas = all_grnas[[input$genome]][[input$horlbeck]],
                                                   lib = lib),
                                 genome = input$genome,
                                 horlbeck = input$horlbeck)
      
      removeModal()
      
      p4
      
    })
  
  } )
  
  # initialize table containing the selected gRNAs
  collect_input <- reactiveValues(selected_gRNAs = tibble())

  # add gRNAs that are near the location where the user clicks
  observeEvent(input$plot1_click, {
    
    collect_input$selected_gRNAs <- bind_rows(collect_input$selected_gRNAs,
                                              nearPoints(format_gRNA_table(gn = input$gene,
                                                                           genome = input$genome,
                                                                           horlbeck = input$horlbeck,
                                                                           tss = all_tss[[input$genome]],
                                                                           grnas = all_grnas[[input$genome]][[input$horlbeck]],
                                                                           lib = lib),
                                                         input$plot1_click)) %>%
      # make sure there are no duplicate entries
      distinct_all() %>% 
      dplyr::select(-start)
    
  })
  
  # add gRNAs that are inside the rectangle that the user draws
  observeEvent(input$plot1_brush, {
    
    collect_input$selected_gRNAs <-bind_rows(collect_input$selected_gRNAs,
                                             brushedPoints(format_gRNA_table(gn = input$gene,
                                                                             genome = input$genome,
                                                                             horlbeck = input$horlbeck,
                                                                             tss = all_tss[[input$genome]],
                                                                             grnas = all_grnas[[input$genome]][[input$horlbeck]],
                                                                             lib = lib),
                                             input$plot1_brush)) %>%
      # make sure there are no duplicate entries
      distinct_all() %>% 
      dplyr::select(-start)
    
  })
  
  # clear table if the user clicks "Clear"
  observeEvent( input$clear, {  
    
    collect_input$selected_gRNAs <- tibble()
    
  } )
  
  # show table containing the selected gRNAs
  output$brush_info <- DT::renderDataTable({
    
    DT::datatable(collect_input$selected_gRNAs, 
                  filter = "top",
                  class = 'cell-border stripe',
                  rownames = FALSE,extensions = 'Buttons',
                  options = list(dom = 'Bfrtip', buttons = c('csv', 'excel')))
    
  })
  
}


### run the application ###

shinyApp(ui = ui, server = server)
