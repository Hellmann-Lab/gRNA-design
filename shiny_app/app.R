#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#


### (if there's an error message with Rle/match, the package GenomicRanges needs to be unloaded -> best to restart the R session) ###

library(shiny)
library(shinyBS)
library(tidyverse)
library(plyranges)
library(DT)
library(Gviz)
library(plotly)
library(universalmotif)
library(cowplot)


# Setting list of options for shiny app
# option_list_shiny = list(make_option(c("--helper_functions_path"),
#                                type = "character",
#                                help = "Path to sourcecode funcions",
#                                metavar = "character",
#                                default = '/data/share/htp/perturb-seq/gRNA_design_workflow/shiny_app/helper_functions.R'),
#                    make_option(c("--paramter_path"),
#                                type = "character",
#                                help = "Path to Input paramters YAML file",
#                                metavar = "character",
#                                default = '/data/share/htp/perturb-seq/gRNA_design_workflow/shiny_app/input_chimp.yaml')
# )
# 
# # options parser
# opt_parser_shiny = OptionParser(option_list = option_list_shiny)
# opt_shiny = parse_args(opt_parser_shiny)

### Load functions ###
#source(opt_shiny$helper_functions_path)

### Read user paramters ###
#input_parameters <- read_yaml(opt_shiny$paramter_path)

### Set working directory ### 
#setwd(input_parameters$shiny_working_dir)

### Load input data ###

# list of autosomal TFs with well-supported TSSs in both species and at least 1 annotated motif
#TF_list <- readRDS(input_parameters$TF_list_path)

# genomes
genomes <- list.files("data_files")

# gene models
all_gene_models <- lapply(genomes, function(genome) {
  #'List of gene models
  #'
  #'Returns a list of gene models by reading the given path for each species
  #'
  #'@param genome Genome
  
  gene_model_names_spec <- gsub(".rds", "", list.files(glue::glue("data_files/{genome}/gene_models")))
  
  gene_models_spec <- lapply(gene_model_names_spec, function(gene_model_name) {
    
    readRDS(glue::glue("data_files/{genome}/gene_models/{gene_model_name}.rds"))
    
  })
  
  names(gene_models_spec) <- gene_model_names_spec
  
  return(gene_models_spec)
  
 })
names(all_gene_models) <- genomes

# TSS genomic ranges objects
all_tss <- lapply(genomes, function(genome) {
  #'TSS for each species
  #'
  #'Returns a list of TSS for each species by reading the given path
  #'
  #'@param genome Genome
  
  tss_names_spec <- gsub(".rds", "", list.files(glue::glue("data_files/{genome}/tss")))
  
  tss_spec <- lapply(tss_names_spec, function(tss_name) {
    
    readRDS(glue::glue("data_files/{genome}/tss/{tss_name}.rds"))
    
  })
  
  names(tss_spec) <- tss_names_spec
  
  return(tss_spec)
  
})
names(all_tss) <- genomes

# designed gRNAs
all_grnas <- lapply(genomes, function(genome) {
  #'gRNAS for each species
  #'
  #'Returns a list of gRNAS for each species by reading the given path
  #'
  #'@param genome Genome

  readRDS(glue::glue("data_files/{genome}/gRNAs/gRNAs.rds"))

})
names(all_grnas)<- genomes

# ATAC
all_atacBW <- lapply(genomes, function(genome) {
  #'ATAC-Seq peaks for each species
  #'
  #'Returns a list of ATAC-Seq peaks for each species by reading the given path
  #'
  #'@param genome Genome
  
  atac_spec <- list.files(glue::glue("data_files/{genome}/atac"), recursive=T, full.names = T)
  names_atac_spec <- str_split(atac_spec, "\\/|\\.", simplify = T)[,5]
  names(atac_spec) <- names_atac_spec
  
  return(atac_spec)
  
})
names(all_atacBW) <- genomes

# nanopore reads
all_nanoporeBAM <- lapply(genomes, function(genome) {
  #'Nanopore expression for each species
  #'
  #'Returns a list of Nanopore expression data for each species by reading the given path
  #'
  #'@param genome Genome
  
  if (file.exists("data_files/{genome}/nanopore_reads/nanopore.bam")){
    glue::glue("data_files/{genome}/nanopore_reads/nanopore.bam")
  } else {NA}
  
  
})
names(all_nanoporeBAM) <- genomes


### Set up user interface ###

ui <- fluidPage(
  
  shiny::tags$style("ul { padding-left: 3em;}"),
  
  titlePanel("gRNA-design"),
  br(),
  
  sidebarLayout(
    
    # side panel to set input parameters
    sidebarPanel(
      
      width = 2,
      
      # genome
      radioButtons(
        "genome",
        label = "Genome",
        choiceNames = list.files("data_files"),
        choiceValues = list.files("data_files"),
        selected = input_parameters$genome
      ),
      
      # TF (drop-down menu with only those TFs that passed the cutoffs above)
      selectizeInput(
        "gene",
        label = paste0("Gene (n = ", length(TF_list), " above cutoffs)"),
        choices = c("", TF_list)
      ),
      
      # TSS (TSS n selected TF)
      selectizeInput(
        "tss_sel",
        label = paste0("Select TSS in gene"),
        choices = c("all")
      ),
      
      
      # show the whole gene body or zoom in on TSS?
      radioButtons(
        "gviz_focus",
        label = "Show on Gviz",
        choiceNames = c("entire gene body", "TSS only"),
        choiceValues = c("gene_body", "tss"),
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
      
      # gRNA design
      tabPanel("gRNA design",
        
        # Gviz plot       
        fluidRow(column(width = 12,
                        br(),
                        h4("Genomic region"))),
        
        fluidRow(column(
          width = 12,
          plotOutput("gviz_plot", height = 900, width = "100%")
        )),
        
        fluidRow(column(width = 12,
                        br(),
                        actionButton('save_gviz_plot', "Save Plot"))),
        
        # Scores and positions of the gRNAs
        fluidRow(column(
          width = 12,
          br(),
          h4("Scores & positions of the gRNAs"),
          shiny::tags$br(style = "line-height: 10px"),
          plotOutput(
            "gRNA_selection",
            height = 450, 
            width = "100%",
            # Equivalent to: click = clickOpts(id = "plot_click")
            click = "plot1_click",
            brush = brushOpts(id = "plot1_brush")
          )
        )),
        
        fluidRow(column(width = 12,
                        br(),
                        actionButton('save_grna_plot', "Save Plot"))),
        
        # table containing the selected gRNAs
        fluidRow(column(
          width = 12,
          height = 3,
          br(),
          h4("Selected gRNAs"),
          dataTableOutput("brush_info",  width = "100%", height = "3%"),
          br(),
          actionButton("clear", "Clear selected")
        )),
        
        fluidRow(column(width = 12,
                        br(),
                        actionButton('save_grna_table', "Save as Plot"))),
        
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

server <- function(input, output, session) {
  
  
  # store the user's last action (will be used to remove outdated results as soon as the user changes one of the parameters)
  values <- reactiveValues()
  values$lastAction <- NULL
  observe({if ( input$gene != 0 | input$genome != 0 | input$offset != 0 | input$gviz_focus != 0) {values$lastAction <- "data"
  
  # Update value of TSS selection
  updateSelectInput(session, 'tss_sel', choices = c(all_tss[[input$genome]]$final %>% as_tibble() %>% 
                         filter(gene_name == input$gene) %>% 
                         pull(tss_id), 'all'), selected = 'all')}})
  observe({if (input$go != 0) {values$lastAction <- "plot"}})
  
  # upon pressing "Go"...
  observeEvent( input$go, {  
    
    # .p3 <- reactive(make_gRNA_gviz( gn = input$gene, 
    #                                 tss_sel = input$tss_sel,
    #                                 gen = input$genome,
    #                                 offset = input$offset,
    #                                 gene_models = all_gene_models[[input$genome]],
    #                                 tss = all_tss[[input$genome]],
    #                                 grnas = all_grnas[[input$genome]],
    #                                 atacBW = all_atacBW[[input$genome]],
    #                                 nanoporeBAM = all_nanoporeBAM[[input$genome]],
    #                                 gviz_focus = input$gviz_focus))
    
    #...show Gviz plot of the genomic region around the TSS of the selected gene (GENCODE annotation, nanopore data, positions of the designed gRNAs, human and cynomolgus ATAC-seq data)
    output$gviz_plot <- renderPlot({
      
      # remove the plot if the user changes some of the parameters afterwards
      if (is.null(values$lastAction)) return(NULL)
      if (values$lastAction=="data") return(NULL) 
      if (input$gene == "") return(NULL)
      if (input$tss_sel == "") return(NULL)
      
      # show modal dialog while waiting
      showModal(modalDialog("Creating Gviz plot", footer=NULL))
      
      p3 <- make_gRNA_gviz( gn = input$gene, 
                            tss_sel = input$tss_sel,
                            gen = input$genome,
                            offset = input$offset,
                            gene_models = all_gene_models[[input$genome]],
                            tss = all_tss[[input$genome]],
                            grnas = all_grnas[[input$genome]],
                            atacBW = all_atacBW[[input$genome]],
                            nanoporeBAM = all_nanoporeBAM[[input$genome]],
                            gviz_focus = input$gviz_focus)
      
      removeModal()
      
      p3
      
    })
    
    observeEvent(input$save_gviz_plot, {
      
      # Save the Gviz as a SVG file
      svg(file="gviz.svg", width=10, height=7)
      p3 <- make_gRNA_gviz( gn = input$gene, 
                            tss_sel = input$tss_sel,
                            gen = input$genome,
                            offset = input$offset,
                            gene_models = all_gene_models[[input$genome]],
                            tss = all_tss[[input$genome]],
                            grnas = all_grnas[[input$genome]],
                            atacBW = all_atacBW[[input$genome]],
                            nanoporeBAM = all_nanoporeBAM[[input$genome]],
                            gviz_focus = input$gviz_focus)
      dev.off()

    })

    
    # show detailed plot about the positions and predicted activity scores of the gRNAs
    output$gRNA_selection <- renderPlot({
      
      # remove the plot if the user changes some of the parameters afterwards
      if (is.null(values$lastAction)) return(NULL)
      if (values$lastAction=="data") return(NULL) 
      if (input$gene == "") return(NULL)
      if (input$tss_sel == "") return(NULL)
      
      showModal(modalDialog("Creating gRNA selection plot", footer=NULL))
      
      p4 <- gRNA_selection_plot(format_gRNA_table(gn = input$gene, 
                                                  tss_sel = input$tss_sel,
                                                  genome = input$genome,
                                                  grnas = all_grnas[[input$genome]]),
                                 genome = input$genome)
      
      removeModal()
      
      p4
      
    })
    
    observeEvent(input$save_grna_plot, {
      
      # Save the ggplot as a SVG file
      #svg(file="grna.svg", width=10, height=5)
      p4 <- gRNA_selection_plot(format_gRNA_table(gn = input$gene, 
                                                  tss_sel = input$tss_sel,
                                                  genome = input$genome,
                                                  grnas = all_grnas[[input$genome]]),
                                genome = input$genome)
      
      ggsave('grna.svg', p4, width=15, height=6, units = "in", dpi = 300)
      #dev.off()
      
    })
  
  } )
  
  # initialize table containing the selected gRNAs
  collect_input <- reactiveValues(selected_gRNAs = tibble())
  
  # add gRNAs that are near the location where the user clicks
  observeEvent(input$plot1_click, {
    
    collect_input$selected_gRNAs <- bind_rows(collect_input$selected_gRNAs,
                                              nearPoints(format_gRNA_table(gn = input$gene,
                                                                           tss_sel = input$tss_sel,
                                                                           genome = input$genome,
                                                                           grnas = all_grnas[[input$genome]]),
                                                         input$plot1_click)) 
    
  })
  
  # add gRNAs that are inside the rectangle that the user draws
  observeEvent(input$plot1_brush, {
    
    collect_input$selected_gRNAs <-bind_rows(collect_input$selected_gRNAs,
                                             brushedPoints(format_gRNA_table(gn = input$gene,
                                                                             tss_sel = input$tss_sel,
                                                                             genome = input$genome,
                                                                             grnas = all_grnas[[input$genome]]),
                                                           input$plot1_brush))
    
  })
  
  # clear table if the user clicks "Clear"
  observeEvent( input$clear, {  
    
    collect_input$selected_gRNAs <- tibble()
    
  } )
  
  # show table containing the selected gRNAs
  output$brush_info <- DT::renderDataTable({
    
    DT::datatable(collect_input$selected_gRNAs %>% distinct() %>% dplyr::select(-start), 
                  filter = "top",
                  class = 'cell-border stripe',
                  rownames = FALSE,extensions = 'Buttons',
                  options = list(dom = 'Bfrtip', buttons = c('csv', 'excel'), 
                                 width = 30,
                                 height = 3,
                                 initComplete = JS(
                                   "function(settings, json) {",
                                   
                                   "$(this.api().table().header()).css({'font-size': '140%', 'text-align': 'center'});",
                                   "$(this.api().table().container()).css({'font-size': '130%', 'text-align': 'center'});", 
                                   "}")))
    
  })
  
  observeEvent(input$save_grna_table, {
    
    # Save the table as a PDF
    dt_table <- DT::datatable(collect_input$selected_gRNAs %>% distinct() %>% dplyr::select(-start), 
                              filter = "top",
                              class = 'cell-border stripe',
                              rownames = FALSE,extensions = 'Buttons',
                              options = list(dom = 'Bfrtip', buttons = c('csv', 'excel'), 
                                             width = 30,
                                             height = 3,
                                             initComplete = JS(
                                               "function(settings, json) {",
                                               
                                               "$(this.api().table().header()).css({'font-size': '140%', 'text-align': 'center'});",
                                               "$(this.api().table().container()).css({'font-size': '130%', 'text-align': 'center'});", 
                                               "}")))
    
    html <- "grna_table.html"
    saveWidget(dt_table, html)
    # webshot::webshot(html,
    #                  file = "grna_table.png")
    
    })
  
}

### run the application ###
shinyApp(ui = ui, server = server)
