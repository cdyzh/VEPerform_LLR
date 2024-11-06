# ui.R
library(shiny)
library(shinyjs)
library(dplyr)
library(yogiroc)
library(rsconnect)
library(rmarkdown)
library(knitr)
library(tinytex)
library(httr)
library(jsonlite)
library(shinycssloaders)
library(DT)

# Define UI for the application
ui <- fluidPage(
  
  useShinyjs(),  # Initialize shinyjs for resetting files after error
  
  # Introductory Page with Navigation
  navbarPage(
    title = "VEPerform",
    id = "navbar",  # Set an ID for the navbarPage
    
    # Introduction Tab
    tabPanel("Introduction",
             fluidRow(
               column(6, offset = 3,
                      h2("Welcome to VEPerform"),
                      p("This tool allows you to evaluate the performance of variant effect predictors for your favorite genes."),
                      p("You can either use our pre-existing dataset, fetch data live, or upload your own dataset."),
                      
                      # Radio buttons for selecting between using an existing dataset or uploading a new one
                      # FIXME: used html formatting for boldfacing this. Make front-end consistent later. 
                      radioButtons("intro_data_source", "Please select an option:",
                                   choiceNames = list(
                                     tags$span(tags$b("Use pre-existing set of reference variants from ClinVar:"), tags$br(), "the quickest way to generate a PRC"),
                                     tags$span(tags$b("Advanced:"), tags$br(), "use your custom reference set and/or fetch the most up-to-date VEP scores and annotations from OpenCRAVAT")
                                     # tags$span(tags$b("Customize and upload your own reference set:"), tags$br(), "if you have your own data and just want to use the PRC generating tool")
                                   ),
                                   choiceValues = list("existing", "fetch"),
                                   selected = "existing"),
                      
                      # Link to API
                      p("You can also interact with the VEPerform API ", 
                        a("here", href = "http://127.0.0.1:4697/", target = "_blank")),
                      
                      # Action button to proceed
                      actionButton("proceed_button", "Proceed", class = "btn-success")
               )
             )
    ),
    
    # Main App (Use Pre-existing Dataset) Tab
    tabPanel("Main App",
             fluidPage(
               # Introductory Text
               tags$h3("Welcome to the Main App"),
               tags$p("This page allows you to evaluate the performance of variant effect predictors for specific genes. Select a gene to start.
                   You can use our pre-existing dataset or upload your own dataset to select a subset of variants for analysis. 
                   After selecting your gene and the scores you’d like to include, generate a Precision-Recall Curve (PRC) to 
                   visualize prediction performance. You can also view and deselect specific variants before plotting. "),
               
               sidebarLayout(
                 sidebarPanel(
                   # Radio buttons to select between uploading a dataset or using an existing one
                   radioButtons("data_source", "Choose Dataset:",
                                choices = list("Use Existing Full Dataset" = "existing", 
                                               "Select Subset of Variants" = "upload"),
                                selected = "existing"),
                   
                   # Conditional panel that shows upload options if the user chooses to upload their own dataset
                   conditionalPanel(
                     condition = "input.data_source == 'upload'",
                     helpText("Using the outputted CSV, delete the rows for variants you would like to remove and keep only the gene and hgvs_pro columns." ),
                     fileInput("file_gene_variant", "Upload Gene and HGVS_Pro CSV File", 
                               accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv"))
                     # actionButton("upload_guide_existing", "Upload Guide", class = "btn-info")
                   ),
                   
                   selectizeInput("Main_gene", "Select Gene Name:", choices = NULL, options = list(maxOptions = 1000)),
                   checkboxInput("Main_common_variant_filter", "Exclude Common Variants (gnomAD AF > 0.005)", value = TRUE),
                   checkboxGroupInput("Main_scores", "Select Scores to Include:",
                                      choices = list("VARITY", "REVEL", "AlphaMissense"),
                                      selected = c("VARITY", "REVEL", "AlphaMissense")),
                   actionButton("Main_plotButton", "Generate PRC Plot"),
                   downloadButton("Main_downloadPlotPNG", "Download PRC Plot as PNG"),
                   downloadButton("Main_downloadPlotPDF", "Download PRC Plot and Metadata as PDF"), 
                   downloadButton("downloadCSV", "Download Variants Used as CSV")
                 ),
                 mainPanel(
                   plotOutput("Main_PRCPlot", width = "600px", height = "600px"),
                   textOutput("Main_ErrorText")
                 )
               )
             )
    ),
    
    tabPanel("Advanced", 
             fluidPage(
               # Introductory Text
               tags$h3("Customizable PRC"),
               tags$p("This option offers maximal input customizability. In this section, you can use your own reference set to generate a PRC. 
               Ensure that your file aligns with the required specifications. You can also supplement your reference set with VEP scores from 
               OpenCRAVAT, or create a new reference set using OpenCRAVAT. 
              (Note: if using OpenCRAVAT, make sure you have this information for each variant: chromosome, position, reference base and alternate base.)"),
               
               sidebarLayout(
                 sidebarPanel(
                   radioButtons("input_type", "Would you like to:",
                                choices = list("Use my own reference set only" = "own",
                                               "Fetch VEP scores from OpenCRAVAT / add to my own reference set" = "fetch")),
                   
                   # File input for Chromosome, Position, Reference_Base, Alternate_Base
                   conditionalPanel(
                     condition = "input.input_type == 'own'",
                     fileInput("file_full", "Upload Full Reference Set CSV File",
                               accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
                     actionButton("upload_guide", "Upload Guide", class = "btn-info")
                     
                   ),
                   
                   
                   conditionalPanel(
                     condition = "input.input_type == 'fetch'",
                     fileInput("file_fetch", "Upload CSV File. Ensure it contains columns chrom, pos, ref_base, and alt_base",
                               accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
                     helpText("Click 'Fetch VEP Data' after uploading your file to fetch data from OpenCRAVAT" ),
                     
                     checkboxGroupInput("Fetch_scores", "Select what to fetch from OpenCRAVAT:",
                                        choices = c(
                                          "ClinVar" = "clinvar", 
                                          "GnomAD" = "gnomad", 
                                          "VARITY_R" = "varity_r", 
                                          "REVEL" = "revel", 
                                          "AlphaMissense" = "alphamissense", 
                                          "VARITY_ER" = "varity_er", 
                                          "MaveDB" = "mavedb", 
                                          "PROVEAN" = "provean",
                                          "SIFT" = "sift",
                                          "PolyPhen-2" = "polyphen2"
                                        ),
                                        selected = c("clinvar", "gnomad")),
                     # Option to run OpenCRAVAT and fetch the scores
                     actionButton("fetchButton", "Fetch VEP Data")
                   ),
                   
                   # selectizeInput("gene", "Select Gene Name:", choices = NULL, options = list(maxOptions = 1000)),
                   
                   # Additional controls for the plot
                   # checkboxInput("common_variant_filter", "Exclude Common Variants (gnomAD AF > 0.005)", value = TRUE),
                   checkboxGroupInput("scores", "Select Scores to Include:", choices = NULL, selected = NULL),
                   actionButton("plotButton", "Generate PRC Plot"),
                   downloadButton("downloadPlotPNG", "Download PRC Plot as PNG"),
                   downloadButton("downloadPlotPDF", "Download PRC Plot and Metadata as PDF"),
                 ),
                 
                 mainPanel(
                   plotOutput("prcPlot", width = "600px", height = "600px"),
                   textOutput("errorText")
                 )
               )
             )
             
    )
  )
)
