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
                                     tags$span(tags$b("Fetch variant data live from OpenCRAVAT:"), tags$br(), "for the most up-to-date scores and annotations"),
                                     tags$span(tags$b("Customize and upload your own reference set:"), tags$br(), "if you have your own data and just want to use the PRC generating tool")
                                   ),
                                   choiceValues = list("existing", "fetch", "upload"),
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
    
    tabPanel("Fetch Live", 
             fluidPage(
               # Introductory Text
               tags$h3("Fetch Data from ClinVar, GnomAD, and OpenCRAVAT"),
               tags$p("This option offers maximal input customizability but takes longer to run. In this section, 
                    you can dynamically create a dataset with custom ClinVar filters and additional 
                    proxy-benign variants from GnomAD. This is then fed into OpenCRAVAT for variant predictor scores."),
               
             sidebarLayout(
               sidebarPanel(
                 radioButtons("input_type", "Select Input Type:",
                              choices = list("Chromosome, Position, Reference_Base, Alternate_Base" = "chrom_pos",
                                             "Transcript ID and HGVSC (in development)" = "hgvsc")),
                 
                 # File input for Chromosome, Position, Reference_Base, Alternate_Base
                 conditionalPanel(
                   condition = "input.input_type == 'chrom_pos'",
                   fileInput("variant_file", "Upload CSV File (Chromosome, Position, Reference_Base, Alternate_Base)",
                             accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv"))
                 ),
                 
                 # File input for Transcript ID and HGVSC
                 conditionalPanel(
                   condition = "input.input_type == 'hgvsc'",
                   fileInput("variant_file_hgvsc", "Upload CSV File (Transcript ID, HGVSC)",
                             accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv"))
                 ),
                 
                 helpText("Click 'Fetch VEP Data' after uploading your file to fetch data from OpenCRAVAT" ),
                 
                 # Option to run OpenCRAVAT and fetch the scores
                 actionButton("fetchButton", "Fetch VEP Data"),
                 
                 # Additional controls for the plot
                 checkboxInput("Fetch_common_variant_filter", "Exclude Common Variants (gnomAD AF > 0.005)", value = TRUE),
                 checkboxGroupInput("Fetch_scores", "Select Scores to Include:",
                                    choices = list(
                                      "VARITY",
                                      "REVEL",
                                      "AlphaMissense"
                                    ),
                                    selected = c("VARITY", "REVEL", "AlphaMissense")
                 ),
                 actionButton("Fetch_plotButton", "Generate PRC Plot"),
                 downloadButton("Fetch_downloadPlotPNG", "Download PRC Plot as PNG"),
                 downloadButton("Fetch_downloadPlotPDF", "Download PRC Plot and Metadata as PDF"),
             ),
               
               mainPanel(
                 plotOutput("Fetch_prcPlot", width = "600px", height = "600px"),
                 textOutput("Fetch_errorText")
               )
             )
             )
             
          ),
    
    # Upload Your Own Dataset Tab
    tabPanel("Custom",
             fluidPage(
               # Introductory Text
               tags$h3("Upload Your Own Dataset"),
               tags$p("In this section, you can upload your own custom reference set. 
                   Ensure your file format aligns with the required specifications. Once uploaded, you can select your desired gene and predictor scores for analysis. 
                   After configuring these settings, generate a Precision-Recall Curve (PRC) to assess variant effect predictor performance."),
             sidebarLayout(
               sidebarPanel(
                 fileInput("file_full", "Upload Full Reference Set CSV File", 
                           accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
                 actionButton("upload_guide", "Upload Guide", class = "btn-info"),
                 
                 selectizeInput("gene", "Select Gene Name:", choices = NULL, options = list(maxOptions = 1000)),
                 # Common variant filter is moved to server logic
                 # checkboxInput("common_variant_filter", "Exclude Common Variants (gnomAD AF > 0.005)", value = TRUE), 
                 # Set choices to NULL initially, to be updated dynamically
                 checkboxGroupInput("scores", "Select Scores to Include:", choices = NULL, selected = NULL),
                 actionButton("plotButton", "Generate PRC Plot"),
                 downloadButton("downloadPlotPNG", "Download PRC Plot as PNG"),
                 downloadButton("downloadPlotPDF", "Download PRC Plot and Metadata as PDF")
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
