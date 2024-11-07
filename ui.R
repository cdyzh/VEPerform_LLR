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
    title = HTML("<strong style='color:blue;'>VEPerform</strong>"),
    id = "navbar",  # Set an ID for the navbarPage
    
    # Main App (Use Pre-existing Dataset) Tab
    tabPanel("VEPerform",
             fluidPage(
               # Introductory Text
               tags$h3("Welcome to VEPerform"),
               tags$p(HTML("
                      <p>How well do variant effect predictors (VEPs) work for your favorite gene? </p>
                      <p>Of course, performance depends on what VEP score threshold is used.</p>
                      <p>VEPerform helps you visualize this by plotting:</p>
                      <ul>
                        <li><strong>Precision</strong> (% of pathogenicity predictions that are correct) versus</li>
                        <li><strong>Recall</strong> (% of known pathogenic variants that were predicted)</li>
                      </ul>
                      <p>Select a gene to begin.</p>
                    ")),
               
               sidebarLayout(
                 sidebarPanel(
                   selectizeInput("Main_gene", "Gene of Interest:", choices = NULL, options = list(maxOptions = 1000)),
                   # Radio buttons to select between uploading a dataset or using an existing one
                   radioButtons("data_source", "Reference Set of Pathogenic/Benign:",
                                choices = list("Default (ClinVar)" = "existing", 
                                               "Your Custom Subset" = "upload"),
                                selected = "existing"),
                   
                   helpText(HTML("
                    <p><strong>NOTE:</strong> Alternatively, upload your own reference set via the <strong style='color:blue;'>Advanced</strong> VEPerform tab.</p>")),
                   
                   # Conditional panel that shows upload options if the user chooses to upload their own dataset
                   conditionalPanel(
                     condition = "input.data_source == 'upload'",
                     helpText("Using the outputted CSV, delete the rows for variants you would like to remove and keep only the gene and hgvs_pro columns." ),
                     fileInput("file_gene_variant", "Upload Gene and HGVS_Pro CSV File", 
                               accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv"))
                     # actionButton("upload_guide_existing", "Upload Guide", class = "btn-info")
                   ),
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
                     downloadButton("download_template_oc", "Download CSV Template"),
                     helpText("Click 'Fetch VEP Data' after uploading your file to fetch data from OpenCRAVAT"),
                     
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
             
    ),
    # About Tab
    tabPanel("About",
             fluidRow(
               column(12, offset = 0,
                      h2("About VEPerform"),
                      p("This tool allows you to evaluate the performance of variant effect predictors for your favorite genes."),
                      p("You can either use our pre-existing dataset, fetch data from OpenCRAVAT, or upload your own dataset."),
                      
                      # Link to API
                      p("You can also interact with the VEPerform API ", 
                        a("here", href = "http://127.0.0.1:4697/", target = "_blank"))
               )
             )
    )
  )
)
