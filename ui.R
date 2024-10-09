# ui.R
library(shiny)
library(shinyjs)
library(dplyr)
library(yogiroc)
library(rsconnect)
library(rmarkdown)
library(knitr)
library(tinytex)

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
                      p("This tool allows you to evaluate the performance of variant effect predictors."),
                      p("You can either use our pre-existing dataset or upload your own dataset to evaluate a subset of variants."),
                      
                      # Radio buttons for selecting between using an existing dataset or uploading a new one
                      radioButtons("intro_data_source", "Please select an option:",
                                   choices = list(
                                     "Use Pre-existing Dataset" = "existing", 
                                     "Upload Your Own Dataset" = "upload"
                                   ),
                                   selected = "existing"),
                      
                      # Link to API
                      p("You can also interact with the VEPerform API ", 
                        a("here", href = "http://localhost:8000", target = "_blank")),
                      
                      # Action button to proceed
                      actionButton("proceed_button", "Proceed", class = "btn-success")
               )
             )
    ),
    
    # Main App (Use Pre-existing Dataset) Tab
    tabPanel("Main App",
             sidebarLayout(
               sidebarPanel(
                 # Radio buttons to select between uploading a dataset or using an existing one
                 radioButtons("data_source", "Choose Dataset:",
                              choices = list("Use Existing Full Dataset" = "existing", 
                                             "Upload Your Own Dataset to Select Subset of Variants" = "upload"),
                              selected = "existing"),
                 
                 # Conditional panel that shows upload options if the user chooses to upload their own dataset
                 conditionalPanel(
                   condition = "input.data_source == 'upload'",
                   fileInput("file_gene_variant", "Upload Gene and HGVS_Pro CSV File", 
                             accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
                   actionButton("upload_guide_existing", "Upload Guide", class = "btn-info")
                 ),
                 
                 selectizeInput("gene", "Select Gene Name:", choices = NULL, options = list(maxOptions = 1000)),
                 checkboxInput("common_variant_filter", "Exclude Common Variants (gnomAD AF > 0.005)", value = TRUE),
                 checkboxGroupInput("scores", "Select Scores to Include:",
                                    choices = list("VARITY", "REVEL", "AlphaMissense"),
                                    selected = c("VARITY", "REVEL", "AlphaMissense")),
                 actionButton("plotButton", "Generate PRC Plot"),
                 downloadButton("downloadPlotPNG", "Download PRC Plot as PNG"),
                 downloadButton("downloadPlotPDF", "Download PRC Plot and Metadata")
               ),
               mainPanel(
                 plotOutput("prcPlot", width = "600px", height = "600px"),
                 textOutput("errorText")
               )
             )
    ),
    
    # Upload Your Own Dataset Tab
    tabPanel("Upload",
             sidebarLayout(
               sidebarPanel(
                 fileInput("file_full", "Upload Full Dataset CSV File", 
                           accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
                 actionButton("upload_guide", "Upload Guide", class = "btn-info"),
                 
                 selectizeInput("gene", "Select Gene Name:", choices = NULL, options = list(maxOptions = 1000)),
                 checkboxInput("common_variant_filter", "Exclude Common Variants (gnomAD AF > 0.005)", value = TRUE),
                 checkboxGroupInput("scores", "Select Scores to Include:",
                                    choices = list("VARITY", "REVEL", "AlphaMissense"),
                                    selected = c("VARITY", "REVEL", "AlphaMissense")),
                 actionButton("plotButton", "Generate PRC Plot"),
                 downloadButton("downloadPlotPNG", "Download PRC Plot as PNG"),
                 downloadButton("downloadPlotPDF", "Download PRC Plot and Metadata")
               ),
               mainPanel(
                 plotOutput("prcPlot", width = "600px", height = "600px"),
                 textOutput("errorText")
               )
             )
    )
  )
)
