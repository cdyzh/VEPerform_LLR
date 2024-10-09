# Load necessary libraries
library(plumber)
library(shiny)
library(dplyr)
library(yogiroc)

#* @apiTitle VEPerform API
#* @apiDescription This API provides access to VEPerform's functionalities.

#* Health check endpoint to verify the API is running
#* @get /health
function() {
  list(status = "API is running")
}

#* Generate a PRC plot based on provided parameters
#* @param gene The gene name
#* @param scores A comma-separated list of scores (e.g., "VARITY,REVEL,AlphaMissense")
#* @param common_variant_filter Should common variants be excluded? (TRUE/FALSE)
#* @post /generate_prc
#* @serializer contentType list(type="image/png")
function(gene, scores = "VARITY,REVEL,AlphaMissense", common_variant_filter = TRUE) {
  scores_list <- unlist(strsplit(scores, ","))
  common_variant_filter <- as.logical(common_variant_filter)
  
  # Load the data
  df <- read.csv("preprocessed_id.csv", stringsAsFactors = FALSE)
  
  # Standardize column names to match Shiny app
  colnames(df) <- c(
    "base__hugo", 
    "base__achange",
    "gnomAD_AF", 
    "VARITY", 
    "AlphaMissense", 
    "REVEL", 
    "classification"
  )
  
  # Apply common variant filter if required
  if (common_variant_filter) {
    df <- df[is.na(df$gnomAD_AF) | df$gnomAD_AF <= 0.005, ]
  }
  
  # Filter data for the selected gene
  prcfiltered <- df %>% filter(base__hugo == gene)
  
  # Generate the PRC plot
  B_org <- sum(prcfiltered$classification == TRUE & rowSums(!is.na(prcfiltered[scores_list])) > 0)
  P_org <- sum(prcfiltered$classification == FALSE & rowSums(!is.na(prcfiltered[scores_list])) > 0)
  
  tryCatch({
    yrobj <- yr2(truth = prcfiltered[['classification']], scores = prcfiltered[scores_list], high = rep(FALSE, length(scores_list)))
    
    # Create the plot and save it to a temporary file
    tmp <- tempfile(fileext = ".png")
    png(tmp, width = 600, height = 600)
    draw.prc(yrobj, lty = c("dashed", "solid", "dashed")[1:length(scores_list)], col = c("purple", "cadetblue2", "orange")[1:length(scores_list)], lwd = 2, balanced = TRUE, main = paste0(gene, " PRCs for ", paste(scores_list, collapse = ", ")))
    abline(h = 90, lty = "dashed")
    legend("left", legend = c(paste("# of Pathogenic and Likely Pathogenic:", P_org), paste("# of Benign and Likely Benign:", B_org)), pch = 15, bty = "n")
    dev.off()
    
    # Return the image as a response
    readBin(tmp, "raw", n = file.info(tmp)$size)
  }, error = function(e) {
    stop("Not enough data to generate PRC plot.")
  })
}

#* Generate a PDF report based on provided parameters
#* @param gene The gene name
#* @param scores A comma-separated list of scores (e.g., "VARITY,REVEL,AlphaMissense")
#* @param common_variant_filter Should common variants be excluded? (TRUE/FALSE)
#* @post /generate_report
#* @serializer contentType list(type="application/pdf")
function(gene, scores = "VARITY,REVEL,AlphaMissense", common_variant_filter = TRUE) {
  scores_list <- unlist(strsplit(scores, ","))
  common_variant_filter <- as.logical(common_variant_filter)
  
  # Load the data
  df <- read.csv("preprocessed_id.csv", stringsAsFactors = FALSE)
  
  # Clean column names to avoid NA or empty names
  df <- df[, !is.na(colnames(df)) & colnames(df) != ""]
  
  # Standardize column names to match Shiny app
  colnames(df) <- c(
    "base__hugo", 
    "base__achange",
    "gnomAD_AF", 
    "VARITY", 
    "AlphaMissense", 
    "REVEL", 
    "classification"
  )
  
  # Apply common variant filter if required
  if (common_variant_filter) {
    df <- df[is.na(df$gnomAD_AF) | df$gnomAD_AF <= 0.005, ]
  }
  
  # Filter data for the selected gene
  prcfiltered <- df %>% filter(base__hugo == gene)
  
  # Generate the report using RMarkdown
  tmp <- tempfile(fileext = ".pdf")
  tryCatch({
    rmarkdown::render(input = "report_template.Rmd",
                      output_file = tmp,
                      params = list(
                        gene_s = gene,
                        selected_scores = scores_list,
                        B_org = sum(prcfiltered$classification == TRUE & rowSums(!is.na(prcfiltered[scores_list])) > 0),
                        P_org = sum(prcfiltered$classification == FALSE & rowSums(!is.na(prcfiltered[scores_list])) > 0),
                        prcfiltered = prcfiltered
                      ),
                      envir = new.env(parent = globalenv()))
    
    # Return the PDF as a response
    readBin(tmp, "raw", n = file.info(tmp)$size)
  }, error = function(e) {
    stop("Not enough data to generate PDF report.")
  })
}
