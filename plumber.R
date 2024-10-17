# Load necessary libraries
library(plumber)
library(shiny)
library(dplyr)
library(yogiroc)
options("plumber.port" = 8000)

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
  
  # Filter data for the selected gene, ignoring case (case-insensitive)
  gene <- tolower(gene)
  prcfiltered <- df %>% filter(tolower(base__hugo) == gene)
  
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

#* Output the reference set used for PRC plot generation in JSON format
#* @param gene The gene name
#* @param scores A comma-separated list of scores (e.g., "VARITY,REVEL,AlphaMissense")
#* @param common_variant_filter Should common variants be excluded? (TRUE/FALSE)
#* @get /get_referenceset
#* @serializer json
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
  
  # Filter data for the selected gene, ignoring case (case-insensitive)
  gene <- tolower(gene)
  prcfiltered <- df %>% filter(tolower(base__hugo) == gene)
  
  # Return the filtered dataframe as JSON
  prcfiltered
}
