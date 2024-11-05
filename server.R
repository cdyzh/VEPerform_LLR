# server.R
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

# Define server logic for the application
server <- function(input, output, session) {
  plot_data <- reactiveVal(NULL)
  selected_variants <- reactiveVal(NULL)  # Store user-selected variants
  
  # Observe the Proceed button to switch tabs
  observeEvent(input$proceed_button, {
    if (input$intro_data_source == "fetch") {
      updateNavbarPage(session, "navbar", selected = "Advanced")
    } else {
      updateNavbarPage(session, "navbar", selected = "Main App")
    }
  })
  
  # Giant if/else block to handle separate logic for Main App and Upload tabs
  observe({
    if (input$navbar == "Main App") {
      # prcdata logic for Main App
      # prcdata logic for Main App
      prcdata <- reactive({
        if (input$data_source == "upload" && !is.null(input$file_gene_variant)) {
          req(input$file_gene_variant)
          df <- read.csv(input$file_gene_variant$datapath, stringsAsFactors = FALSE)
          
          if (ncol(df) != 2) {
            showModal(modalDialog(
              title = "Error",
              "The uploaded dataset must have exactly two columns: base__gene and base__achange.",
              easyClose = TRUE,
              footer = modalButton("Close")
            ))
            reset("file_gene_variant")
            updateFileInput(session, "file_gene_variant", value = NULL)
            return(NULL)
          }
          
          colnames(df) <- c("base__gene", "base__achange")
          full_df <- read.csv("preprocessed.csv", stringsAsFactors = FALSE)
          df <- merge(df, full_df, by = c("base__gene", "base__achange"))
        } else {
          df <- read.table("preprocessed.csv", sep = ',', header = TRUE, stringsAsFactors = FALSE)
        }
        return(df)
      })
      
      # Update gene names
      observe({
        df <- prcdata()
        if (!is.null(df)) {
          gene_names <- unique(df$base__gene)
          updateSelectizeInput(session, "Main_gene", choices = gene_names, selected = gene_names[1], server = TRUE)
        }
      })
      
      # Update scores
      observe({
        req(input$Main_gene)
        df <- prcdata()
        gene_data <- df %>% filter(base__gene == input$Main_gene)
        
        available_scores <- c()
        if (any(!is.na(gene_data$varity_r))) {
          available_scores <- c(available_scores, "VARITY")
        }
        if (any(!is.na(gene_data$alphamissense__pathogenicity))) {
          available_scores <- c(available_scores, "AlphaMissense")
        }
        if (any(!is.na(gene_data$revel__score))) {
          available_scores <- c(available_scores, "REVEL")
        }
        
        updateCheckboxGroupInput(session, "Main_scores", choices = available_scores, selected = available_scores)
      })
      
      # Function to show modal with variant selection
      showVariantSelectionModal <- function(df, gene_selected) {
        
        # Create a modified dataframe for display
        gene_display_df <- df %>% filter(base__gene == gene_selected) %>%
          select(base__gene, base__achange, clinvar, revstat, stars)
        
        
        showModal(modalDialog(
          title = paste("Select Variants for Gene:", gene_selected),
          tags$p("Would you like to de-select some variants? Click on the variants you do not want to include in the PRC. You can also do this later by downloading the csv, removing variants yourself, and uploading it back."),
          DTOutput("variant_table"),
          footer = tagList(
            modalButton("Cancel"),
            actionButton("confirm_selection", "OK")
          ),
          easyClose = FALSE
        ))
        
        # Render the interactive table with checkboxes for the filtered dataframe
        output$variant_table <- renderDT({
          datatable(
            gene_display_df, 
            selection = list(target = "row", selected = 1:nrow(gene_display_df)),  # All rows selected by default
            options = list(pageLength = 10, scrollX = TRUE)
          )
        }, server = TRUE)
      }
      
      # Generate PRC and open modal for variant selection
      observeEvent(input$Main_plotButton, {
        gene_s <- input$Main_gene  # Get the selected gene
        prcfiltered <- prcdata() %>%
          filter(base__gene == gene_s)
        if (!is.null(df)) {
          showVariantSelectionModal(prcfiltered, gene_s)  # Show modal with filtered variants for the selected gene
        }
      })
      
      # Handle the confirmation from the modal
      observeEvent(input$confirm_selection, {
        selected_rows <- input$variant_table_rows_selected  # Get selected rows from the table
        # df <- prcdata()
        gene_s <- input$Main_gene  # Get the selected gene
        prcfiltered <- prcdata() %>%
          filter(base__gene == gene_s)
        
        if (!is.null(prcfiltered) && length(selected_rows) > 0) {
          # Update the selected variants based on the user's choice
          print(selected_rows) #DEBUG
          selected_df <- prcfiltered[selected_rows, ]
          selected_variants(selected_df)
          
          # Proceed with the PRC generation using the selected variants
          gene_s <- input$Main_gene
          print(gene_s) # DEBUG
          
          exclude_common_variants <- input$Main_common_variant_filter
          selected_scores <- input$Main_scores
          print("begin plot with selected variants") # DEBUG
          
          names(selected_df)[names(selected_df) == "varity_r"] <- "VARITY"
          names(selected_df)[names(selected_df) == "alphamissense__pathogenicity"] <- "AlphaMissense"
          names(selected_df)[names(selected_df) == "revel__score"] <- "REVEL"
          names(selected_df)[names(selected_df) == "gnomad__af"] <- "gnomAD_AF"
          
          if (exclude_common_variants) {
            selected_df <- selected_df[is.na(selected_df$gnomAD_AF) | selected_df$gnomAD_AF <= 0.005, ]
          }
          
          selected_df <- selected_df[order(selected_df$clinvar), ]
          
          selected_df <- selected_df %>% mutate(clinvar = ifelse(clinvar == "B/LB", TRUE, FALSE))
          print(selected_df) # DEBUG
          
          B_org <- sum(selected_df$clinvar == TRUE & rowSums(!is.na(selected_df[selected_scores])) > 0)
          P_org <- sum(selected_df$clinvar == FALSE & rowSums(!is.na(selected_df[selected_scores])) > 0)
          
          tryCatch({
            yrobj <- yr2(truth = selected_df[["clinvar"]], scores = selected_df[selected_scores], high = rep(FALSE, length(selected_scores)))
            
            plot_data(list(
              yrobj = yrobj,
              lty_styles = c("dashed", "solid", "dashed")[1:length(selected_scores)],
              col_styles = c("purple", "cadetblue2", "orange")[1:length(selected_scores)],
              gene_s = gene_s,
              selected_scores = selected_scores,
              B_org = B_org,
              P_org = P_org,
              prcfiltered = selected_df  # Save the filtered data for metadata
            ))
            print("plot success with selected variants")  # DEBUG
            output$Main_ErrorText <- renderText("")
            
          }, error = function(e) {
            plot_data(NULL)
            output$Main_ErrorText <- renderText("Not enough data")
          })
          
          output$Main_PRCPlot <- renderPlot({
            plot_info <- plot_data()
            if (!is.null(plot_info)) {
              tryCatch({
                draw.prc(plot_info$yrobj, lty = plot_info$lty_styles, col = plot_info$col_styles, lwd = 2, balanced = TRUE, main = paste0(plot_info$gene_s, " PRCs for ", paste(plot_info$selected_scores, collapse = ", ")))
                abline(h = 90, lty = "dashed")
                legend("left", legend = c(paste("# of Pathogenic and Likely Pathogenic:", plot_info$P_org), paste("# of Benign and Likely Benign:", plot_info$B_org)), pch = 15, bty = "n")
              }, error = function(e) {
                showModal(modalDialog(
                  title = 'Error',
                  'Not enough data',
                  easyClose = TRUE,
                  footer = NULL
                ))
              })
            }
          }, width = 600, height = 600, res = 72)
          
          removeModal()  # Close the modal after the user clicks "OK"
        } else {
          showModal(modalDialog(
            title = "Error",
            "Please select at least one variant.",
            easyClose = TRUE,
            footer = modalButton("Close")
          ))
        }
      })
      
      # Download logic for Main App
      output$Main_downloadPlotPNG <- downloadHandler(
        filename = function() {
          paste("PRC_plot_", input$Main_gene, ".png", sep = "")
        },
        content = function(file) {
          plot_info <- plot_data()
          if (!is.null(plot_info)) {
            png(file, width = 6, height = 6, units = "in", res = 72)
            draw.prc(plot_info$yrobj, lty = plot_info$lty_styles, col = plot_info$col_styles, lwd = 2, balanced = TRUE, main = paste0(plot_info$gene_s, " PRCs for ", paste(plot_info$selected_scores, collapse = ", ")))
            abline(h = 90, lty = "dashed")
            legend("left", legend = c(paste("# of Pathogenic and Likely Pathogenic:", plot_info$P_org), paste("# of Benign and Likely Benign:", plot_info$B_org)), pch = 15, bty = "n")
            dev.off()
          }
        }
      )
      
      output$Main_downloadPlotPDF <- downloadHandler(
        filename = function() {
          paste("PRC_Report_", input$Main_gene, ".pdf", sep = "")
        },
        content = function(file) {
          plot_info <- plot_data()
          if (!is.null(plot_info)) {
            # Generate the PDF report
            rmarkdown::render(input = "report_template.Rmd",
                              output_file = file,
                              params = list(
                                gene_s = plot_info$gene_s,
                                selected_scores = plot_info$selected_scores,
                                B_org = plot_info$B_org,
                                P_org = plot_info$P_org,
                                prcfiltered = plot_info$prcfiltered
                              ),
                              envir = new.env(parent = globalenv()))
          }
        }
      )
      # Download prcfiltered as CSV
      output$downloadCSV <- downloadHandler(
        filename = function() {
          paste("PRC_data_", input$Main_gene, ".csv", sep = "")
        },
        content = function(file) {
          df <- prcdata()  # Get the reactive data
          gene_s <- input$Main_gene
          
          # Filter data for the selected gene
          prcfiltered <- df %>%
            filter(base__gene == gene_s)
          
          # Write the prcfiltered dataframe to CSV
          write.csv(prcfiltered, file, row.names = FALSE)
        }
      )
        
        
    } else {
      
      #logic for advanced - own
        if (input$input_type == "own") {
          prcdata <- reactive({
            req(input$file_full)
            df <- tryCatch({
              read.csv(input$file_full$datapath, stringsAsFactors = FALSE)
            }, error = function(e) {
              showModal(modalDialog(
                title = "Error",
                "There was an error reading the uploaded file. Please ensure it is a valid CSV file.",
                easyClose = TRUE,
                footer = modalButton("Close")
              ))
              return(NULL)
            })
            
            print(colnames(df)) # DEBUG
            
            # Mandatory columns
            colnames(df)[colnames(df) == "gnomad__af"] <- "gnomAD_AF"
            colnames(df)[colnames(df) == "clinvar"] <- "clinvar"
            
            return(df)
          })
          
          # Guide on how to format user-inputted csv
          observeEvent(input$upload_guide, {
            showModal(modalDialog(
              title = "Reference Set Format Information",
              HTML("Please ensure your own reference set is a CSV file.<br><br>
            Mandatory columns:<br>
            <b>gene:</b> Gene name(s)<br>
            <b>clinvar:</b> Variant clinvar (TRUE for Benign, FALSE for Pathogenic)<br><br>
            
            Optional columns: <br>
            <b>gnomad_af:</b> GnomAD allele frequency - for filtering out common variants<br><br>
            
            For any predictors or custom scores that you would like to include, put 'VEP_' before the name of the column. For example:<br>
            <b>VEP_alphamissense:</b> AlphaMissense score<br>
            <b>VEP_custom:</b> Additional scores you would like to include<br><br>"),
              easyClose = TRUE,
              footer = NULL
            ))
          })
          
        } else if (input$input_type == "fetch")  {
          
          # logic for fetch
          prcdata <- reactiveVal(NULL)
          plot_data <- reactiveVal(NULL)
          
          observeEvent(input$fetchButton, {
           # req(input$input_type)
            
              # Initialize empty list to store results
            all_results <- list()
              
            req(input$file_fetch)
            df <- read.csv(input$file_fetch$datapath, stringsAsFactors = FALSE)
            colnames(df) <- trimws(colnames(df))
            # Ensure 'chrom' column in df has the "chr" prefix for consistency with result_df
            df <- df %>%
              mutate(chrom = ifelse(grepl("^chr", as.character(chrom)), as.character(chrom), paste0("chr", chrom)))
            
            # Check for necessary columns
            if (!all(c("chrom", "pos", "ref_base", "alt_base") %in% colnames(df))) {
              output$errorText <- renderText("Error: The uploaded CSV must contain the columns 'chrom', 'pos', 'ref_base', 'alt_base'.")
              return(NULL)
            }
              
            # Show a progress bar while fetching data
            withProgress(message = 'Fetching Variant Data', value = 0, {
              for (i in 1:nrow(df)) {
                
                print(colnames(df)) # DEBUG
                print(df) # DEBUG
                print(nrow(df))
                
                chrom <- df$chrom[i]
                pos <- df$pos[i]
                ref <- df$ref_base[i]
                alt <- df$alt_base[i]
                
                # Construct API URL for chrom-based input
                api_url <- paste0(
                  "https://run.opencravat.org/submit/annotate?",
                  "chrom=", chrom,
                  "&pos=", pos,
                  "&ref_base=", ref,
                  "&alt_base=", alt,
                  "&annotators=clinvar,gnomad,varity_r,revel,alphamissense"
                )
                
                # Make the GET request to OpenCRAVAT
                response <- GET(api_url)
                
                # Parse the JSON response
                result <- fromJSON(content(response, "text"), flatten = TRUE)
                
                #print(result) # DEBUG
                
              
              
                # Handle cases where certain annotations might be missing
                clinvar_sig <- ifelse(!is.null(result$clinvar$sig), result$clinvar$sig, NA)
                
                # Determine the clinvar value (T/F) based on ClinVar significance
                clinvar <- NA
                if (!is.na(clinvar_sig)) {
                  clinvar_sig_lower <- tolower(clinvar_sig)
                  
                  if (grepl("benign", clinvar_sig_lower)) {
                    clinvar <- "B/LB"
                  } else if (grepl("pathogenic", clinvar_sig_lower) && !grepl("conflicting", clinvar_sig_lower)) {
                    clinvar <- "P/LP"
                  }
                }
                
                # Extract into corresponding columns - TODO: be able to add additional predictors
                gene <- ifelse(!is.null(result$crx$hugo), result$crx$hugo, NA)
                achange <- ifelse(!is.null(result$crx$achange), result$crx$achange, NA)
                gnomad_af <- ifelse(!is.null(result$gnomad$af), result$gnomad$af, NA)
                varity_r <- ifelse(!is.null(result$varity_r$varity_r), result$varity_r$varity_r, NA)
                revel_score <- ifelse(!is.null(result$revel$score), result$revel$score, NA)
                alphamissense_path <- ifelse(!is.null(result$alphamissense$am_pathogenicity), result$alphamissense$am_pathogenicity, NA)
                
                # Create a data frame with consistent column names
                result_df <- data.frame(
                  chrom = chrom,
                  pos = pos,
                  ref_base = ref,
                  alt_base = alt,
                  base__gene = gene,
                  base__achange = achange,
                  gnomad__af = gnomad_af,
                  varity_r = varity_r,
                  revel__score = revel_score,
                  alphamissense__pathogenicity = alphamissense_path,
                  clinvar = clinvar,
                  stringsAsFactors = FALSE
                )
                
                # Append to the results list
                all_results[[i]] <- result_df
                
                incProgress(1 / nrow(df))
                showNotification(paste(i, "/", nrow(df), "variants fetched"), duration = 3, type = "message")
                print("here") # DEBUG
                }
                print("done1")
              }
              )
              print("done2")
              
              # Combine all results into a data frame
              variant_data_df <- do.call(rbind, all_results)
              
              # Join variant_data_df with the original df on chrom, pos, ref_base, alt_base
              variant_data_df <- dplyr::left_join(df, variant_data_df, by = c("chrom", "pos", "ref_base", "alt_base"))
              
              # Update reactive variable
              prcdata(variant_data_df)
              
              # Clear any error message
              output$errorText <- renderText("")
              print("fetch complete") # DEBUG
        })
      }
      
      df <- prcdata()
      
      # Initialize reactiveValues to store the state of the gnomad checkbox insertion
      state <- reactiveValues(gnomad_filter_inserted = FALSE)
      
      # Observe changes to the processed data only when prcdata changes
      observe({
        req(input$file_full)
        
        # Use isolate to prevent re-triggering prcdata
        df <- isolate(prcdata())
        if (is.null(df)) return()  # Exit if prcdata is NULL
        
        # Check for gnomAD columns, case-insensitive
        colnames_lower <- tolower(colnames(df))
        gnomad_columns <- grep("gnomad", colnames_lower, value = TRUE)
        
        # If gnomAD columns are found and the checkbox is not yet added
        if (length(gnomad_columns) > 0 && !state$gnomad_filter_inserted) {
          # Remove any existing checkbox to prevent duplicates
          removeUI(selector = "#gnomad_filter_wrapper", immediate = TRUE)
          
          # Insert the checkbox
          insertUI(
            selector = "#upload_guide",
            where = "afterEnd",
            ui = div(id = "gnomad_filter_wrapper",
                     checkboxInput("common_variant_filter", 
                                   "Exclude Common Variants (gnomAD AF > 0.005)", 
                                   value = TRUE)
            )
          )
          state$gnomad_filter_inserted <- TRUE
          
        } else if (length(gnomad_columns) == 0 && state$gnomad_filter_inserted) {
          # Remove the checkbox if no gnomAD columns are found and it’s currently inserted
          removeUI(selector = "#gnomad_filter_wrapper", immediate = TRUE)
          state$gnomad_filter_inserted <- FALSE
        }
      })
      
      
      
      observe({
        req(input$file_full)
        df <- prcdata()  # Assuming prcdata() reads the uploaded CSV
        
        # Detect all columns with the "VEP_" prefix
        predictor_columns <- colnames(df)[grepl("^VEP_", colnames(df))]
        
        # Remove "VEP_" prefix from predictor column names
        predictor_columns <- gsub("^VEP_", "", predictor_columns)
        
        # Update checkboxGroupInput with detected predictor columns
        updateCheckboxGroupInput(session, "scores", choices = predictor_columns, selected = predictor_columns)
      })
      
      observeEvent(input$plotButton, {
        df <- prcdata()
        df <- df[!is.na(df$clinvar), ]
        print(df) # DEBUG
        
        if (is.null(df) || nrow(df) == 0) {
          output$errorText <- renderText("Not enough rows to generate the PRC plot.")
          return()
        }
        
        gene_s <- df$base__gene[1]
        exclude_common_variants <- input$common_variant_filter
        selected_scores <- input$scores
        
        # Rename columns to match what you need for PRC plotting
        names(df)[names(df) == "varity_r"] <- "VARITY"
        names(df)[names(df) == "alphamissense__pathogenicity"] <- "AlphaMissense"
        names(df)[names(df) == "revel__score"] <- "REVEL"
        names(df)[names(df) == "gnomad__af"] <- "gnomAD_AF"
        
        if (exclude_common_variants) {
          df <- df[is.na(df$gnomAD_AF) | df$gnomAD_AF <= 0.005, ]
        }
        
        # Remove "P-" prefix from predictor column names
        colnames(df) <- gsub("^VEP_", "", colnames(df))
        
        df <- df[order(df$clinvar), ]
        
        print(df) # DEBUG
        
        # Filter for the selected scores
        prcfiltered <- df %>%
          filter(rowSums(!is.na(df[selected_scores])) > 0) %>%
          mutate(clinvar = ifelse(clinvar == "B/LB", TRUE, FALSE))
        
        B_org <- sum(prcfiltered$clinvar == TRUE & rowSums(!is.na(prcfiltered[selected_scores])) > 0)
        P_org <- sum(prcfiltered$clinvar == FALSE & rowSums(!is.na(prcfiltered[selected_scores])) > 0)
        
        tryCatch({
          yrobj <- yr2(truth = prcfiltered[["clinvar"]], scores = prcfiltered[selected_scores], high = rep(FALSE, length(selected_scores)))
          
          # Store plot information in plot_info
          plot_data(list(
            yrobj = yrobj,
            lty_styles = c("dashed", "solid", "dashed")[1:length(selected_scores)],
            col_styles = c("purple", "cadetblue2", "orange")[1:length(selected_scores)],
            gene_s = gene_s,
            selected_scores = selected_scores,
            B_org = B_org,
            P_org = P_org,
            prcfiltered = prcfiltered  # Save the filtered data for metadata
          ))
          
          output$errorText <- renderText("")
          
        }, error = function(e) {
          plot_data(NULL)
          output$errorText <- renderText("Not enough data to generate PRC plot.")
        })
        
        output$prcPlot <- renderPlot({
          plot_info <- plot_data()
          if (!is.null(plot_info)) {
            tryCatch({
              draw.prc(plot_info$yrobj, lty = plot_info$lty_styles, col = plot_info$col_styles, lwd = 2, balanced = TRUE, main = paste0(plot_info$gene_s, " PRCs for ", paste(plot_info$selected_scores, collapse = ", ")))
              abline(h = 90, lty = "dashed")
              legend("left", legend = c(paste("# of Pathogenic and Likely Pathogenic:", plot_info$P_org), paste("# of Benign and Likely Benign:", plot_info$B_org)), pch = 15, bty = "n")
            }, error = function(e) {
              showModal(modalDialog(
                title = 'Error',
                'Not enough data',
                easyClose = TRUE,
                footer = NULL
              ))
            })
          }
        }, width = 600, height = 600, res = 72)
      })
      
      output$downloadPlotPNG <- downloadHandler(
        filename = function() {
          paste("PRC_plot_", Sys.Date(), ".png", sep = "")
        },
        content = function(file) {
          plot_info <- plot_data()
          if (!is.null(plot_info)) {
            png(file, width = 6, height = 6, units = "in", res = 72)
            draw.prc(plot_info$yrobj, lty = plot_info$lty_styles, col = plot_info$col_styles, lwd = 2, balanced = TRUE, main = paste0(plot_info$gene_s, " PRCs for ", paste(plot_info$selected_scores, collapse = ", ")))
            abline(h = 90, lty = "dashed")
            legend("left", legend = c(paste("# of Pathogenic and Likely Pathogenic:", plot_info$P_org), paste("# of Benign and Likely Benign:", plot_info$B_org)), pch = 15, bty = "n")
            dev.off()
          }
        }
      )
      
      output$downloadPlotPDF <- downloadHandler(
        filename = function() {
          paste("PRC_Report_", Sys.Date(), ".pdf", sep = "")
        },
        content = function(file) {
          plot_info <- plot_data()
          if (!is.null(plot_info)) {
            # Generate the PDF report
            rmarkdown::render(input = "report_template.Rmd",
                              output_file = file,
                              params = list(
                                gene_s = plot_info$gene_s,
                                selected_scores = plot_info$selected_scores,
                                B_org = plot_info$B_org,
                                P_org = plot_info$P_org,
                                prcfiltered = plot_info$prcfiltered
                              ),
                              envir = new.env(parent = globalenv()))
        }
      }
    )
    }
  })
}
  