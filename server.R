# server.R
library(shiny)
library(shinyjs)
library(dplyr)
library(yogiroc)
library(maveLLR)
library(rsconnect)
library(rmarkdown)
library(knitr)
library(tinytex)
library(httr)
library(jsonlite)
library(DT)
library(ggplot2)
library(shinycssloaders)


server <- function(input, output, session) {
  plot_data <- reactiveVal(NULL)
  selected_variants <- reactiveVal(NULL)  # Store user-selected variants
  threshold_data <- reactiveVal(NULL) # for storing PRC thresholds - REFACTOR: include as part of prc data? - jumptag
  
  # Giant if/else block to handle separate logic for Main App and Advanced
  observe({
    if (input$navbar == "Basic") {
      # logic for Main App
      prcdata <- reactive({
        df <- read.table("preprocessed.csv", sep = ',', header = TRUE, stringsAsFactors = FALSE)
        return(df)
      })
      
      # Update gene names based on what exists
      observe({
        df <- prcdata()
        if (!is.null(df)) {
          gene_names <- unique(df$base__gene)
          updateSelectizeInput(session, "Main_gene", choices = gene_names, selected = character(0), server = TRUE)
        }
      })
      
      # Update scores - if a VEP is all NA, checkbox will not show up
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
      
      # walktag - for LLR VUS plots
      full_data <- reactive({
        req(input$Main_gene)  # Ensure a gene is selected
        
        # Read and filter full.csv for the selected gene
        df_full <- read.csv("full.csv", stringsAsFactors = FALSE)
        df_full$clinvar[df_full$clinvar == "N/A"] <- NA # Treat as actual NA
        df_filtered <- df_full %>% filter(base__gene == input$Main_gene & !is.na(clinvar))
        
        return(df_filtered)
      })
      
      
      # Show variant selection display
      showVariantSelectionModal <- function(df, gene_selected) {
        
        # Modified dataframe for display
        gene_display_df <- df %>% filter(base__gene == gene_selected) %>%
          select(base__gene, base__achange, clinvar, revstat, stars)
        
        showModal(modalDialog(
          title = paste("Select Variants for Gene:", gene_selected),
          tags$p("Would you like to de-select some variants? Click on the variants you do not want to include in the PRC. You can also do this later by downloading the CSV, removing variants yourself, and uploading it back."),
          DTOutput("variant_table"),
          footer = tagList(
            modalButton("Cancel"),
            actionButton("confirm_selection", "OK")
          ),
          easyClose = FALSE
        ))
        
        # Render the interactive table
        output$variant_table <- renderDT({
          datatable(
            gene_display_df, 
            selection = list(target = "row", selected = 1:nrow(gene_display_df)),  # All rows selected by default
            options = list(pageLength = 10, scrollX = TRUE)
          )
        }, server = TRUE)
      }
      
      # Generate PRC with user selection
      observeEvent(input$Main_plotButton, {
        gene_s <- input$Main_gene
        prcfiltered <- prcdata() %>%
          filter(base__gene == gene_s)
        if (!is.null(df)) {
          showVariantSelectionModal(prcfiltered, gene_s)
        }
      })
      
      # After confirmation from the modal
      observeEvent(input$confirm_selection, {
        # Wrap everything in a withProgress block
        withProgress(message = "Generating plots...", value = 0, {
          
          # Keep the existing “Please wait” modal
          showModal(modalDialog(
            title = "Please wait",
            "Generating plots...",
            footer = NULL,
            easyClose = FALSE
          ))
          
          # Step 1: Select rows
          #incProgress(0.1, detail = "Selecting rows")
          selected_rows <- input$variant_table_rows_selected
          gene_s <- input$Main_gene
          prcfiltered <- prcdata() %>%
            filter(base__gene == gene_s)
          
          if (!is.null(prcfiltered) && length(selected_rows) > 0) {
            # Step 2: Prepare selected data
            #incProgress(0.3, detail = "Preparing selected data")
            selected_df <- prcfiltered[selected_rows, ]
            selected_variants(selected_df)
            
            # Check that there is at least one P/LP and one B/LB - walktag
            if (!(any(selected_df$clinvar == "P/LP") && any(selected_df$clinvar == "B/LB"))) {
              removeModal() # remove "Please wait" modal if it's open
              showModal(modalDialog(
                title = "Error",
                "You must select at least one P/LP and one B/LB variant.",
                easyClose = TRUE,
                footer = modalButton("Close")
              ))
              return()
            }
            
            gene_s <- input$Main_gene
            exclude_common_variants <- input$Main_common_variant_filter
            selected_scores <- input$Main_scores
            
            # Renaming for consistency
            names(selected_df)[names(selected_df) == "VEP_varity_r"] <- "VARITY"
            names(selected_df)[names(selected_df) == "VEP_alphamissense__pathogenicity"] <- "AlphaMissense"
            names(selected_df)[names(selected_df) == "VEP_revel__score"] <- "REVEL"
            names(selected_df)[names(selected_df) == "gnomad__af"] <- "gnomAD_AF"
            
            # Common variant filter
            if (exclude_common_variants) {
              selected_df <- selected_df[is.na(selected_df$gnomAD_AF) | selected_df$gnomAD_AF <= 0.005, ]
            }
            selected_df <- selected_df[order(selected_df$clinvar), ]
            
            # Convert label to T/F
            selected_df <- selected_df %>% mutate(clinvar = ifelse(clinvar == "P/LP", TRUE, FALSE))
            
            # Count # of P and B
            P_org <- sum(selected_df$clinvar == TRUE & rowSums(!is.na(selected_df[selected_scores])) > 0)
            B_org <- sum(selected_df$clinvar == FALSE & rowSums(!is.na(selected_df[selected_scores])) > 0)
            
            tryCatch({
              #incProgress(0.5, detail = "Generating PRC plot")
              yrobj <- yr2(truth = selected_df[["clinvar"]], scores = selected_df[selected_scores], high = rep(TRUE, length(selected_scores)))
              
              # Added for threshold calculation - jumptag
              thresh_ranges <- calculate_thresh_range(yrobj, x = 0.9, balanced = TRUE)
              
              # Store thresholds for rendering
              threshold_data(thresh_ranges)
              
              # end threshold
              
              plot_data(list(
                yrobj = yrobj,
                lty_styles = c("dashed", "solid", "dashed")[1:length(selected_scores)],
                col_styles = c("purple", "cadetblue2", "orange")[1:length(selected_scores)],
                gene_s = gene_s,
                selected_scores = selected_scores,
                B_org = B_org,
                P_org = P_org,
                common_variant_filter = exclude_common_variants,
                prcfiltered = selected_df  # Save the filtered data for metadata
              ))
              output$Main_ErrorText <- renderText("")
              
            }, error = function(e) {
              plot_data(NULL)
              output$Main_ErrorText <- renderText("Error generating the plot.")
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
                    'Not enough data - must have at least one pathogenic and benign. Try selecting less VEPs or unselecting the common variant filter.
                  You can also add your own annotations (see Advanced Tab or download VUS).',
                    easyClose = TRUE,
                    footer = NULL
                  ))
                })
              }
            }, width = 600, height = 600, res = 72)
            
            # Render the threshold table - jumptag
            output$thresholdTableUI <- renderUI({
              req(threshold_data())
              
              tagList(
                h5("Threshold at 90% Balanced Precision"),  # Add the title only if the table exists
                tableOutput("thresholdTable")  
              )
            })
            
            output$thresholdTable <- renderTable({
              req(threshold_data())
              
              df <- threshold_data()
              
              # Round numeric columns to 3 decimals if numeric
              df <- as.data.frame(lapply(df, function(x) {
                if (is.numeric(x)) {
                  format(round(x, 3), nsmall = 3)
                } else {
                  x
                }
              }))
              df
            }, rownames = FALSE)
            
            # LLR Plot - walktag
            llr_tabs <- list()
            llr_scores <- intersect(selected_scores, c("VARITY", "REVEL", "AlphaMissense"))
            
            incProgress(0.5, detail = "Reading full LLR data (this might take a while)")
            
            # Retrieve filtered full data
            full_filtered <- full_data()

            # Renaming for consistency - walktag
            names(full_filtered)[names(full_filtered) == "VEP_varity_r"] <- "VARITY"
            names(full_filtered)[names(full_filtered) == "VEP_alphamissense__pathogenicity"] <- "AlphaMissense"
            names(full_filtered)[names(full_filtered) == "VEP_revel__score"] <- "REVEL"
            names(full_filtered)[names(full_filtered) == "gnomad__af"] <- "gnomAD_AF"
            
            for (score_idx in seq_along(llr_scores)) {
              score <- llr_scores[score_idx]
              incProgress(0.5+0.5/length(llr_scores), detail = paste("Generating LLR plot for:", score))
              
              posScores <- na.omit(setNames(
                selected_df[selected_df$clinvar == TRUE, score],
                selected_df[selected_df$clinvar == TRUE, "base__achange"]
              ))
              negScores <- na.omit(setNames(
                selected_df[selected_df$clinvar == FALSE, score],
                selected_df[selected_df$clinvar == FALSE, "base__achange"]
              ))
              
              if (length(posScores) > 0 & length(negScores) > 0) {
                llrObj <- buildLLR.kernel(posScores, negScores)
                
                full_filtered_copy <- full_filtered
                full_filtered_copy$llr <- llrObj$llr(full_filtered_copy[[score]]) # full_filtered_copy has llr values
                
                # Define breakpoints and labels for the categories
                breaks <- c(-Inf, -1.27, -0.32, 0.32, 0.64, 1.27, 2.54, Inf)
                labels <- c(
                  "benign_strong",   # (-Inf, -1.27]
                  "benign_support",  # (-1.27, -0.32]
                  "none",            # (-0.32, 0.32]
                  "patho_support",   # (0.32, 0.64]
                  "patho_moderate",  # (0.64, 1.27]
                  "patho_strong",    # (1.27, 2.54]
                  "patho_vstrong"    # (2.54, Inf)
                )
                
                # Categorize llr based on the defined bins
                full_filtered_copy$category <- cut(full_filtered_copy$llr, breaks=breaks, labels=labels, include.lowest=TRUE)
                
                full_filtered_copy <- full_filtered_copy %>% filter(!is.na(category)) # filter NA - due to score and llr being NA
                
                # Define thresholds and search range for crossings
                llrTs <- llrThresholds(optiLLR(0.1))
                x_range <- range(c(posScores, negScores))
                
                # Find where the LLR function crosses the thresholds
                crossings_df <- findLLRcrossings(llrObj$llr, llrTs, x_range)
                
                # Create the stacked bar plot of clinvar vs category
                # Assign colors as requested:
                # dark blue for benign_strong, red for patho_vstrong, grey for none
                # other categories lighter shades in between
                category_colors <- c(
                  "benign_strong"   = "darkblue",
                  "benign_support"  = "lightblue",
                  "none"            = "grey",
                  "patho_support"   = "#FFC0CB",   # light pink
                  "patho_moderate"  = "#FF9999",   # medium pink
                  "patho_strong"    = "#FF6666",   # darker pink
                  "patho_vstrong"   = "red"
                )
                
                plot_stack_id <- paste0("Main_Stacked_", score)
                
                # Preserve current iteration's variables, or else it would be same plot
                local({
                  # Make local copies of the variables to avoid referencing loop variables
                  score_copy <- score
                  posScores_copy <- posScores
                  negScores_copy <- negScores
                  llrObj_copy <- llrObj
                  full_filtered_copy_local <- full_filtered_copy
                  crossings_df_copy <- crossings_df
                  
                  plot_id <- paste0("Main_LLRPlot_", score_copy)
                  table_id <- paste0("Main_LLRCrossingsTable_", score_copy)
                  
                  # -----------------------------
                  # 1) DEFINE DOWNLOAD BUTTON IDS
                  # -----------------------------
                  download_llr_png_id <- paste0("download_", score_copy, "_llr_png")
                  download_bar_svg_id <- paste0("download_", score_copy, "_bar_svg")
                  download_both_pdf_id <- paste0("download_", score_copy, "_both_pdf")
                  download_csv_id      <- paste0("download_", score_copy, "_csv")
                  
                  output[[plot_id]] <- renderPlot({
                    tryCatch({
                      drawDensityLLR_fixedRange(
                        full_filtered_copy_local[[score_copy]], 
                        llrObj_copy$llr, 
                        llrObj_copy$posDens, 
                        llrObj_copy$negDens, 
                        posScores_copy, 
                        negScores_copy
                      )
                    }, error = function(e) {
                      showModal(modalDialog(
                        title = 'Error',
                        'Not enough data',
                        easyClose = TRUE,
                        footer = NULL
                      ))
                    })
                  }, width = 500, height = 600, res = 72)
                  
                  # Render the stacked bar chart
                  output[[plot_stack_id]] <- renderPlot({
                    # Make sure we have factor ordering if needed
                    full_filtered_copy_local$clinvar <- factor(full_filtered_copy_local$clinvar, 
                                                               levels = c("P/LP", "B/LB", "VUS", "Conflicting"))
                    ggplot(full_filtered_copy_local, aes(x=clinvar, fill=category)) +
                      geom_bar(position=position_stack(reverse=TRUE)) + # reverse to put pathogenic at top like the LLR
                      scale_fill_manual(values=category_colors) +
                      theme_minimal() +
                      labs(x="ClinVar Category", y="Count", fill="LLR Category") +
                      theme(axis.text.x = element_text(angle=45, hjust=1))
                  }, width = 300, height = 500, res = 72)
                  
                  # Render the table of LLR threshold crossings
                  output[[table_id]] <- renderTable({
                    crossings_df_copy
                  })
                  
                  # -----------------------------------
                  # 3) DEFINE DOWNLOAD HANDLERS (NEW)
                  # -----------------------------------
                  
                  # 3a) Download LLR (PNG)
                  output[[download_llr_png_id]] <- downloadHandler(
                    filename = function() { paste0(score_copy, "_LLRPlot.png") },
                    content = function(file) {
                      # Use a PNG device and re-draw the same LLR plot
                      png(file, width = 500, height = 600, res = 72)
                      drawDensityLLR_fixedRange(
                        full_filtered_copy_local[[score_copy]],
                        llrObj_copy$llr,
                        llrObj_copy$posDens,
                        llrObj_copy$negDens,
                        posScores_copy,
                        negScores_copy
                      )
                      dev.off()
                    }
                  )
                  
                  # 3b) Download bar plot (SVG)
                  output[[download_bar_svg_id]] <- downloadHandler(
                    filename = function() {
                      paste0(score_copy, "_BarPlot.svg")
                    },
                    content = function(file) {
                      # Open the SVG device
                      svg(file, width = 5, height = 5)
                      
                      # Construct the ggplot as usual
                      p <- ggplot(full_filtered_copy_local, aes(x=clinvar, fill=category)) +
                        geom_bar(position=position_stack(reverse=TRUE)) +
                        scale_fill_manual(values=category_colors) +
                        theme_minimal() +
                        labs(x="ClinVar Category", y="Count", fill="LLR Category") +
                        theme(axis.text.x = element_text(angle=45, hjust=1))
                      
                      # Print the ggplot to render onto the SVG device
                      print(p)
                      
                      # Close the device
                      dev.off()
                    }
                  )
                  
                  # 3c) Download both (PDF)
                  output[[download_both_pdf_id]] <- downloadHandler(
                    filename = function() { paste0(score_copy, "_LLR_and_Bar.pdf") },
                    content = function(file) {
                      # Create a 2-page PDF: first page = LLR plot, second page = bar plot
                      pdf(file, width = 6, height = 6)
                      
                      # Page 1: LLR
                      drawDensityLLR_fixedRange(
                        full_filtered_copy_local[[score_copy]],
                        llrObj_copy$llr,
                        llrObj_copy$posDens,
                        llrObj_copy$negDens,
                        posScores_copy,
                        negScores_copy
                      )
                      
                      # Page 2: bar plot
                      plot.new()  
                      p <- ggplot(full_filtered_copy_local, aes(x=clinvar, fill=category)) +
                        geom_bar(position=position_stack(reverse=TRUE)) +
                        scale_fill_manual(values=category_colors) +
                        theme_minimal() +
                        labs(x="ClinVar Category", y="Count", fill="LLR Category") +
                        theme(axis.text.x = element_text(angle=45, hjust=1))
                      print(p)
                      
                      dev.off()
                    }
                  )
                  
                  # 3d) Download CSV of variants w/ LLR + category
                  output[[download_csv_id]] <- downloadHandler(
                    filename = function() { paste0(score_copy, "_LLR_data.csv") },
                    content = function(file) {
                      # This CSV includes all variants (including VUS/conflicting) 
                      # that remain in full_filtered_copy_local, along with LLR and category
                      write.csv(full_filtered_copy_local, file, row.names = FALSE)
                    }
                  )
                  
                  # -----------------------------------
                  # 4) BUILD THE UI FOR THIS TAB
                  # -----------------------------------
                  llr_tabs[[length(llr_tabs) + 1]] <<- tabPanel(
                    score_copy,
                    fluidRow(
                      column(
                        5,
                        # LLR plot
                        plotOutput(plot_id, width = "500px", height = "600px")
                      ),
                      column(
                        4, offset = 3,
                        # Bar plot
                        plotOutput(plot_stack_id, width = "300px", height = "500px"),
                        
                        # Download buttons stacked vertically
                        tags$div(
                          style = "margin-top: 40px; text-align: left;",
                          downloadButton(download_llr_png_id,   "Download LLR (PNG)"),
                          tags$br(),
                          downloadButton(download_bar_svg_id,   "Download bar plot (SVG)"),
                          tags$br(),
                          downloadButton(download_both_pdf_id,  "Download both (PDF)"),
                          tags$br(),
                          downloadButton(download_csv_id,       "Download all variants with LLRs (CSV)")
                        )
                      )
                    ),
                    tableOutput(table_id)
                  )
                })
              }
            }
            
            if (length(llr_tabs) > 0) {
              output$Main_LLRTabs <- renderUI({
                do.call(tabsetPanel, c(id="llr_tabs", llr_tabs))
              })
            } else {
              output$Main_LLRTabs <- renderUI({
                "No LLR plots available"
              })
            }
            
            # Final step: completed
            incProgress(1.0, detail = "Done")
            removeModal()
            
            if (P_org < 11 || B_org < 11) {
              showModal(modalDialog(
                title = "Warning",
                "There are fewer than 11 Pathogenic/Likely Pathogenic or Benign/Likely Benign variants. 
           Please use extra caution when interpreting these plots.",
                easyClose = TRUE,
                footer = modalButton("Close")
              ))
            }
            
          } else {
            removeModal() # In case wait modal is still visible
            showModal(modalDialog(
              title = "Error",
              "Please select at least two variants.",
              easyClose = TRUE,
              footer = modalButton("Close")
            ))
          }
        })
      })
      
      output$Main_download_buttons <- renderUI({
        req(plot_data())  # Only show if plot_data has been generated
        
        tagList(
          actionButton("Main_helpButton", "Plot Explanation", class = "btn-info"),
          helpText(HTML("<span style='color:black;'><strong>Download Options: </strong></span>")),
          downloadButton("Main_downloadPlotPNG", "Download PRC Plot as PNG"),
          downloadButton("Main_downloadPlotPDF", HTML("Download PRC Plot and <br>Metadata as PDF")), 
          downloadButton("Main_downloadCSV", HTML("Download Variants Used as <br>CSV")),
          div(
            style = "display: inline-flex; align-items: center;",
            downloadButton("Main_downloadCSV_VUS", HTML("Download CSV with VUS")),
            actionLink("Main_helpButton_VUS", label = NULL, icon = icon("question-circle"), style = "margin-left: 5px;")
          )
        )
      })
      
      # Download PNG logic for Main App
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
      
      # Download PDF with metadata logic for Main App
      output$Main_downloadPlotPDF <- downloadHandler(
        filename = function() {
          paste("PRC_Report_", input$Main_gene, ".pdf", sep = "")
        },
        content = function(file) {
          plot_info <- plot_data()
          if (!is.null(plot_info)) {
            rmarkdown::render(input = "report_template.Rmd",
                              output_file = file,
                              params = list(
                                gene_s = plot_info$gene_s,
                                selected_scores = plot_info$selected_scores,
                                common_variant_filter = plot_info$common_variant_filter,
                                B_org = plot_info$B_org,
                                P_org = plot_info$P_org,
                                prcfiltered = plot_info$prcfiltered
                              ),
                              envir = new.env(parent = globalenv()))
          }
        }
      )
      
      # Download variants used as CSV
      output$Main_downloadCSV <- downloadHandler(
        filename = function() {
          paste("PRC_data_", input$Main_gene, ".csv", sep = "")
        },
        
        content = function(file) {
          selected_df <- selected_variants()
          if (!is.null(selected_df)) {
            write.csv(selected_df, file, row.names = FALSE)
          }
        }
      )
      
      # Download CSV with VUS
      output$Main_downloadCSV_VUS <- downloadHandler(
        filename = function() {
          paste("PRC_data_VUS_", input$Main_gene, ".csv", sep = "")
        },
        content = function(file) {
          df_VUS <- read.csv("full.csv", stringsAsFactors = FALSE)
          gene_s <- input$Main_gene
          
          prcfiltered <- df_VUS %>%
            filter(base__gene == gene_s) %>%
            arrange(clinvar) # Order alphabetical
          
          write.csv(prcfiltered, file, row.names = FALSE)
        }
      )
      
      # VUS explanation help button
      observeEvent(input$Main_helpButton_VUS, {
        showModal(modalDialog(
          title = "What does this data include?",
          HTML("In addition to P/LP or B/LB variants, this download includes VUS and Conflicting variants. 
          VUS stands for Variants of Uncertain Significance, and Conflicting means conflicting interpretations of pathogenicity. 
          These variants are not included in the PRC generation as they are not classified as Pathogenic or Benign. 
          However, by downloading this CSV, you can override the VUS or Conflicting annotations and reupload it as a custom reference set in Advanced Mode."),
          easyClose = TRUE,
          footer = modalButton("Close")
        ))
      }, ignoreInit = TRUE)
      
      # Plot explanation help button
      observeEvent(input$Main_helpButton, {
        showModal(modalDialog(
          title = "Plot Explanation",
          HTML("<p><b>What is balanced precision?</b></p>
            <p>Balanced precision is useful in situations where the class distribution of the test set is imbalanced. In other words, it is the precision that would have been expected had the proportion of 'P/LP' examples been balanced (equal to 50%).</p>"),
          p("Definition ", 
            a("here", href = "https://doi.org/10.1016/j.ajhg.2021.08.012", target = "_blank")),
          HTML("<p><b>What are R90BP and AUBPRC?</b></p>
          <p><b>R90BP:</b> The recall achieved at a stringent (90%) balanced precision threshold </p>
          <p><b>AUBPRC:</b> The area under the BPRC curve</p>"),
          easyClose = TRUE,
          footer = modalButton("Close")
        ))
      })
      
    } else {
      # Logic for Advanced
      rv <- reactiveValues(
        state = reactiveValues(gnomad_filter_inserted = FALSE)
      )
      
      # Define reactive values
      rv$prcdata_own <- reactiveVal(NULL)
      rv$prcdata_fetch <- reactiveVal(NULL)
      rv$plot_data <- reactiveVal(NULL)
      rv$variant_data_df <- reactiveVal(NULL)
      
      # Observe 'input$input_type' to reset reactive values and inputs when user switches option
      observeEvent(input$input_type, {
        rv$plot_data(NULL)
        rv$prcdata_fetch(NULL)
        rv$prcdata_own(NULL)
        rv$gnomad_filter_inserted <- FALSE
        
        updateCheckboxGroupInput(session, "scores", choices = NULL, selected = NULL)
        removeUI(selector = "#gnomad_filter_wrapper", immediate = TRUE)
      })
      
      # Logic for 'own' input_type
      observeEvent(input$file_full, {
        if (input$input_type == "own") {
          df <- tryCatch({
            req(input$file_full)
            df <- read.csv(input$file_full$datapath, stringsAsFactors = FALSE)
            
            # Check for necessary columns
            if (!all(c("base__gene", "clinvar") %in% colnames(df))) {
              showModal(modalDialog(
                title = "Error",
                "The uploaded CSV must contain the columns 'base__gene' and 'clinvar'.",
                easyClose = TRUE,
                footer = modalButton("Close")
              ))
              return(NULL)
            }
            
            # Consistency for mandatory columns
            colnames(df)[colnames(df) == "gnomad__af"] <- "gnomAD_AF"
            colnames(df)[colnames(df) == "clinvar"] <- "clinvar"
            
            df
          }, error = function(e) {
            showModal(modalDialog(
              title = "Error",
              "There was an error reading the uploaded file. Please ensure it is a valid CSV file.",
              easyClose = TRUE,
              footer = modalButton("Close")
            ))
            NULL
          })
          rv$prcdata_own(df)  # Set the reactiveVal
        }
      })
      
      # Guide on how to format user-inputted csv
      observeEvent(input$upload_guide, {
        showModal(modalDialog(
          title = "Reference Set Format Information",
          HTML("Please ensure your own reference set is a CSV file.<br><br>
            Mandatory columns:<br>
            <b>base__gene:</b> Gene name<br>
            <b>clinvar:</b> Variant interpretation (B/LB or P/LP)<br><br>
            
            Optional columns: <br>
            <b>gnomad_af:</b> GnomAD allele frequency - for filtering out common variants<br><br>
            
            For any predictors or custom scores that you would like to include, put 'VEP_' before the name of the column. For example:<br>
            <b>VEP_alphamissense:</b> AlphaMissense score<br>
            <b>VEP_custom:</b> Additional scores you would like to include<br><br>"),
          easyClose = FALSE,
          footer = tagList(
            downloadButton("download_template_own", "Download CSV Template"),
            modalButton("Close")
          )
        ))
      }, ignoreInit = TRUE)
      
      # Download template for own reference set
      output$download_template_own <- downloadHandler(
        filename = function() {
          "reference_set_template.csv"
        },
        content = function(file) {
          file.copy("own_template.csv", file)
        }
      )
      
      # Download OC reference template
      output$download_template_oc <- downloadHandler(
        filename = function() {
          "reference_set_template.csv"
        },
        content = function(file) {
          file.copy("oc_template.csv", file)
        }
      )
      
      # Logic for 'fetch' input_type
      observeEvent(input$fetchButton, {
        if (input$input_type == "fetch") {
          # Initialize empty list to store results
          all_results <- list()
          
          tryCatch({
            req(input$file_fetch)
            df <- read.csv(input$file_fetch$datapath, stringsAsFactors = FALSE)
            
            # Ensure 'chrom' column in df always has the "chr" prefix for consistency
            df <- df %>%
              mutate(chrom = ifelse(grepl("^chr", as.character(chrom)), as.character(chrom), paste0("chr", chrom)))
            
            # Show progress while fetching data
            withProgress(message = 'Fetching Variant Data', value = 0, {
              for (i in 1:nrow(df)) {
                
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
                  "&annotators=", paste(sub("varity_er", "varity_r", tolower(input$Fetch_scores)), collapse = ",") 
                )
                
                response <- GET(api_url)
                
                result <- fromJSON(content(response, "text"), flatten = TRUE)
                
                # Extract gene and achange into corresponding columns
                gene <- ifelse(!is.null(result$crx$hugo), result$crx$hugo, NA)
                achange <- ifelse(!is.null(result$crx$achange), result$crx$achange, NA)
                
                # Initialize result_list for required fields
                result_list <- list(
                  chrom = chrom,
                  pos = pos,
                  ref_base = ref,
                  alt_base = alt
                )
                
                # Add gene and achange columns only if they don't already exist in df
                if (!"base__gene" %in% colnames(df)) {
                  result_list[["base__gene"]] <- gene
                }
                if (!"base__achange" %in% colnames(df)) {
                  result_list[["base__achange"]] <- achange
                }
                
                # Mapping possible scores from API result to data
                score_mappings <- list(
                  "gnomad" = list(
                    column_name = "gnomAD_AF",
                    value = if (!is.null(result$gnomad$af)) result$gnomad$af else NA
                  ),
                  "varity_r" = list(
                    column_name = "VEP_VARITY_R",
                    value = if (!is.null(result$varity_r$varity_r)) result$varity_r$varity_r else NA
                  ),
                  "varity_er" = list(
                    column_name = "VEP_VARITY_ER",
                    value = if (!is.null(result$varity_r$varity_er)) result$varity_r$varity_er else NA
                  ),
                  "revel" = list(
                    column_name = "VEP_REVEL",
                    value = if (!is.null(result$revel$score)) result$revel$score else NA
                  ),
                  "alphamissense" = list(
                    column_name = "VEP_AlphaMissense",
                    value = if (!is.null(result$alphamissense$am_pathogenicity)) result$alphamissense$am_pathogenicity else NA
                  ),
                  "provean" = list(
                    column_name = "VEP_Provean",
                    value = if (!is.null(result$provean$rankscore)) result$provean$rankscore else NA
                  ),
                  "sift" = list(
                    column_name = "VEP_SIFT",
                    value = if (!is.null(result$sift$rankscore)) result$sift$rankscore else NA
                  ),
                  "polyphen2" = list(
                    column_name = "VEP_PolyPhen2",
                    value = if (!is.null(result$polyphen2$hvar_rank)) result$polyphen2$hvar_rank else NA
                  ),
                  # Convert clinvar to P or B, ignore VUS or conflicting
                  "clinvar" = list(
                    column_name = "clinvar",
                    value = if (!is.null(result$clinvar$sig)) {
                      clinvar_sig <- tolower(result$clinvar$sig)
                      if (grepl("benign", clinvar_sig)) {
                        "B/LB"
                      } else if (grepl("conflicting", clinvar_sig)) {
                        "Conflicting"
                      } else if (grepl("pathogenic", clinvar_sig)) {
                        "P/LP"
                      } else if (grepl("uncertain", clinvar_sig)) {
                        "VUS"
                      } else {
                        NA
                      }
                    } else {
                      NA
                    }
                  )
                )
                
                for (score in input$Fetch_scores) {
                  if (score %in% names(score_mappings)) {
                    mapping <- score_mappings[[score]]
                    if (mapping$column_name %in% colnames(df)) {
                      # Update column if it already exists
                      df[[mapping$column_name]][i] <- mapping$value
                    } else {
                      # Add as new column
                      result_list[[mapping$column_name]] <- mapping$value
                    }
                  }
                }
                
                result_df <- as.data.frame(result_list, stringsAsFactors = FALSE)
                all_results[[i]] <- result_df
                
                # Show progress bar
                incProgress(1 / nrow(df))
                showNotification(paste(i, "/", nrow(df), "variants fetched"), duration = 3, type = "message")
              }
            })
            
            # Combine all results into a data frame
            variant_data_df <- do.call(rbind, all_results)
            
            variant_data_df <- dplyr::left_join(df, variant_data_df, by = c("chrom", "pos", "ref_base", "alt_base"))
            
            output$errorText <- renderText("")
            rv$prcdata_fetch(variant_data_df)  # Set the reactiveVal
            rv$variant_data_df <- variant_data_df

          }, error = function(e) {
            showModal(modalDialog(
              title = "Error",
              "An error occurred while fetching data.",
              easyClose = TRUE,
              footer = modalButton("Close")
            ))
            rv$prcdata_fetch(NULL)
          })
        }
      })
      
      
      # Define 'prcdata' reactive depending on where it comes from (own or fetched)
      prcdata <- reactive({
        if (input$input_type == "own") {
          rv$prcdata_own()
        } else if (input$input_type == "fetch") {
          rv$prcdata_fetch()
        } else {
          NULL
        }
      })
      
      # Observe when prcdata changes to add gnomAD checkbox
      observeEvent(prcdata(), {
        df <- prcdata()
        if (is.null(df)) return()
        
        # Check for gnomAD columns, case-insensitive
        colnames_lower <- tolower(colnames(df))
        gnomad_columns <- grep("gnomad", colnames_lower, value = TRUE)
        
        # If gnomAD columns are found and the checkbox is not yet added
        if (length(gnomad_columns) > 0 && isFALSE(rv$gnomad_filter_inserted)) {
          # Remove any existing checkbox to prevent duplicates
          removeUI(selector = "#gnomad_filter_wrapper", immediate = TRUE)
          
          # Insert the checkbox
          insertUI(
            selector = "#scores",
            where = "beforeBegin",
            ui = div(id = "gnomad_filter_wrapper",
                     checkboxInput("common_variant_filter", 
                                   "Exclude Common Variants (gnomAD AF > 0.005)", 
                                   value = TRUE)
            )
          )
          rv$gnomad_filter_inserted <- TRUE
          
        } else if (length(gnomad_columns) == 0 && isTRUE(rv$gnomad_filter_inserted)) {
          # Remove the checkbox if no gnomAD columns are found and it is still there
          removeUI(selector = "#gnomad_filter_wrapper", immediate = TRUE)
          rv$gnomad_filter_inserted <- FALSE
        }
        
        # Remove "VEP_" prefix from predictor column names for display
        predictor_columns <- colnames(df)[grepl("^VEP_", colnames(df))]
        predictor_columns <- gsub("^VEP_", "", predictor_columns)
        # Update checkboxGroupInput with trimmed predictor columns
        updateCheckboxGroupInput(session, "scores", choices = predictor_columns, selected = predictor_columns)
      })
      
      # Observe 'plotButton'
      observeEvent(input$plotButton, {
        df <- prcdata()
        if (is.null(df) || nrow(df) == 0) {
          output$errorText <- renderText("Not enough rows to generate the PRC plot.")
          return()
        }
        
        gene_s <- ifelse("base__gene" %in% colnames(df), df$base__gene[1], "Custom Gene")
        exclude_common_variants <- input$common_variant_filter
        selected_scores <- input$scores
        
        if (isTRUE(exclude_common_variants) && "gnomAD_AF" %in% colnames(df)) { # walktag - added isTRUE
          df <- df[is.na(df$gnomAD_AF) | df$gnomAD_AF <= 0.005, ]
        }
        
        # Remove "VEP_" prefix from predictor column names
        colnames(df) <- gsub("^VEP_", "", colnames(df))
        
        df <- df[order(df$clinvar), ]
        
        prcfiltered <- df %>%
          filter(rowSums(!is.na(df[selected_scores])) > 0) %>%
          mutate(clinvar = ifelse(clinvar == "P/LP", TRUE, ifelse(clinvar == "B/LB", FALSE, NA)))
        
        prcfiltered <- prcfiltered[!is.na(prcfiltered$clinvar), ]
        
        P_org <- sum(prcfiltered$clinvar == TRUE)
        B_org <- sum(prcfiltered$clinvar == FALSE)
        
        tryCatch({
          yrobj <- yr2(truth = prcfiltered[["clinvar"]], scores = prcfiltered[selected_scores], high = rep(TRUE, length(selected_scores)))
          
          # Extra colors for lines more than 3
          # Exclude white and near-white colors
          available_colors <- grDevices::colors()[!grepl("white|ivory|seashell|snow|honeydew|azure|aliceblue|mintcream|ghostwhite", grDevices::colors(), ignore.case = TRUE)]
          
          # Generate random colors for additional predictors if more than three
          num_scores <- min(length(selected_scores), 3)  # Only take a maximum of three predictors for specific styles
          extra_colors <- if (length(selected_scores) > 3) {
            sample(available_colors, length(selected_scores) - 3)
          } else {
            NULL
          }
          
          # Assign styles for the first three predictors, followed by additional random colors if needed
          rv$plot_data(list(
            yrobj = yrobj,
            lty_styles = c("dashed", "solid", "dashed")[1:num_scores],
            col_styles = c("purple", "cadetblue2", "orange")[1:num_scores],
            gene_s = gene_s,
            selected_scores = selected_scores,
            B_org = B_org,
            P_org = P_org,
            common_variant_filter = exclude_common_variants,
            prcfiltered = prcfiltered,
            extra_colors = extra_colors  # Store extra colors for additional predictors if any
          ))
          
          output$errorText <- renderText("")
          
        }, error = function(e) {
          rv$plot_data(NULL)
          output$errorText <- renderText("Error generating plot.")
        })
        
        output$prcPlot <- renderPlot({
          plot_info <- rv$plot_data()
          if (!is.null(plot_info)) {
            tryCatch({
              colors_to_use <- c(plot_info$col_styles, plot_info$extra_colors)
              draw.prc(
                plot_info$yrobj,
                lty = c(plot_info$lty_styles, rep(c("solid", "dashed"), length(plot_info$extra_colors))),
                col = colors_to_use,
                lwd = 2,
                balanced = TRUE,
                main = paste0(plot_info$gene_s, " PRCs for ", paste(plot_info$selected_scores, collapse = ", "))
              )
              abline(h = 90, lty = "dashed")
              legend("left", legend = c(paste("# of Pathogenic and Likely Pathogenic:", plot_info$P_org), paste("# of Benign and Likely Benign:", plot_info$B_org)), pch = 15, bty = "n")
            }, error = function(e) {
              showModal(modalDialog(
                title = 'Error',
                'Not enough data - must have at least one pathogenic and benign.',
                easyClose = TRUE,
                footer = NULL
              ))
            })
          }
        }, width = 600, height = 600, res = 72)
        
        if (P_org < 11 || B_org < 11) {
          showModal(modalDialog(
            title = "Warning",
            "There are fewer than 11 Pathogenic/Likely Pathogenic or Benign/Likely Benign variants. 
           Please use extra caution when interpreting these plots.",
            easyClose = TRUE,
            footer = modalButton("Close")
          ))
        }
        
      })
      
      # Download logic 
      output$download_buttons <- renderUI({
        req(rv$plot_data())  # Only show if plot has been generated
        
        tagList(
          actionButton("helpButton", "Plot Explanation", class = "btn-info"),
          helpText(HTML("<span style='color:black;'><strong>Download Options: </strong></span>")),
          downloadButton("downloadPlotPNG", "Download PRC Plot as PNG"),
          downloadButton("downloadPlotPDF", HTML("Download PRC Plot and <br>Metadata as PDF")),
          downloadButton("downloadCSV", HTML("Download Variants Used as <br>CSV")),
          if (is.data.frame(rv$variant_data_df) && nrow(rv$variant_data_df) > 0) { # Only if getting data from fetch
            div(
              style = "display: inline-flex; align-items: center;",
              downloadButton("downloadCSV_VUS", HTML("Download CSV with VUS")),
              actionLink("helpButton_VUS", label = NULL, icon = icon("question-circle"), style = "margin-left: 5px;")
            )
          } else {
            NULL
          }
        )
      })
      
      output$downloadPlotPNG <- downloadHandler(
        filename = function() {
          paste("PRC_plot_", Sys.Date(), ".png", sep = "")
        },
        content = function(file) {
          plot_info <- rv$plot_data()
          if (!is.null(plot_info)) {
            num_scores <- length(plot_info$selected_scores)
            
            available_colors <- grDevices::colors()[!grepl("white|ivory|seashell|snow|honeydew|azure|aliceblue|mintcream|ghostwhite", grDevices::colors(), ignore.case = TRUE)]
            
            col_styles <- if (num_scores <= 3) {
              c("purple", "cadetblue2", "orange")[1:num_scores]
            } else {
              sample(available_colors, num_scores)
            }
            
            lty_styles <- rep(c("solid", "dashed"), length.out = num_scores)
            
            png(file, width = 6, height = 6, units = "in", res = 72)
            draw.prc(
              plot_info$yrobj,
              lty = lty_styles,
              col = col_styles,
              lwd = 2,
              balanced = TRUE,
              main = paste0(plot_info$gene_s, " PRCs for ", paste(plot_info$selected_scores, collapse = ", "))
            )
            abline(h = 90, lty = "dashed")
            legend(
              "left",
              legend = c(
                paste("# of Pathogenic and Likely Pathogenic:", plot_info$P_org),
                paste("# of Benign and Likely Benign:", plot_info$B_org)
              ),
              pch = 15,
              bty = "n"
            )
            dev.off()
          }
        }
      )
      
      
      output$downloadPlotPDF <- downloadHandler(
        filename = function() {
          paste("PRC_Report_", Sys.Date(), ".pdf", sep = "")
        },
        content = function(file) {
          plot_info <- rv$plot_data()
          if (!is.null(plot_info)) {
            rmarkdown::render(input = "report_template_custom.Rmd",
                              output_file = file,
                              params = list(
                                gene_s = plot_info$gene_s,
                                selected_scores = plot_info$selected_scores,
                                common_variant_filter = plot_info$common_variant_filter,
                                B_org = plot_info$B_org,
                                P_org = plot_info$P_org,
                                prcfiltered = plot_info$prcfiltered
                              ),
                              envir = new.env(parent = globalenv()))
          }
        }
      )
      
      output$downloadCSV <- downloadHandler(
        filename = function() {
          paste("PRC_data_", Sys.Date(), ".csv", sep = "")
        },
        content = function(file) {
          plot_info <- rv$plot_data()
          if (!is.null(plot_info) && !is.null(plot_info$prcfiltered) && nrow(plot_info$prcfiltered) > 0) {
            prcfiltered <- plot_info$prcfiltered
            prcfiltered <- prcfiltered %>%
              mutate(clinvar = ifelse(clinvar == TRUE, "P/LP", "B/LB"))
            write.csv(prcfiltered, file, row.names = FALSE)
          }
        }
      )
      
      output$downloadCSV_VUS <- downloadHandler(
        filename = function() {
          paste("PRC_data_VUS", Sys.Date(), ".csv", sep = "")
        },
        content = function(file) {
          req(rv$variant_data_df)
          write.csv(rv$variant_data_df, file, row.names = FALSE)
        }
      )
      
      observeEvent(input$helpButton_VUS, {
        showModal(modalDialog(
          title = "What does this data include?",
          HTML("In addition to P/LP or B/LB variants, this download includes VUS and Conflicting variants. 
          VUS stands for Variants of Uncertain Significance, and Conflicting means conflicting interpretations of pathogenicity. 
          These variants are not included in the PRC generation as they are not classified as Pathogenic or Benign. 
          However, by downloading this CSV, you can override the VUS or Conflicting annotations and reupload it as a custom reference set in Advanced Mode."),
          easyClose = TRUE,
          footer = modalButton("Close")
        ))
      }, ignoreInit = TRUE)
      
      # ObserveEvent for helpButton
      observeEvent(input$helpButton, {
        showModal(modalDialog(
          title = "Plot Explanation",
          HTML("<p><b>What is balanced precision?</b></p>
      <p>Balanced precision is useful in situations where the class distribution is imbalanced. In other words, it is the precision that would have been expected had the proportion of positive examples been balanced (equal to 50%).</p>"),
          p("Definition ", 
            a("here", href = "https://doi.org/10.1016/j.ajhg.2021.08.012", target = "_blank")),
          HTML("<p><b>What are R90BP and AUBPRC?</b></p>
      <p><b>R90BP:</b> The recall achieved at a stringent (90%) balanced precision threshold </p>
      <p><b>AUBPRC:</b> The area under the BPRC curve</p>"),
          easyClose = TRUE,
          footer = modalButton("Close")
        ))
      })
    }
    
  })
}
