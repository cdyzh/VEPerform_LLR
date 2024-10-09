# server.R
library(shiny)
library(shinyjs)
library(dplyr)
library(yogiroc)
library(rsconnect)
library(rmarkdown)
library(knitr)
library(tinytex) # test

# Define server logic for the application
server <- function(input, output, session) {
  plot_data <- reactiveVal(NULL)
  
  # Observe the Proceed button to switch tabs
  observeEvent(input$proceed_button, {
    if (input$intro_data_source == "upload") {
      updateNavbarPage(session, "navbar", selected = "Upload")
    } else {
      updateNavbarPage(session, "navbar", selected = "Main App")
    }
  })
  
  # Giant if/else block to handle separate logic for Main App and Upload tabs
  observe({
    if (input$navbar == "Main App") {
      # prcdata logic for Main App
      prcdata <- reactive({
        if (input$data_source == "upload" && !is.null(input$file_gene_variant)) {
          req(input$file_gene_variant)
          df <- read.csv(input$file_gene_variant$datapath, stringsAsFactors = FALSE)
          
          if (ncol(df) != 2) {
            showModal(modalDialog(
              title = "Error",
              "The uploaded dataset must have exactly two columns: base__hugo and base__achange.",
              easyClose = TRUE,
              footer = modalButton("Close")
            ))
            reset("file_gene_variant")
            updateFileInput(session, "file_gene_variant", value = NULL) # TODO: error handling
            return(NULL)
          }
          
          colnames(df) <- c("base__hugo", "base__achange")
          full_df <- read.csv("preprocessed_id.csv", stringsAsFactors = FALSE)
          df <- merge(df, full_df, by = c("base__hugo", "base__achange"))
        } else {
          df <- read.table("preprocessed_id.csv", sep = ',', header = TRUE, stringsAsFactors = FALSE)
        }
        return(df)
      })
      
      observe({
        df <- prcdata()
        if (!is.null(df)) {
          gene_names <- unique(df$base__hugo)
          updateSelectizeInput(session, "Main_gene", choices = gene_names, selected = gene_names[1], server = TRUE)
        }
      })
      
      observe({
        req(input$Main_gene)
        df <- prcdata()
        gene_data <- df %>% filter(base__hugo == input$gene)
        
        available_scores <- c()
        if (any(!is.na(gene_data$varity_r__varity_r))) {
          available_scores <- c(available_scores, "VARITY")
        }
        if (any(!is.na(gene_data$alphamissense__am_pathogenicity))) {
          available_scores <- c(available_scores, "AlphaMissense")
        }
        if (any(!is.na(gene_data$revel__score))) {
          available_scores <- c(available_scores, "REVEL")
        }
        
        updateCheckboxGroupInput(session, "scores", choices = available_scores, selected = available_scores)
      })
      
      # Plotting logic for Main App
      observeEvent(input$Main_plotButton, {
        df <- prcdata()
        gene_s <- input$Main_gene
        exclude_common_variants <- input$Main_common_variant_filter
        selected_scores <- input$Main_scores
        print("Main App: begin plot") # DEBUG
        
        names(df)[names(df) == "varity_r__varity_r"] <- "VARITY"
        names(df)[names(df) == "alphamissense__am_pathogenicity"] <- "AlphaMissense"
        names(df)[names(df) == "revel__score"] <- "REVEL"
        names(df)[names(df) == "gnomad__af"] <- "gnomAD_AF"
        
        if (exclude_common_variants) {
          df <- df[is.na(df$gnomAD_AF) | df$gnomAD_AF <= 0.005, ]
        }
        
        df <- df[order(df$classification), ]
        
        prcfiltered <- df %>%
          filter(base__hugo == gene_s)
        print(prcfiltered) #DEBUG
        
        B_org <- sum(prcfiltered$classification == TRUE & rowSums(!is.na(prcfiltered[selected_scores])) > 0)
        P_org <- sum(prcfiltered$classification == FALSE & rowSums(!is.na(prcfiltered[selected_scores])) > 0)
        
        tryCatch({
          yrobj <- yr2(truth = prcfiltered[["classification"]], scores = prcfiltered[selected_scores], high = rep(FALSE, length(selected_scores)))
          
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
          print("Main App plot success") #DEBUG
          output$Main_ErrorText <- renderText("")
          
        }, error = function(e) {
          plot_data(NULL)
          output$Main_ErrorText <- renderText("Not enough data")
        })
        
        output$Main_PRCPlot <- renderPlot({
          plot_info <- plot_data()
          print("Main App rendering") # DEBUG
          print(plot_info)
          if (!is.null(plot_info)) {
            tryCatch({
              draw.prc(plot_info$yrobj, lty = plot_info$lty_styles, col = plot_info$col_styles, lwd = 2, balanced = TRUE, main = paste0(plot_info$gene_s, " PRCs for ", paste(plot_info$selected_scores, collapse = ", ")))
              abline(h = 90, lty = "dashed")
              legend("left", legend = c(paste("# of Pathogenic and Likely Pathogenic:", plot_info$P_org), paste("# of Benign and Likely Benign:", plot_info$B_org)), pch = 15, bty = "n")
            }, error = function(e) {
              showModal(modalDialog(
                title = 'Error',
                'Not enough data 2',
                easyClose = TRUE,
                footer = NULL
              ))
            })
          }
        }, width = 600, height = 600, res = 72)
      })
      
      # Download logic for Main App
      output$Main_downloadPlotPNG <- downloadHandler(
        filename = function() {
          paste("PRC_plot_", input$gene, ".png", sep = "")
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
          paste("PRC_Report_", input$gene, ".pdf", sep = "")
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
      
    } else {
      # prcdata logic for Upload
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
        
        if (ncol(df) != 6) {
          showModal(modalDialog(
            title = "Error",
            "The uploaded dataset does not have the required number of columns. Ensure your file has exactly these columns: base__hugo, gnomad__af, varity_r__varity_r, alphamissense__am_pathogenicity, revel__score, classification.",
            easyClose = TRUE,
            footer = modalButton("Close")
          ))
          reset("file_full")
          return(NULL)
        }
        
        colnames(df) <- c("base__hugo", "gnomad__af", "varity_r__varity_r", 
                          "alphamissense__am_pathogenicity", "revel__score", "classification")
        return(df)
      })
      
      # Plotting logic for Upload
      observeEvent(input$plotButton, {
        df <- prcdata()
        gene_s <- input$gene
        exclude_common_variants <- input$common_variant_filter
        selected_scores <- input$scores
        print("Upload: begin plot") # DEBUG
        
        names(df)[names(df) == "varity_r__varity_r"] <- "VARITY"
        names(df)[names(df) == "alphamissense__am_pathogenicity"] <- "AlphaMissense"
        names(df)[names(df) == "revel__score"] <- "REVEL"
        names(df)[names(df) == "gnomad__af"] <- "gnomAD_AF"
        
        if (exclude_common_variants) {
          df <- df[is.na(df$gnomAD_AF) | df$gnomAD_AF <= 0.005, ]
        }
        
        df <- df[order(df$classification), ]
        
        prcfiltered <- df %>%
          filter(base__hugo == gene_s)
        print(prcfiltered) #DEBUG
        
        B_org <- sum(prcfiltered$classification == TRUE & rowSums(!is.na(prcfiltered[selected_scores])) > 0)
        P_org <- sum(prcfiltered$classification == FALSE & rowSums(!is.na(prcfiltered[selected_scores])) > 0)
        
        tryCatch({
          yrobj <- yr2(truth = prcfiltered[["classification"]], scores = prcfiltered[selected_scores], high = rep(FALSE, length(selected_scores)))
          
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
          print("Upload plot success") #DEBUG
          output$errorText <- renderText("")
          
        }, error = function(e) {
          plot_data(NULL)
          output$errorText <- renderText("Not enough data")
        })
        
        output$prcPlot <- renderPlot({
          plot_info <- plot_data()
          print("Upload rendering") # DEBUG
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
      
      # Download logic for Upload
      output$downloadPlotPNG <- downloadHandler(
        filename = function() {
          paste("PRC_plot_", input$gene, ".png", sep = "")
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
          paste("PRC_Report_", input$gene, ".pdf", sep = "")
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
