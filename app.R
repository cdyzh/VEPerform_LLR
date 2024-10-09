library(shiny)

# Load the UI and server
source("ui.R")
source("server.R")

# Run the application
shinyApp(ui = ui, server = server)
