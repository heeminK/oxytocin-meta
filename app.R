library(shiny)
library(metafor)
library(DT)

# Define UI for the app
ui <- fluidPage(
  titlePanel("Forest Plot"),
  sidebarLayout(
    sidebarPanel(
      fileInput("file1", "Choose CSV File",
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv")),
      checkboxInput("use_title", "Use custom title", FALSE),
      textInput("title", "Title", "Forest Plot"),
      sliderInput("size", "Font Size",
                  min = 8, max = 24, value = 12, step = 1),
      sliderInput("ci_level", "Confidence Interval Level",
                  min = 0.9, max = 0.99, value = 0.95, step = 0.01),
      selectInput("method", "Method",
                  choices = c("DL", "REML", "ML"), selected = "DL"),
      checkboxInput("show_weights", "Show Weights", TRUE),
      checkboxInput("show_pooled", "Show Pooled Effect Size", TRUE)
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Data Table", dataTableOutput("data_table")),
        tabPanel("Forest Plot", plotOutput("forest_plot"))
      )
    )
  )
)

# Define server logic
server <- function(input, output) {
  
  # Read CSV file
  data <- reactive({
    data <- read.csv(input$file1$datapath)
    return(data)
  })
  
  # Create data table
  output$data_table <- renderDataTable({
    data()
  })
  
  # Create forest plot
  output$forest_plot <- renderPlot({
    # Meta-analysis model
    res <- rma(yi, vi, data = data(), method = input$method)
    # Forest plot
    forest(res, 
           slab = data()$study, 
           xlim = c(-5, 5),
           at = log(c(0.05, 0.25, 1, 4)),
           atransf = exp,
           ilab = cbind(data()$n.e, data()$n.c),
           ilab.xpos = c(-1.5, 1.5),
           main = ifelse(input$use_title, input$title, NULL),
           fontsize = input$size,
           mlab = "",
           level = input$ci_level,
           showweights = input$show_weights,
           showpooled = input$show_pooled)
  })
  
}

# Run the app
shinyApp(ui = ui, server = server)