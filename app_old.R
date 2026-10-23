# Load packages
library(shiny)
library(metafor)
library(DT)
library(dplyr)
library(shinyWidgets)
# library(devtools)
source("https://raw.githubusercontent.com/heeminK/oxytocin-meta/analysis/analysis_func.R")
#source("shiny_func.R")
# Define UI for the app
ui <- fluidPage(titlePanel("Forest Plot"),
                sidebarLayout(
                  sidebarPanel(
                    width = 3,
                    fileInput(
                      "file1",
                      "Choose CSV File",
                      accept = c("text/csv",
                                 "text/comma-separated-values,text/plain",
                                 ".csv")
                    ),
                    checkboxInput("use_title", "Use custom title", FALSE),
                    textInput("title", "Title", "Forest Plot"),
                    sliderInput(
                      "eq_level",
                      "Equivalence Bound Level",
                      min = 0.10,
                      max = 0.25,
                      value = 0.10,
                      step = 0.10
                    ),
                    awesomeRadio(
                      inputId = "mod_level",
                      label = "Moderator analysis",
                      choices = c(
                        "Between vs. within",
                        "Administration",
                        "Focal vs. Non-focal",
                        "Executive function"
                      )
                    )
                  ),
                  mainPanel(
                    width = 9,
                    tabsetPanel(
                      tabPanel("Data Table", dataTableOutput("data_table")),
                      tabPanel("Forest Plot", plotOutput("forest_plot")),
                      tabPanel("Funnel Plot", plotOutput("funnel_plot")),
                      tabPanel("Moderator", plotOutput("forest_moderator"))
                  ))
                ))

# Define server logic
server <- function(input, output) {
  # Read CSV file
  data <- reactive({
    if (is.null(input$file1)) {
      # Default CSV file
      data <-
        read.csv("https://raw.githubusercontent.com/heeminK/oxytocin-meta/main/ES.csv")
      
      # Data pre-processing
      data <- data |>
        filter(!is.na(yi))
      data$N_PLC[is.na(data$N_PLC)] <- 0
      data <- data |> arrange(-yi)
    } else {
      
      data <- read.csv(input$file1$datapath)
    }
    return(data)
  })
  
  # Create data table
  output$data_table <- renderDataTable({
    subset(data(), select = -c(X, STUDY_ID)) 
  })
  
  # Create forest plot
  output$forest_plot <- renderPlot({
    # Meta-analysis model
    res <- rma.mv(
      yi = yi,
      V = vi, 
      random = ~ 1 | STUDY_ID/X, 
      data = data(),
      slab = paste0(gsub( ", .*$", "", AUTHORS), " et al. (", YEAR, ")")
    )
    
    # Set the plot margins using par()
    par(mar = c(3.5, 1, 1, 1))
    
    # Forest plot
    forest_CIs <- with(
      data(),
      metafor::forest(
        yi,
        vi,
        xlim = c(-7, 4),
        ylim = c(-3, nrow(data()) + 3),
        efac = 1,
        slab = paste0(gsub(", .*$", "", AUTHORS), " et al. (", YEAR, ")"),
        ilab = EXP_OUTCOME_MEASURE,
        ilab.xpos = -4,
        at = c(-2,-1, 0, 1, 2),
        xlab = "Hedges' g",
        order = yi,
        cex = 1.15,
        header = "Author(s) (Year)",
        fontsize = input$size,
      )
    )
    op <- par(cex=1.15, font=2)
    text(-4, res$k+2, "Outcome measure")
    par(op)
    par(new=TRUE)
    with(
      data(),
      metafor::forest(
        yi,
        vi,
        xlim = c(-7, 4),
        ylim = c(-3, nrow(data()) + 3),
        level = 90,
        lwd = 3,
        efac = NA,
        annotate = FALSE,
        refline = NA,
        slab = NA,
        xaxt = "n",
        xlab = "",
        lty = c("solid", "blank"),
        order = yi,
        fontsize = input$size
      )
    )
    rect(
      -1 * input$eq_level,
      0,
      input$eq_level,
      20,
      col = rgb(
        red = 0.2,
        blue = 0.2,
        green = 0.2,
        alpha = 0.3
      ),
      border = "transparent"
    )
    
    rect(
      -1 * input$eq_level,
      -1.5,
      input$eq_level,
      -0.5,
      col = rgb(
        red = 0.2,
        blue = 0.2,
        green = 0.2,
        alpha = 0.3
      ),
      border = "transparent"
    )
    
    abline(h = 0)
    
    segments(
      res$ci.lb,
      -1,
      res$ci.ub,
      -1,
      lty = 1,
      col = "red",
      lwd = 1.5
    )
    segments(-0.028,
             -1,
             0.227,
             -1,
             lty = 1,
             col = "red",
             lwd = 3)
    
    points(
      x = res$b,
      y = -1,
      col = "red",
      pch = 23,
      bg = "white",
      cex = 1.5
    )
    points(
      x = res$ci.lb,
      y = -1,
      col = "red",
      pch = "|",
      cex = 1.25
    )
    points(
      x = res$ci.ub,
      y = -1,
      col = "red",
      pch = "|",
      cex = 1.25
    )
    text(4,-1, pos = 2, cex = 1.15, est_lab(res))
    text(-7,-1, pos = 4, cex = 1.15, mlabfun("Multilevel RE Model", res))
    
  })
  
  # Create contour-enhanced funnel plot
  output$funnel_plot <- renderPlot({
    res <- rma.mv(
      yi = yi,
      V = vi, 
      random = ~ 1 | STUDY_ID/X, 
      data = data(),
      slab = paste0(gsub( ", .*$", "", AUTHORS), " et al. (", YEAR, ")")
    )
    par(mar = c(5, 5, 1, 1))
    funnel.rma(
      res,
      level = c(90, 95, 99),
      back = "white",
      shade = c("light grey", "#5C6BC0", "#9FA8DA"),
      legend = FALSE,
      xlab = "Hedge's g",
      col = case_when(
        data()$MISSING_DATA == TRUE ~ 'red',
        data()$MISSING_DATA == FALSE ~ 'black'
      ),
      cex = 2,
      digits = 2L,
      cex.lab = 2,
      cex.axis = 2
    )
    legend(
      "topright",
      inset = 0,
      bg = "white",
      legend = c(
        expression("0.10 < p" <= "1.00"),
        expression("0.05 < p" <= "0.10"),
        expression("0.01 < p" <= "0.05"),
        expression("0.00 < p" <= "0.01")
      ),
      bty = "o",
      fill = c("light grey", "#5C6BC0", "#9FA8DA", "white"),
      box.col = "black",
      cex = 2
    )
  })
  
  # Forest plots for moderator analysis
  output$forest_moderator <- renderPlot({
    res <- rma.mv(
      yi = yi,
      V = vi, 
      random = ~ 1 | STUDY_ID/X, 
      data = data(),
      slab = paste0(gsub( ", .*$", "", AUTHORS), " et al. (", YEAR, ")")
    ) 
    
    par(mar = c(3.5, 1, 1, 1))
    if (input$mod_level == "Between vs. within") {
    # Between vs within
      generate_mod_fp(data(), "EXP_DESIGN", res)
    } else if (input$mod_level == "Administration") {
    # intranasal vs other
      generate_mod_fp(data(), "OT_ROUTE", res)
    } else if (input$mod_level == "Focal vs. Non-focal") {
    # primary vs secondary
      generate_mod_fp(data(), "EX_NONSOCIAL", res)
    } else if (input$mod_level == "Executive function") {
      # working memory vs flexibility vs inhibition
      generate_mod_fp(data(), "EXP_EF_DOMAIN", res)
      }
  })
}

# Run the app
shinyApp(ui = ui, server = server)
