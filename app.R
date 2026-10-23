# Load packages
# (Shinylive detects the packages from these library() calls)
library(shiny)
library(metafor)
library(DT)
library(dplyr)
library(shinyWidgets)
library(stringr)

# Local files must sit in the same folder as app.R
source("analysis_func.R")

# Define UI for the app
ui <- fluidPage(
  titlePanel("Effects of oxytocin administration on non-social executive functions"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      radioGroupButtons(
        inputId = "eq_level",
        label = "Equivalence bound level",
        choices = c("0.10", "0.20", "0.25")
      ),
      radioGroupButtons(
        inputId = "contour",
        label = "Display for funnel plot",
        choices = c("Default", "Contour-enhanced")
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
      ),
      h5(
        "This app presents the data and results associated with a meta-analysis described in the following paper:",
        br(), br(),
        "Kang et al. (2022). Effects of exogenous oxytocin administration on non-social executive function in humans: A preregistered systematic review and meta-analysis protocol.",
        br(), br(),
        tags$a(href = "https://doi.org/10.31219/osf.io/8fzdy",
               "https://doi.org/10.31219/osf.io/8fzdy", target = "_blank")
      )
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        tabPanel("Data Table", DT::DTOutput("data_table")),
        tabPanel("Codebook", DT::DTOutput("codebook")),
        tabPanel("Forest Plot", plotOutput("forest_plot", height = "600px")),
        tabPanel("Funnel Plot", plotOutput("funnel_plot")),
        tabPanel("Moderator", plotOutput("forest_moderator", height = "600px"))
      )
    )
  )
)

# Define server logic
server <- function(input, output) {
  # Read CSV file
  data <- reactive({
    data <- read.csv("ES.csv")

    # Data pre-processing
    data <- data |>
      filter(!is.na(yi))
    data$N_PLC[is.na(data$N_PLC)] <- 0
    # Change the direction of effects if EXP_DIR_BENEFICIAL == "NEG"
    data <- data |>
      mutate(yi = ifelse(
        EXP_DIR_BENEFICIAL == "NEG",
        (-1) * yi,
        yi)
      )
    data <- data |> arrange(-yi)

    return(data)
  })

  # Create a data table for study characteristics
  output$data_table <- DT::renderDT({
    subset(data(), select = -c(X, STUDY_ID))
  })

  # Create a codebook data table
  output$codebook <- DT::renderDT({
    read.csv("codebook.csv")
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
    par(mar = c(3.5, 1, 0, 1))

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
        header = "Author(s) (Year)"
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
        order = yi
      )
    )
    rect(
      -1 * as.double(input$eq_level),
      0,
      as.double(input$eq_level),
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
      -1 * as.double(input$eq_level),
      -1.5,
      as.double(input$eq_level),
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
    par(mar = c(8.5, 5, 1.5, 1))

    if (input$contour == "Contour-enhanced") {

      funnel.rma(
        res,
        xlim = c(-1.5, 1.5),
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
      mtext(expression(paste("Black dots represent data that were directly obtained from the articles, \nand red dots represent data that were indirectly obtained using alternative methods.  \nThe vertical line represents the overall effect size (Hedges\u2019 g = 0.10).")),
            side = 1,
            line = 7,
            adj = 0,
            font = 1,
            cex = 1.15
      )
    } else {
      funnel.rma(
        res,
        xlim = c(-1.5, 1.5),
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
      mtext(expression(paste("Black dots represent data that were directly obtained from the articles, \nand red dots represent data that were indirectly obtained using alternative methods. \nThe vertical line represents the overall effect size (Hedges\u2019 g = 0.10).")),
            side = 1,
            line = 7,
            adj = 0,
            font = 1,
            cex = 1.15)
    }
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

    par(mar = c(3.5, 0, 0, 0))
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
