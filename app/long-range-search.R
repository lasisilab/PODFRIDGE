# shiny-based Web app for long range familial search

# 1. Load required libraries
library(shiny)

# 2. Define UI for the application
ui <- fluidPage(
  titlePanel("Long Range Familial Search"),

  tabsetPanel(
    tabPanel(
      title = "About",
      tags$h3("About this Application"),
      tags$p("This application calculates the probability of identifying relatives in a database based on user-defined parameters."),
      tags$p("It uses a mathematical model to estimate the likelihood of finding a relative given the size of the population, the size of the database, and the degree of cousinship."),
      tags$p("For more information, please refer to the documentation or contact the developer.")
    ),

    tabPanel(
      title = "Identification Probability",
      tags$h3("Calculate the Probability of Identifying Relatives in a Database"),
      sidebarLayout(
        sidebarPanel(
          h4("Search Parameters"),
          numericInput(
            inputId = "population",
            label   = "Population Size (int):",
            value   = 3e8,
            min     = 0,
            step    = 1e7
          ),
          numericInput(
            inputId = "database",
            label   = "Database Size (int):",
            value   = 3e7,
            min     = 0,
            step    = 1e6
          ),
          numericInput(
            inputId = "degree",
            label  = "Degree of cousinship (int, view siblings as the 0-th cousins):",
            value  = 2,
            min    = 0,
            step   = 1
          ),
          numericInput(
            inputId = "relatives",
            label   = "Number of Relatives (int):",
            value   = 10,
            min     = 0,
            step    = 1
            ),
          tags$hr(),

          h4("IBD Segment Search Parameters"),
          numericInput(
            inputId = "min_num_seg",
            label   = "Minimum Number of IBD Segments (int):",
            value   = 2,
            min     = 0,
            step    = 1
          ),
          numericInput(
            inputId = "min_length_seg",
            label   = "Mininum Length of IBD Segments (cM):",
            value   = 6,
            min     = 0,
            step    = 1
          ),
        ),

        mainPanel(
          h4("Parameter Summary"),
          tableOutput("search_params"),
          tags$hr(),

          h4("The probability of identifying at least a relative in the database is:"),
          verbatimTextOutput("probability")
        )
      )
    ),

    tabPanel(
      title = "Simulations",
      tags$h3("Simulating generalized population")
    )
  )
)

# 3. Define custom functions

# calculate the probability of observing a match
genetic_match <- function(g, m = 6, num_chrs = 22, genome_size = 35, min_num_seg = 2) {
  m = m / 100 # Convert m from cM to fraction
  f = exp(-2 * g * m) / 2^(2 * g - 2) # calculate the probability of sharing a detectable IBD segment
  pr = 1 - pbinom(min_num_seg - 1, num_chrs + genome_size * 2 * g, f) # calculate probability
  return(pr)
}

# 4. Define server logic
server <- function(input, output, session) {

  # Reactive expression to calculate the probability
  probability <- reactive({
    N <- input$population
    K <- input$database
    g <- input$degree + 1
    R <- input$relatives

    if (N <= 0 || K <= 0 || g < 0 || R < 0) {
      return(NA)
    }

    # Calculate the probability of identifying a relative in the database
    p_match <- genetic_match(g = g, m = input$min_length_seg, min_num_seg = input$min_num_seg) * (K / N)
    p_identify <- 1 - (1 - p_match)^R
    return(p_identify)
  })

  # Output the search parameters
  output$search_params <- renderTable({
    df <- data.frame(
      Parameter = c("Population Size", "Database Size", "Degree of Cousinship", "Number of Relatives",
                    "Minimum Number of Segments", "Minimum Length of Segments (cM)"),
      Value = c(input$population, input$database, input$degree, input$relatives,
                input$min_num_seg, input$min_length_seg)
    )
  })

  # Output the calculated probability
  output$probability <- renderPrint({
    prob <- probability()
    if (is.na(prob)) {
      "Please enter valid parameters."
    } else {
      paste("Identification Probability:", round(prob, 4))
    }
  })
}

# 5. Run the application
shinyApp(ui = ui, server = server)

