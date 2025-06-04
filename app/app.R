# app.R

library(shiny)

# ─────────────────────────────────────────────────────────────────────────────
# 1. TWO PURE R FUNCTIONS (shared by both tabs/tasks)
# ─────────────────────────────────────────────────────────────────────────────

# Function #1: transforms a numeric vector by applying multiplier and offset
transform_vector <- function(x, multiplier, offset) {
  return(x * multiplier + offset)
}

# Function #2: computes the mean of a numeric vector
compute_mean_scalar <- function(x) {
  if (length(x) == 0) {
    return(NA_real_)
  }
  return(mean(x, na.rm = TRUE))
}


# ─────────────────────────────────────────────────────────────────────────────
# 2. UI: fluidPage with a tabsetPanel containing two tabs (two tasks)
# ─────────────────────────────────────────────────────────────────────────────

ui <- fluidPage(
  titlePanel("App with Two Different Tasks"),

  # Tabbed interface: Task 1 and Task 2
  tabsetPanel(

    # ─────────────────────────────────────────────────────────────────────────
    # Task 1: CSV upload + transform one column
    # ─────────────────────────────────────────────────────────────────────────
    tabPanel(
      title = "Task 1: CSV → Transform & Mean",
      sidebarLayout(
        sidebarPanel(
          # (A) CSV upload
          fileInput(
            inputId  = "file1_csv",
            label    = "1) Upload a CSV file",
            accept   = c("text/csv", "text/comma-separated-values,text/plain", ".csv")
          ),
          tags$hr(),

          # (B) Two scalar inputs
          numericInput(
            inputId = "file1_multiplier",
            label   = "2) Multiplier (scalar):",
            value   = 1,
            step    = 0.1
          ),
          numericInput(
            inputId = "file1_offset",
            label   = "3) Offset (scalar):",
            value   = 0,
            step    = 0.1
          ),
          tags$hr(),

          # (C) Dynamic selector for numeric columns in the uploaded data
          uiOutput("file1_column_selector")
        ),

        mainPanel(
          # Preview first 10 rows
          h4("Data Preview (first 10 rows)"),
          tableOutput("file1_data_table"),
          tags$hr(),

          # Summary of the ORIGINAL column
          h4("Original Column Summary"),
          verbatimTextOutput("file1_original_summary"),
          tags$hr(),

          # Summary of the TRANSFORMED column
          h4("Transformed Column Summary"),
          verbatimTextOutput("file1_transformed_summary"),
          tags$hr(),

          # Display ONE scalar: the mean of the transformed column
          h4("Mean of Transformed Column (single scalar)"),
          textOutput("file1_mean_transformed")
        )
      )
    ),

    # ─────────────────────────────────────────────────────────────────────────
    # Task 2: Manual numeric-vector input → transform & median
    # ─────────────────────────────────────────────────────────────────────────
    tabPanel(
      title = "Task 2: Manual Vector → Transform & Median",
      sidebarLayout(
        sidebarPanel(
          # (D) TextArea for comma-separated numbers
          textAreaInput(
            inputId = "manual_vector",
            label   = "1) Enter a comma-separated numeric vector:",
            value   = "1, 2, 3, 4, 5",
            rows    = 3,
            placeholder = "e.g. 10, 20, 30, 40"
          ),
          tags$hr(),

          # (E) Re‐use the same scalar inputs for Task 2
          numericInput(
            inputId = "file2_multiplier",
            label   = "2) Multiplier (scalar):",
            value   = 1,
            step    = 0.1
          ),
          numericInput(
            inputId = "file2_offset",
            label   = "3) Offset (scalar):",
            value   = 0,
            step    = 0.1
          )
        ),

        mainPanel(
          # Show the parsed “original” vector
          h4("Original Vector (parsed)"),
          verbatimTextOutput("file2_original_vector"),
          tags$hr(),

          # Show transformed vector
          h4("Transformed Vector"),
          verbatimTextOutput("file2_transformed_vector"),
          tags$hr(),

          # Compute & display the MEDIAN of the transformed vector (single scalar)
          h4("Median of Transformed Vector (single scalar)"),
          textOutput("file2_median_transformed")
        )
      )
    )
  )
)


# ─────────────────────────────────────────────────────────────────────────────
# 3. SERVER logic
# ─────────────────────────────────────────────────────────────────────────────

server <- function(input, output, session) {

  # ─────────────────────────────────────────────────────────────────────────
  # —— TASK 1: handle CSV upload, column transform, and mean
  # ─────────────────────────────────────────────────────────────────────────

  # 1A. Reactive: read the uploaded CSV for Task 1
  file1_data <- reactive({
    req(input$file1_csv)
    read.csv(input$file1_csv$datapath, header = TRUE, stringsAsFactors = FALSE)
  })

  # 1B. Generate a selectInput of numeric columns for Task 1
  output$file1_column_selector <- renderUI({
    df <- file1_data()
    numeric_cols <- names(df)[sapply(df, is.numeric)]
    if (length(numeric_cols) == 0) {
      helpText("No numeric columns found in the uploaded file.")
    } else {
      selectInput(
        inputId  = "file1_chosen_col",
        label    = "4) Select a numeric column to transform:",
        choices  = numeric_cols,
        selected = numeric_cols[1]
      )
    }
  })

  # 1C. Show first 10 rows of the uploaded data
  output$file1_data_table <- renderTable({
    head(file1_data(), n = 10)
  }, striped = TRUE, hover = TRUE)

  # 1D. Original column summary
  output$file1_original_summary <- renderPrint({
    req(input$file1_chosen_col)
    df <- file1_data()
    summary(df[[input$file1_chosen_col]])
  })

  # 1E. Reactive: apply transform_vector() to the chosen column
  file1_transformed_col <- reactive({
    req(input$file1_chosen_col)
    df <- file1_data()
    original <- df[[input$file1_chosen_col]]
    m <- as.numeric(input$file1_multiplier)
    b <- as.numeric(input$file1_offset)
    transform_vector(original, m, b)
  })

  # 1F. Transformed column summary
  output$file1_transformed_summary <- renderPrint({
    summary(file1_transformed_col())
  })

  # 1G. Compute and output a single scalar (mean) for Task 1
  output$file1_mean_transformed <- renderText({
    mean_val <- compute_mean_scalar(file1_transformed_col())
    paste0("Mean = ", round(mean_val, digits = 3))
  })


  # ─────────────────────────────────────────────────────────────────────────
  # —— TASK 2: parse manual vector, transform, and compute median
  # ─────────────────────────────────────────────────────────────────────────

  # 2A. Reactive: parse the comma-separated input into a numeric vector
  file2_original_vec <- reactive({
    # Split on commas, remove whitespace, and convert to numeric
    req(input$manual_vector)
    str_vals <- strsplit(input$manual_vector, ",")[[1]]
    nums <- as.numeric(trimws(str_vals))
    # Remove any non-numeric or NA entries
    nums <- nums[!is.na(nums)]
    return(nums)
  })

  # 2B. Show the parsed “original” vector
  output$file2_original_vector <- renderPrint({
    vec <- file2_original_vec()
    if (length(vec) == 0) {
      "No valid numeric entries found."
    } else {
      vec
    }
  })

  # 2C. Apply transform_vector() to that manual vector
  file2_transformed_vec <- reactive({
    vec <- file2_original_vec()
    m   <- as.numeric(input$file2_multiplier)
    b   <- as.numeric(input$file2_offset)
    transform_vector(vec, m, b)
  })

  # 2D. Show the transformed vector
  output$file2_transformed_vector <- renderPrint({
    file2_transformed_vec()
  })

  # 2E. Compute and display a single scalar (MEDIAN) for Task 2
  output$file2_median_transformed <- renderText({
    transformed <- file2_transformed_vec()
    if (length(transformed) == 0) {
      "Median = NA (no data)"
    } else {
      med_val <- median(transformed, na.rm = TRUE)
      paste0("Median = ", round(med_val, digits = 3))
    }
  })
}

# ─────────────────────────────────────────────────────────────────────────────
# 4. RUN THE APP
# ─────────────────────────────────────────────────────────────────────────────

shinyApp(ui = ui, server = server)
