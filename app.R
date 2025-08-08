# Gene Visualization: Winter Lab Database of bulk RNA-seq Data 1/22/2025
# Tyler Therron, MS
# Deborah Winter, PhD
# Version 4.3.1.9

#  ======================================. Libraries ==================

library(shiny)
library(dplyr)
library(tidyr)
library(ggplot2) # Not used!
library(plotly)
library(readr)
library(shiny.fluent)
library(htmlwidgets) 
library(webshot) 
library(future)
library(DT)
library(future.apply)
library(RColorBrewer)
library(logger)
library(ragg)

options(
  shiny.maxRequestSize = 1000 * 1024^2,
  future.globals.maxSize = 12 * 1024^3,  # Allow up to 8 GiB for globals
  warn = -1,
  shiny.reactlog = TRUE
)

#  ======================================. Global Section for help functions and data ==================

source("./Sourced_Functions/PreDefined_Functions.R")
source("./Sourced_Functions/HM_PreDefined_Functions.R")
source("./Sourced_Functions/HM_Plotting_Functions.R")
source("./Sourced_Functions/Heatmap_Integration_Function.R")

#  ======================================. Define UI ==================
# Configure logger
log_file_path <- "/dev/null"
# log_file_path <- "/Users/ttm3567/Documents/January2025/GeneExpressionVisualization_v4.3.1.8_app.log"
# log_file_path <- "/srv/shiny-server/internal_server_GEV_app.log"

# if (!dir.exists(dirname(log_file_path))) {
#   dir.create(dirname(log_file_path), recursive = TRUE)
# }
# log_appender(appender_file(log_file_path))
# log_threshold(DEBUG)
# log_layout(layout_glue_generator(format = '{time} {level} [{namespace}] {fn}(): {msg}'))

# ui <- fluidPage(
#   # Include custom styles in the head section
#   tags$head(
#     tags$style(HTML("
#             .grid-container {
#               display: grid;
#               grid-template-columns: repeat(2, 1fr);
#               grid-gap: 20px;
#               padding: 20px;
#             }
#             .grid-item {
#               background-color: #f8f9fa;  /* Light grey background */
#               padding: 10px;
#               border-radius: 5px;
#               box-shadow: 0 2px 4px rgba(0,0,0,0.1);  /* Subtle shadow for depth */
#               overflow: hidden; /* Ensures the content fits in the grid area */
#             }
#         "))
#   ),
#   titlePanel("Winter Lab Database of Bulk RNA-seq"),
#   tabsetPanel(
#     id = "main_tabs",
#     tabPanel(
#       "Introduction",
#       h2("Welcome to the Winter Lab Database of Bulk RNA-seq!"),
#       br(),
#       p(
#         style = "font-size: 20px; margin-left: 40px; margin-right: 40px;",
#         "This application contains published and unpublished Bulk RNA-seq datasets from mice.
#       The primary purpose of this tool is to explore normalized gene expression levels among experimental groups from preloaded or user-loaded Bulk RNA-seq datasets."
#       ),
#       br(),
#       br(),
#       fluidRow(tags$div(
#         style = "display: flex; justify-content: center; margin-left: 40px; margin-right: 40px;",
#         DT::dataTableOutput("Dataset_Info")
#       )),
#       br(),
#       br(),
#       fluidRow(
#         p(
#           style = "font-size: 20px; margin-left: 40px; margin-right: 40px;",
#           "For instructions on how to properly format datatables and create groupfiles, click the link below:"
#         ),
#         tags$a(
#           style = "font-size: 20px; margin-left: 40px; margin-right: 40px;",
#           id = "newapp",
#           href = "https://winterlab-webapps.fsm.northwestern.edu/Template_GroupFileMaker_v3.1.2_App_20240521/",
#           "Template Data and Experimental Group Formatting App."
#         ),
#         tags$script(
#           HTML(
#             "$(document).on('click', '#newapp', function (e) {
#                       e.preventDefault();
#                       window.open($(this).attr('href'));
#                   });"
#           )
#         ),
#         br(),
#         br(),
#         br(),
#         p(
#           style = "font-size: 20px; margin-left: 40px; margin-right: 40px;",
#           "Below is a table with definitions of all the acronyms used in this web application."
#         ),
#         tags$div(
#           style = "display: flex; justify-content: left; margin-left: 40px; margin-right: 40px;",
#           tableOutput("legend_dataset")
#         ),
#         br(),
#         br()
#       )
#     ),
#     tabPanel(
#       "Data Visualization",
#       br(),
#       sidebarLayout(
#         sidebarPanel(
#           div(
#             p(
#               "To visualize gene expression in mouse datasets.\n Select one or two gene(s) of interest, and \nplot an interactive boxplot for gene expression data"
#             )
#           ),
#           
#           # checkboxInput("select_all", "Select All Preloaded Data"),
#           actionButton("select_all_button", "Select All"),
#           actionButton("deselect_all_button", "Deselect All"), 
#           br(),
#           checkboxGroupInput("gene_expression_data", "Preloaded Data:", choices = NULL),
#           checkboxInput("uploaded_data", "Uploaded Data:"),
#           div(
#             p(
#               "To download a plot as a PNG image, click the camera icon in the menu at the top-right of the rendered plot.",
#               style = "font-style: italic; font-size: smaller;"
#             )
#           ),
#           conditionalPanel(
#             condition = "input.uploaded_data == true",
#             fileInput(
#               'uploaded_data_file1', 'Choose CSV File',
#               accept = c(
#                 'text/csv',
#                 'text/comma-separated-values,text/plain',
#                 '.csv'
#               )
#             ),
#             fileInput(
#               "groups_File1", "Upload sample groups",
#               multiple = FALSE,
#               accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")
#             )
#           ),
#           conditionalPanel(
#             condition = "input.uploaded_data == true",
#             tags$a(id = "newapp", href = "https://winterlab-webapps.fsm.northwestern.edu/Template_GroupFileMaker_v3.1.2_App_20240521/", "Create Formatted Experimental Group File."),
#             tags$script(
#               HTML(
#                 "$(document).on('click', '#newapp', function (e) {
#                       e.preventDefault();
#                       window.open($(this).attr('href'));
#                   });"
#               )
#             )
#           ),
#           radioButtons("cpm_or_fpkm", "CPM or FPKM preloaded data repository:", c("CPM" = "cpm", "FPKM" = "fpkm"), selected = "cpm"),
#           radioButtons("colSel_pri", "Select Gene Nomenclature:", c("Ensembl ID" = "one", "Gene Symbol" = "two"), selected = "two"),
#           radioButtons("gene_select_method", "Upload a gene list or select individual genes to visualize:", c("Individual Genes" = "gene", "Gene List CSV" = "csv"), selected = "gene"),
#           conditionalPanel(
#             condition = "input.gene_select_method == 'csv'",
#           fileInput("gene_file", "Upload CSV of Genes", 
#                     accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv"))),
#           conditionalPanel(
#             condition = "input.gene_select_method == 'gene'",
#           uiOutput('select_pri')),
#           width = 3
#         ),
#         mainPanel(br(),
#                   tags$div(
#                     style = "overflow-y: scroll; height: 1200px;",
#           uiOutput("dataset_tabs")), br(),
#           div(
#             p(
#               "(The boxplots are in a scrollable window.)",
#               style = "font-style: italic; font-size: smaller;"
#             )
#           )
#         )
#       )
#     ),
#     tabPanel( "Gene Expression Heatmaps",
#               # This used to be in a sidebarLayout with the selectInputs and conditionalPanels in a sidebarPanel with width 3. The mainPanel had width 9. 
#               br(),
#               selectInput("ExpGrp_or_Sample", "Visualize by Experimental Groups or Samples?", choices = c("Group", "Sample"), selected = "Group"),
#               selectInput("cluster_hm_method", "Select Method to Cluster Heatmap Genes by:", choices = c("No Clustering" = "none",
#                                                                                                          "K-Means" = "kmeans",
#                                                                                                          "Hierarchical" = "hierarchical"), selected = "none"),
#               conditionalPanel( # CPM files - this is where the groupmaker app link will be
#                 condition = "input.cluster_hm_method == 'kmeans'",
#                 numericInput("number_o_clusters", "Number of K-Means Clusters:", value = 2, min = 0, step = 1)),
#               conditionalPanel( # CPM files - this is where the groupmaker app link will be
#                 condition = "input.cluster_hm_method == 'hierarchical'",
#                 selectInput("dendrogram", "Dendrogram?", choices = c("No", "Yes"), selected = "No")),
#               
#               tags$div(
#                 style = "overflow-y: scroll; height: 1200px;",
#               uiOutput("heatmap_tabs")), br(),
#               div(
#                 p(
#                   "(The heatmap is in a scrollable window.)",
#                   style = "font-style: italic; font-size: smaller;"
#                 )
#               ))
#   )
# ) # UI script

ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      .grid-container {
        display: grid;
        grid-template-columns: repeat(2, 1fr);
        grid-gap: 20px;
        padding: 20px;
      }
      .grid-item {
        background-color: #f8f9fa;
        padding: 10px;
        border-radius: 5px;
        box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        overflow: hidden;
      }
    "))
  ),
  
  titlePanel("Winter Lab Database of Bulk RNA-seq"),
  
  sidebarLayout(
    # 👉 Shared Sidebar
    sidebarPanel(
      div(
        p("To visualize gene expression in mouse datasets...") # shortened for brevity
      ),
      actionButton("select_all_button", "Select All"),
      actionButton("deselect_all_button", "Deselect All"),
      checkboxGroupInput("gene_expression_data", "Preloaded Data:", choices = NULL),
      checkboxInput("uploaded_data", "Uploaded Data:"),
      conditionalPanel(
        condition = "input.uploaded_data == true",
        fileInput('uploaded_data_file1', 'Choose CSV File'),
        fileInput("groups_File1", "Upload sample groups")
      ),
      conditionalPanel(
        condition = "input.uploaded_data == true",
        tags$a(id = "newapp", href = "https://winterlab-webapps.fsm.northwestern.edu/Template_GroupFileMaker_v3.1.2_App_20240521/", "Create Formatted Experimental Group File."),
        tags$script(HTML("
          $(document).on('click', '#newapp', function (e) {
              e.preventDefault();
              window.open($(this).attr('href'));
          });
        "))
      ),
      radioButtons("cpm_or_fpkm", "CPM or FPKM:", c("CPM" = "cpm", "FPKM" = "fpkm")),
      radioButtons("colSel_pri", "Gene Nomenclature:", choices = c("Ensembl ID" = "one", "Gene Symbol" = "two"), selected = "two"), # <- this sets "Gene Symbol" as default
      radioButtons("gene_select_method", "Select Genes:", c("Individual Genes" = "gene", "Gene List CSV" = "csv")),
      conditionalPanel("input.gene_select_method == 'csv'",
                       fileInput("gene_file", "Upload CSV of Genes")),
      conditionalPanel("input.gene_select_method == 'gene'",
                       uiOutput('select_pri')),
      # submitButton("plot"),
      width = 3
    ),
    
    # 👇 Dynamic Main Panel based on selected tab
    mainPanel(
      tabsetPanel(id = "main_tabs",
                  tabPanel("Introduction"),
                  tabPanel("Data Visualization"),
                  tabPanel("Gene Expression Heatmaps")
      ),
      conditionalPanel(
        condition = "input.main_tabs == 'Introduction'",
        h2("Welcome to the Winter Lab Database of Bulk RNA-seq!"),
        br(),
        p(style = "font-size: 20px; margin-left: 40px; margin-right: 40px;",
          "This application contains published and unpublished Bulk RNA-seq datasets from mice..."),
        br(), br(),
        fluidRow(tags$div(style = "display: flex; justify-content: center; margin-left: 40px; margin-right: 40px;",
                          DT::dataTableOutput("Dataset_Info"))),
        br(), br(),
        fluidRow(
          p(style = "font-size: 20px; margin-left: 40px; margin-right: 40px;",
            "For instructions..."),
          tags$a(style = "font-size: 20px; margin-left: 40px; margin-right: 40px;", id = "newapp", 
                 href = "https://winterlab-webapps.fsm.northwestern.edu/Template_GroupFileMaker_v3.1.2_App_20240521/",
                 "Template Data and Experimental Group Formatting App."),
          br(), br(), br(),
          p(style = "font-size: 20px; margin-left: 40px; margin-right: 40px;",
            "Below is a table with definitions of all the acronyms..."),
          tags$div(style = "display: flex; justify-content: left; margin-left: 40px; margin-right: 40px;",
                   tableOutput("legend_dataset")),
          br(), br()
        )
      ),
      conditionalPanel(
        condition = "input.main_tabs == 'Data Visualization'",
        tags$div(style = "overflow-y: scroll; height: 1200px;", uiOutput("dataset_tabs")),
        br(),
        div(p("(The boxplots are in a scrollable window.)", style = "font-style: italic; font-size: smaller;"))
      ),
      conditionalPanel(
        condition = "input.main_tabs == 'Gene Expression Heatmaps'",
        br(),
        
        # This used to be in a sidebarLayout with the selectInputs and conditionalPanels in a sidebarPanel with width 3. The mainPanel with everything else had width 9. 
        selectInput("ExpGrp_or_Sample", "Visualize by Groups or Samples?", choices = c("Group", "Sample")),
        selectInput("cluster_hm_method", "Clustering Method:", choices = c("No Clustering" = "none", "K-Means" = "kmeans", "Hierarchical" = "hierarchical")),
        conditionalPanel("input.cluster_hm_method == 'kmeans'",
                         numericInput("number_o_clusters", "K-Means Clusters:", value = 2, min = 0, step = 1)),
        conditionalPanel("input.cluster_hm_method == 'hierarchical'",
                         selectInput("dendrogram", "Dendrogram?", choices = c("No", "Yes"))),
        
        tags$div(style = "overflow-y: scroll; height: 1200px;", uiOutput("heatmap_tabs")),
        br(),
        div(p("(The heatmap is in a scrollable window.)", style = "font-style: italic; font-size: smaller;"))
      )
    )
  )
)
# ======= SERVER ====================================================

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  
  # log_info("Session started: {session$token}")
  # log_info("R version: {R.Version()$version.string}")
  # log_info("ggdendro version: {packageVersion('ggdendro')}")
  # log_info("ggplot2 version: {packageVersion('ggplot2')}")
  # log_info("dendextend version: {packageVersion('dendextend')}")

  selected_genes <- reactiveValues(genes = NULL)
  previous_selected_genes <- reactiveVal(character()) 
  
  # Reactive value to track previous dataset selection
  last_selected_datasets <- reactiveVal(NULL)
  
  RenPlo <- reactiveValues(data = list())  # Stores already processed datasets
  non_processed_df_list <- reactiveValues(names = character())  # To be processed
  Removed_DF <- reactiveValues(names = character())
  
  gene_file_genes <- reactive({
    req(input$gene_file)
    
    # Read the uploaded file and extract the gene symbols
    gene_data <- read.csv(input$gene_file$datapath, header = TRUE, stringsAsFactors = FALSE)
    
    # Assuming the gene symbols are in the first column
    gene_vector <- gene_data[, 1]  # Ensure it is a character vector and lowercase
    
    # Return the vector of gene symbols
    return(gene_vector)
  })
  
  processed_datasets <- reactiveValues(data = list(), groups = list())
  
  # take 3
  combined_data <- reactive({
    req(length(selected_genes$genes) > 0)
    
    # Preloaded dataset logic
    datasets <- setdiff(input$gene_expression_data, "User Uploaded Data")
    preloaded_result <- load_datasets(
      datasets,
      input$cpm_or_fpkm,
      selected_genes$genes,
      selNum_pri()
    )
    
    # Only pull user-uploaded data from the reactive functions
    # (NOT from processed_datasets$data), so it don't cause a feedback loop.
    user_uploaded_result <- list(data = list(), groups = list())
    if (input$uploaded_data) {
      user_data <- uploaded_expr_data()       # direct from reactive
      user_groups <- uploaded_group_data()    # direct from reactive
      
      if (!is.null(user_data)) {
        filter_col <- if (selNum_pri() == 1) "Symbol" else "Gene_Symbol"
        user_data <- user_data[user_data[[filter_col]] %in% selected_genes$genes, , drop = FALSE]
        user_uploaded_result$data[["User Uploaded Data"]] <- user_data
      }
      if (!is.null(user_groups)) {
        user_uploaded_result$groups[["User Uploaded Data"]] <- user_groups
      }
    }
    
    # Merge them
    result = list(
      data   = c(preloaded_result$data,  user_uploaded_result$data),
      groups = c(preloaded_result$groups, user_uploaded_result$groups)
    )
    return(result)
  })
  
  # Observe when combined_data changes and update processed_datasets
  observe({
    new_data <- combined_data()
    req(new_data)  # Make sure it's not null

    processed_datasets$data <- new_data$data
    processed_datasets$groups <- new_data$groups
  })

  # take 3 user upload data
  # Reactive for expression data
  uploaded_expr_data <- reactive({
    req(input$uploaded_data_file1)
    read.csv(input$uploaded_data_file1$datapath, check.names = FALSE)
  })
  
  # Reactive for group data
  uploaded_group_data <- reactive({
    req(input$groups_File1, uploaded_expr_data())
    
    expr_data <- uploaded_expr_data()
    group_df <- read.csv(input$groups_File1$datapath)
    
    # Validate that the sample column matches the expression data columns
    validation <- validate_uploaded_data(expr_data, group_df)
    
    # If validation fails, show a notification and stop
    if (!validation$valid) {
      showNotification(validation$message, type = "error")
      return(NULL)
    }
    
    transpose_groupfile(group_df)
  })
  
  # Once both are non-null, store them in processed_datasets
  observe({
    # Ensure both reactivities are valid
    req(uploaded_expr_data(), uploaded_group_data())
    
    processed_datasets$data[["User Uploaded Data"]]   <- uploaded_expr_data()
    processed_datasets$groups[["User Uploaded Data"]] <- uploaded_group_data()
    
    print("User data & group files both read, stored in processed_datasets.")
  })
  
  datasetSelected <- reactiveVal(NULL)
  
  lastDataset <- reactiveValues(data = NULL)
  
  observe({
    lastDataset$data <- input$gene_expression_data
  })
  
  observeEvent(input$cpm_or_fpkm, {
    print("input$cpm_or_fpkm")
    datasets <- datasets_info("./FPKM_expression_files", "./Primary")
    
    if (input$cpm_or_fpkm == "cpm") {
      cpm_path <- "./Primary"
      datasetSelected("cpm")
      cpm_choice_list <- c(sort_numeric(list.files(path = cpm_path, full.names = FALSE, recursive = FALSE)))
      updateCheckboxGroupInput(session, "gene_expression_data", choices = cpm_choice_list, selected = if (is.null(lastDataset$data) || length(lastDataset$data) == 0) { cpm_choice_list[1] } else { lastDataset$data })
    } else if (input$cpm_or_fpkm == "fpkm") {
      fpkm_path <- "./FPKM_expression_files"
      datasetSelected("fpkm")
      fpkm_choice_list <- c(sort_numeric(list.files(path = fpkm_path, full.names = FALSE, recursive = FALSE)))
      updateCheckboxGroupInput(session, "gene_expression_data", choices = fpkm_choice_list, selected = if (is.null(lastDataset$data) || length(lastDataset$data) == 0) { fpkm_choice_list[1] } else { lastDataset$data })
    }
  })
  
  # II
  observeEvent(input$select_all_button, {
    print("input$select_all_button")
    req(input$cpm_or_fpkm)
    cpm_path <- "./Primary"
    fpkm_path <- "./FPKM_expression_files"
    cpm_choice_list <- c(sort_numeric(list.files(path = cpm_path, full.names = FALSE, recursive = FALSE)))
    fpkm_choice_list <- c(sort_numeric(list.files(path = fpkm_path, full.names = FALSE, recursive = FALSE)))
    
    choices <- if (input$cpm_or_fpkm == "cpm") cpm_choice_list else fpkm_choice_list
    
    updateCheckboxGroupInput(session, "gene_expression_data", selected = choices)
  })
  
  # Deselect all datasets when clicking "Deselect All"
  # observeEvent(input$deselect_all_button, {
  #   updateCheckboxGroupInput(session, "gene_expression_data", selected = character(0))
  # })
  observeEvent(input$deselect_all_button, {
    print("input$deselect_all_button")
    # Update the UI to clear selections
    updateCheckboxGroupInput(session, "gene_expression_data", selected = character(0))
    
    isolate({
      # Move all datasets from RenPlo to Removed_DF
      removed_datasets <- names(RenPlo$data)
      print("Removed datasets:")
      print(removed_datasets)
      
      if (length(removed_datasets) > 0) {
        print("length(removed_datasets) > 0")
        Removed_DF$names <- unique(c(Removed_DF$names, removed_datasets))  # Track removed datasets
        RenPlo$data <- list()  # Clear all plots
        
        print("Deselect All clicked: Removed datasets moved to Removed_DF")
      }
    })
  })
  
  #II
  observeEvent(input$gene_expression_data, {
    req(input$cpm_or_fpkm)
    print("input$gene_expression_data")
    
    cpm_path <- "./Primary"
    fpkm_path <- "./FPKM_expression_files"
    
    cpm_choice_list <- c(sort_numeric(list.files(path = cpm_path, full.names = FALSE, recursive = FALSE)))
    fpkm_choice_list <- c(sort_numeric(list.files(path = fpkm_path, full.names = FALSE, recursive = FALSE)))
    
    selected_datasets <- input$gene_expression_data
    
    # Determine available choices based on CPM/FPKM selection
    choices <- if (input$cpm_or_fpkm == "cpm") cpm_choice_list else fpkm_choice_list
    
    # ✅ Only update "Select All" checkbox if the user actually selected/deselected all manually
    if (length(selected_datasets) == length(choices)) {
      updateCheckboxInput(session, "select_all", value = TRUE)
    } else if (length(selected_datasets) == 0) {
      updateCheckboxInput(session, "select_all", value = FALSE)
    } 
    # ❌ No automatic change to "select_all" when just unchecking one dataset
  })
  
  observe({
    datasets <- datasets_info("./FPKM_expression_files", "./Primary")
    if (!is.null(lastDataset$data) && length(lastDataset$data) > 0 && any(lastDataset$data %in% datasets$unavailable)) {
      updateRadioButtons(session, "cpm_or_fpkm", selected = "cpm", choices = c(CPM = "cpm"))
    } else {
      updateRadioButtons(session, "cpm_or_fpkm", selected = input$cpm_or_fpkm, choices = c(CPM = "cpm", FPKM = "fpkm"))
    }
  })
  
  observeEvent(input$cpm_or_fpkm, {
    if (input$cpm_or_fpkm == "cpm") {
      datasetSelected("cpm")
    } else if (input$cpm_or_fpkm == "fpkm") {
      datasetSelected("fpkm")
    }
  })
  
  selNum_pri <- reactive(switch(input$colSel_pri, one = 1, two = 2))
  
  title_pri <- reactive({
    get_data_title(input$gene_expression_data, input$uploaded_data_title)
  })
  
  choiceList_pri <- throttle(reactive(gene_choices_for_plots[, selNum_pri()]), 2000)
  
  # II
  observeEvent(input$gene_expression_data, {
    req(input$cpm_or_fpkm)

    selected_datasets <- input$gene_expression_data
    current_selected_genes <- selected_genes$genes  # Get currently selected genes
    previous_genes <- previous_selected_genes()  # Get previously selected genes

    # Identify newly selected datasets that haven't been processed yet
    new_datasets <- setdiff(selected_datasets, names(RenPlo$data))

    # Identify datasets that have been removed
    removed_datasets <- setdiff(names(RenPlo$data), selected_datasets)

    # Identify datasets that need updating (existing datasets where new genes were added)
    datasets_needing_update <- names(RenPlo$data)[length(setdiff(current_selected_genes, previous_genes)) > 0]

    # Identify datasets that were previously removed but are now re-selected
    reselected_datasets <- intersect(selected_datasets, Removed_DF$names)

    # Merge the new datasets, reselected datasets, and those needing an update
    datasets_to_process <- unique(c(new_datasets, reselected_datasets, datasets_needing_update))

    # Update tracking lists
    non_processed_df_list$names <- unique(c(non_processed_df_list$names, datasets_to_process))
    Removed_DF$names <- setdiff(Removed_DF$names, reselected_datasets)  # Remove re-selected datasets from Removed_DF

    # Remove plots from RenPlo for deselected datasets
    for (dataset in removed_datasets) {
      RenPlo$data[[dataset]] <- NULL
    }

    # Update the previously selected genes tracker
    previous_selected_genes(current_selected_genes)
  })
  
  # III
  preprocessing_selected_genes_datasets <- reactive({
    
    print("inside preprocessing_selected_genes_datasets")
    print(selected_genes$genes)
    req(processed_datasets$data, selected_genes$genes)
    
    # Check if any previously processed datasets need updates due to newly selected genes
    affected_datasets <- names(RenPlo$data)
    datasets_with_new_genes <- Filter(function(dataset) {
      missing_genes <- setdiff(selected_genes$genes, names(RenPlo$data[[dataset]]))
      return(length(missing_genes) > 0)
    }, affected_datasets)
    
    # Merge new datasets + datasets with new genes
    datasets_to_process <- unique(c(non_processed_df_list$names, datasets_with_new_genes))
    
    datasets_to_process <- datasets_to_process[!is.na(datasets_to_process)]
    processed_datasets$data
    
    # If nothing to process, return empty list
    #if (length(datasets_to_process) == 0) {
    #  return(list())
    #}
    
    # These lines update the plots whether one is added or removed
    valid_datasets <- processed_datasets$data[datasets_to_process]
    # Process everything
    valid_datasets <- processed_datasets$data
    valid_datasets <- valid_datasets[!sapply(valid_datasets, is.null)]
    
    results_list <- future_lapply(valid_datasets, function(dataset) {
      future_lapply(selected_genes$genes, function(gene) {
        gene_data <- dataset[dataset[[selNum_pri()]] == gene, ]
        process_gene_data(gene_data, gene, selNum_pri)
      })
    })
    
    named_results_list <- assign_names_to_sublist(results_list, selected_genes$genes)
    
    return(named_results_list)
  })
  
  
  
  # II
  formattedData_pri <- reactive({
    print("Entered formattedData_pri")
    results_list <- preprocessing_selected_genes_datasets()
    print(results_list)

    # log_info("Line 506 results_list")
    # log_info("data in results_list Data: {paste(capture.output(print(results_list)), collapse = '\n')}")
    
    # Only format the datasets that were just processed
    #groups_to_process <- processed_datasets$groups[names(results_list)]
    groups_to_process <- processed_datasets$groups

    group_and_data <- future_mapply(function(dataset_genes_list, groupfile_list) {
      format_and_merge_UI_chosen_genes(dataset_genes_list, groupfile_list)
    }, results_list, groups_to_process, SIMPLIFY = FALSE)
    
    # log_info("Line 513 adding new genes by filtering the chosen dataset")
    # log_info("data in observe Event Formatted Pri Data: {paste(capture.output(print(group_and_data)), collapse = '\n')}")

    return(group_and_data)
  })
  
  cpm_fpkm_labeler_Yaxis <- reactive({
    if (input$cpm_or_fpkm == "cpm") {
      label <- "CPM"
    } else if (input$cpm_or_fpkm == "fpkm") {
      label <- "FPKM"
    }
    return(label)
  })
  
  observe({
    if (input$gene_select_method == "gene") {
      # log_debug("Updating selected_genes$genes with input$mygene_pri: ", paste(input$mygene_pri, collapse=", "))
      selected_genes$genes <- input$mygene_pri
    } else {
      # log_debug("Updating selected_genes$genes with gene_file_genes()")
      selected_genes$genes <- gene_file_genes()
    }

  })
  
  # GOOD - original gene choice
  # This part spazzes out when more than one input is selected at once
  output$select_pri = renderUI({
    print("Entered output$select_pri")
    choiceList_pri <- choiceList_pri()
    
    selectizeInput("mygene_pri",
                  label = h5("Select Gene of Interest to Visualize:"),
                  multiple = TRUE,
                  selected = if (is.null(input$mygene_pri) || length(input$mygene_pri) == 0) { "" }
                              else { selected_genes$genes },
                  choices = choiceList_pri,
                  options = list(create = TRUE, placeholder = 'Type here or choose from the list'))
    
    
    
  })
  
  # What genes are selected?
  observeEvent(input$mygene_pri, {
    print("input$mygene_pri changed to:")
    print(input$mygene_pri)
  })
  
  output$dataset_tabs <- renderUI({
    print("output$dataset_tabs")
    print("selected_genes$genes")
    print(selected_genes$genes)
    
    # There must be at least 1 selected gene to open a new dataset tab
    req(length(selected_genes$genes) > 0)
    # Only include datasets that are in RenPlo and not in Removed_DF
    active_datasets <- setdiff(names(RenPlo$data), Removed_DF$names)
    
    dataset_tabs <- lapply(active_datasets, function(dataset) {
      tabPanel(
        title = dataset,
        uiOutput(paste0("gene_plots_", gsub("[^a-zA-Z0-9]", "_", dataset)))
      )
    })
    
    do.call(tabsetPanel, dataset_tabs)
  })
  
  # II
  observeEvent(formattedData_pri(), {
    print("Are we getting here?")
    formatted_data <- formattedData_pri()

    # Require formatted data
    req(formatted_data)

    # Ensure we are processing ALL datasets that need rendering
    datasets_to_process <- names(formatted_data)  # Process ALL datasets with new genes

    print("Datasets to process:")
    print(datasets_to_process)

    withProgress(message = 'Rendering plots...', value = 0, {
      total_datasets <- length(datasets_to_process)
      dataset_count <- 0

      # ✅ Ensure UI outputs are always created for all genes in active datasets
      lapply(datasets_to_process, function(dataset_name) {
        output[[paste0("gene_plots_", gsub("[^a-zA-Z0-9]", "_", dataset_name))]] <- renderUI({
          gene_plots <- lapply(names(formatted_data[[dataset_name]]), function(gene_name) {
            plot_data <- formatted_data[[dataset_name]][[gene_name]]

            if (!("Count" %in% names(plot_data))) {
              print(paste("Skipping gene", gene_name, "in dataset", dataset_name, "- missing 'Count' column"))

              # Show a message instead of an empty tab with italics
              return(div(
                class = "missing-gene-message",
                p(style = "color: grey40; font-size: 16px; text-align: center;",
                  HTML(paste("<em>Gene, <strong style='color: blue;'>", gene_name, "</strong> , cannot be visualized in dataset because it is absent from the data.</em>"))
                )
              ))
            }

            plot_output_id <- paste0(
              "plotly_plot_",
              gsub("[^a-zA-Z0-9]", "_", dataset_name), "_",
              gsub("[^a-zA-Z0-9]", "_", gene_name)
            )

            div(
              plotlyOutput(outputId = plot_output_id, height = "800px", width = "100%"),
              class = "grid-item"
            )
          })

          gene_plots <- Filter(Negate(is.null), gene_plots)

          div(
            class = "grid-container",
            do.call(tagList, gene_plots)
          )
        })

        # ✅ Fix: Pass dataset_name explicitly
        dataset_count <<- dataset_count + 1
        incProgress(1 / total_datasets, detail = paste("Preparing plots for", dataset_name))
      })

      # ✅ Second loop: Render plots for ALL datasets
      # lapply(datasets_to_process, function(dataset_name) {
      for (dataset_name in datasets_to_process) {
        # ⚠️ Fix: Check if this dataset already has genes stored
        existing_genes <- names(RenPlo$data[[dataset_name]])
        new_genes <- setdiff(names(formatted_data[[dataset_name]]), existing_genes)

        if (length(new_genes) == 0) {
          return()  # Skip if no new genes
        }

        print(paste("New genes detected in dataset:", dataset_name, "->", paste(new_genes, collapse = ", ")))

        lapply(new_genes, function(gene_name) {
          plot_data <- formatted_data[[dataset_name]][[gene_name]]

          # log_info("Line 631 adding new genes by filtering the chosen dataset")
          # log_info("data in observe Event Formatted Pri Data: {paste(capture.output(print(plot_data)), collapse = '\n')}")

          plot_output_id <- paste0(
            "plotly_plot_",
            gsub("[^a-zA-Z0-9]", "_", dataset_name), "_",
            gsub("[^a-zA-Z0-9]", "_", gene_name)
          )

          output[[plot_output_id]] <- renderPlotly({
            plot_the_data(plot_data, dataset_name, cpm_fpkm_labeler_Yaxis(), gene_name)
          })

          # ✅ Store newly rendered genes inside RenPlo
          if (!is.null(plot_data)) {
            if (is.null(RenPlo$data[[dataset_name]])) {
              RenPlo$data[[dataset_name]] <- list()
            }
            RenPlo$data[[dataset_name]][[gene_name]] <- plot_data
          }
        })

        # ✅ Fix: Pass dataset_name explicitly
        dataset_count <<- dataset_count + 1
        incProgress(1 / total_datasets, detail = paste("Rendering plots for", dataset_name))
      }#)
    }) # End of withProgress

    # ✅ Remove processed datasets from `non_processed_df_list`
    non_processed_df_list$names <- setdiff(non_processed_df_list$names, datasets_to_process)

    print("Updated RenPlo data after rendering:")
    print(names(RenPlo$data))
  })
  
  output$Dataset_Info <- renderDataTable({
    datatable(
      dataset_info,
      rownames = FALSE,
      escape = FALSE,
      options = list(
        autoWidth = FALSE,
        columnDefs = list(list(targets = "_all", className = "dt-center")),
        pageLength = -1,
        dom = 't'
      )
    )
  }, server = FALSE)
  
  output$legend_dataset <- renderTable(legend_info, bordered = TRUE, striped = TRUE, rownames = FALSE)
  
  # HM tab reactive functions ---------------------
  
  ## REACTIVE heatmaps
  averaged_data_list <- reactive({
    # Create a progress bar
    withProgress(message = "Averaging data...", value = 0, {
      # Get the original data and groups
      
      GEtable_list <- processed_datasets$data
      grps_list <- processed_datasets$groups
      selected_genes <- selected_genes$genes
      
      total_datasets <- length(GEtable_list)
      progress_increment <- 1 / total_datasets  # Calculate the increment for each dataset
      
      # Initialize an empty list to store the results
      averaged_data <- list()
      
      # Loop over the datasets, updating the progress bar
      for (i in seq_along(GEtable_list)) {
        dataset_name <- names(GEtable_list)[i]
        GEtable <- GEtable_list[[i]]
        grps <- grps_list[[i]]
        
        tryCatch({
          # Process the data
          grps <- groupfile_reformat(grps)
          filtered_GEtable <- filter_genes(GEtable, selected_genes)
          GEtable_avg <- exp_group_avg_gene_expression_calculator(filtered_GEtable, grps)
          
          # Store the result
          averaged_data[[dataset_name]] <- GEtable_avg
        }, error = function(e) {
          # Handle the error (e.g., log it or show a message)
          showNotification(paste("Error processing", dataset_name, ":", e$message), type = "error")
        })
        
        # Update the progress bar
        incProgress(progress_increment, detail = paste("Processing", dataset_name))
      }
      
      return(averaged_data)
    })
  })
  

  combined_HM <- reactive({
    if (input$main_tabs != "Gene Expression Heatmaps") {
      return(NULL)  # Exit the reactive function if not on Heatmap tab
    }
    
    clustering_method <- input$cluster_hm_method
    cluster_number <- input$number_o_clusters
    dendro <- input$dendrogram
    GEtable_list <- processed_datasets$data
    grps_list <- processed_datasets$groups
    ExpGrp_or_Sample <- input$ExpGrp_or_Sample
    selected_genes <- selected_genes$genes
    
    # Use averaged data if "Group" is selected
    if (ExpGrp_or_Sample == "Group") {
      data_list <- averaged_data_list()
      if (is.null(data_list)) {
        # Data not yet averaged or "Group" not selected
        return(NULL)
      }
    } else {
      data_list <- GEtable_list
    }
    
    print(" grps_list from within the combined HM function")
    print(grps_list)
    # Use Map to apply the helper function over data_list and grps_list
    # browser()
    heatmap_list <- Map(process_heatmap_data, 
                        data_list, grps_list, names(data_list), 
                        MoreArgs = list(clustering_method = clustering_method, 
                                        cluster_number = cluster_number, 
                                        dendro = dendro, 
                                        selected_genes = selected_genes, 
                                        ExpGrp_or_Sample = ExpGrp_or_Sample))
    # this is where my ERROR TRACKING is left off 20250116
    # browser()
    return(heatmap_list)
  })
  
  
  # Generate the heatmap tabs dynamically using combined_HM()
  output$heatmap_tabs <- renderUI({
    
    # Initialize the list of datasets to include the preloaded datasets
    datasets <- input$gene_expression_data
    
    # Include the uploaded dataset if the checkbox is checked
    if (input$uploaded_data == TRUE) {
      datasets <- c(datasets, "User Uploaded Data")
    }
    
    # Create tabs for each dataset
    heatmap_tabs <- lapply(datasets, function(dataset) {
      tabPanel(
        title = dataset,
        plotlyOutput(outputId = paste0("heatmap_plot_", gsub("[^a-zA-Z0-9]", "_", dataset)), height = "600px")
      )
    })
    
    # Combine all the tabs into one tabset panel
    do.call(tabsetPanel, heatmap_tabs)
  })
  
  # CHUNK LOAD
  observeEvent(combined_HM(), {
    
    heatmap_list <- combined_HM()
    chunk_size <- 5  # Adjust this based on available memory
    dataset_names <- names(heatmap_list)
    
    print("dataset_names")
    print(dataset_names) # this is NULL
    
    # browser()
    # Split dataset names into chunks
    if (!is.null(dataset_names)) {
      print("dataset_names is not null!")
      chunks <- split(dataset_names, ceiling(seq_along(dataset_names) / chunk_size))
      
      for (chunk_index in seq_along(chunks)) {
        chunk <- chunks[[chunk_index]]
        
        # Process each chunk
        lapply(chunk, function(dataset_name) {
          result <- heatmap_list[[dataset_name]]
          plot_output_id <- paste0("heatmap_plot_", gsub("[^a-zA-Z0-9]", "_", dataset_name))
          
          output[[plot_output_id]] <- renderPlotly({
            req(result)
            
            if (result$error) {
              # Render error annotation
              plot_ly() %>%
                layout(
                  annotations = list(
                    list(
                      x = 0.5, y = 0.5,
                      text = result$message,
                      showarrow = FALSE,
                      xref = "paper", yref = "paper",
                      xanchor = "center", yanchor = "middle",
                      font = list(size = 16, color = "grey20", family = "Arial Italic")
                    )
                  )
                )
            } else {
              # Render heatmap
              if (!is.null(result$message)) {
                result$heatmap %>%
                  layout(
                    annotations = list(
                      list(
                        x = 0.5, y = -0.25,
                        text = result$message,
                        showarrow = FALSE,
                        xref = "paper", yref = "paper",
                        xanchor = "center", yanchor = "top",
                        font = list(size = 14, color = "grey40", family = "Arial Italic")
                      )
                    ),
                    margin = list(b = 200)
                  )
              } else {
                result$heatmap
              }
            }
          })
        })
      }
    }
  })
  
  session$onSessionEnded(function() {
    # log_info("Session ended: {session$token}")
  })
} # SERVER script bracket

# Run the application 
shinyApp(ui = ui, server = server)
