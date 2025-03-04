# analysis  模块

# analysis_module.R
# Optimized Shiny module for sequencing data analysis

# library(shiny)
# library(shinyBS)  # For popovers
# library(DT)       # For interactive tables
# library(promises)  # For async execution
# library(future)    # For non-blocking tasks
# library(shinyWidgets)  # For progress bars
# library(dplyr)

# UI Module
mod_analysis_ui <- function(id) {
  ns <- NS(id)  # Namespace for module
  tabPanel(
    title = h4('Analyze Dataset'), value = ns('analysis_fastq'),
    
    fluidRow(
      column(2, 
             div(style="display: inline-block;vertical-align:right; ",h3('Experimential data')),
             div(style="display: inline-block;vertical-align:right;", bsButton(ns("hint_sampleinput"), label = "", icon = icon("circle-info"),style="", size = "medium")),
             bsPopover(
               id = ns("hint_sampleinput"), title = "",
               content = "Click on a row to edit/delete, then modify and click the corresponding button.",
               placement = "right", trigger = "hover", options = list(container = "body")
             ),
             textInput(ns("chip_id"), "Flowcell ID"),
             textInput(ns("sample_id"), "Sample ID"),
             numericInput(ns("nucli_conce"), "Nucleic Acid Concentration (ng/μl)", 0),
             numericInput(ns("adaptor"), "Adaptor", 0),
             numericInput(ns("dna_conce"), "Library Concentration (ng/μl)", 0),
             textInput(ns("seq_id"), "Sequencing File"),
             actionButton(ns("add"), "Add"),
             actionButton(ns("edit"), "Edit"),
             actionButton(ns("deleteRows"), "Delete")
      ),
      column(10, 
             shiny::tags$br(),
             h4("NGS Info"),
             dataTableOutput(ns("TBL1"))
      )
    ),
    hr(),
    
    fluidRow(
      column(4, 
             actionButton(ns("exc"), "Confirm and Analyze", icon = icon('magnifying-glass'), 
                          style = "color: #fff; background-color: #337ab7;padding:60px 60px; 
                          font-size:24px;margin-top:50px;display:block;margin-left: auto;margin-right:auto;"),
             br(), br(),
             

      ),
      column(1,offset=0),
      
      column(5, 
             # Sub-tab panel inside the main panel
             tabsetPanel(
               id = ns("exc_analysis"),
               tabPanel(h3("Processing Log", style = "font-size:22px;"),
                        div(style="height:500px;width:1200px;overflow-y:auto;border:1px solid #ddd;padding:10px",verbatimTextOutput(ns("each_step"))),
                        shiny::tags$style(type='text/css','#each_step{color:black;font-size:15px;}')
                        #verbatimTextOutput(ns("log_file"))
               )
             )
      ),
      column(2,offset = 0)
    )
  ) # End of tabPanel "Analyze Dataset"
}



# Server Module

mod_analysis_server <- function(id, input_data, datafr1,datafr2) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    # ✅ If input_data is NULL, initialize an empty dataframe
    if (is.null(input_data) || nrow(input_data) == 0) {
      input_data <- data.frame(
        Chip_id = character(),
        Sample_id = character(),
        Extracted_NAC = numeric(),
        Adaptor = numeric(),
        Library_NAC = numeric(),
        Rawfastq_id = character(),
        stringsAsFactors = FALSE
      )
    }
    
    # ✅ Reactive Values to Store Data and Analysis State
    # ✅ Initialize reactiveValues with input_data or previous_data
    #rv <- reactiveValues(df = if (!is.null(input_data) && nrow(input_data) > 0) input_data else datafr1)
    rv <- reactiveValues(df = input_data, row_selected = NULL)
    datafr1_rv <- reactiveValues(df= datafr1)
    results_path <- reactiveValues(path = NULL)
    file_num <- reactiveValues(ndf = NULL)
    #pgress <- reactiveValues(each_step ="",fpp ="",log_text="")
    pgress <- reactiveValues(each_step ="")
    log_text <- reactiveVal("")
    folder_counts <- reactiveVal(NULL)
    # ✅ Add New Sample (Fixing NULL Input Issues)
    observeEvent(input$add, {
      req(input$sample_id, input$seq_id)  # Ensure inputs are not empty
      
      new_entry <- data.frame(
        Chip_id = input$chip_id,
        Sample_id = input$sample_id,
        Extracted_NAC = input$nucli_conce,
        Adaptor = input$adaptor,
        Library_NAC = input$dna_conce,
        Rawfastq_id = input$seq_id,
        stringsAsFactors = FALSE
      )
      
      # ✅ Ensure reactivity is triggered by reassigning `rv$df`
      rv$df <- rbind(rv$df, new_entry)
      showNotification("Sample added successfully!", type = "message")
    })
    
    
    # ✅ Delete Selected Row (Fixing NULL Selection Issues)
    observeEvent(input$deleteRows, {
      if (!is.null(input$TBL1_rows_selected) && length(input$TBL1_rows_selected) > 0) {
        rv$df <- rv$df[-as.numeric(input$TBL1_rows_selected), ]
        showNotification("Sample deleted!", type = "message")
      }
    })
    
    # ✅ Edit Selected Row (Fixing NULL Selection Issues)
    observeEvent(input$edit, {
      if (!is.null(input$TBL1_rows_selected) && length(input$TBL1_rows_selected) > 0) {
        cols_to_edit <- c('chip_id', 'sample_id', 'nucli_conce', 'adaptor', 'dna_conce', 'seq_id')
        colnms <- c('Chip_id', 'Sample_id', 'Extracted_NAC', 'Adaptor', 'Library_NAC', 'Rawfastq_id')
        rv$row_selected <- input$TBL1_rows_selected
        
        walk2(cols_to_edit, colnms, ~{
          rv$df[input$TBL1_rows_selected, ..2] <<- input[[..1]]
        })
        showNotification("Sample updated successfully!", type = "message")
      }
    })
    
    # ✅ Render Table (Fixing Empty Data Issues)
    output$TBL1 <- renderDataTable({
      req(rv$df)  # Ensure `rv$df` exists
      datatable(rv$df, selection = 'single', options = list(dom = 'frtip', pageLength = 10)) %>%
        formatStyle(0, fontWeight = 'bold')
    })
    
    
    
    # ✅ Confirm and Start Analysis (Fixing NULL Input Issues)
    observeEvent(input$exc, {
      if (nrow(rv$df) == 0) {
        showNotification("No data available for analysis.", type = "error")
        
        return()
      }

      # Merge new data into existing data (Update sequencinginfo.csv.csv)
      #  print(rbind(datafr1, as.data.frame(newlines)))
      # updated_data <- unique(rbind(datafr1, as.data.frame(newlines)))
      # str(datafr1_rv$df)
      # str(rv$df)
      # print(setdiff(colnames(datafr1_rv$df),colnames(rv$df)))
      updated_data <- unique(rbind(datafr1_rv$df, rv$df))

      write.csv(updated_data, seqinfo_file, row.names = FALSE)  # ✅ Save updated table

      showNotification("Table seqinfo updated successfully!", type = "message")


      # Merge new data into existing data (Update clinicalinfo.csv)
      sampleclinical_df <- data.frame("Sample_id" = unique(rv$df$Sample_id),"Name"= "","Gender"="",
                                      "Age"="","Sampling_day"="","Test_day"="","Tel"="","Infections"="","Requesting_physician"="",
                                      "Sample_type"="","Sample_status"="","Requesting_apartment"="","Requesting_hospital"="","Location"="",
                                      "Symptoms"="","Concerning_pathogens"="","Pretreatment"="","Pathogens"="")
      print(sampleclinical_df)
      clin_update_data <- unique(rbind(datafr2,sampleclinical_df))

      write.csv(clin_update_data, clininfo_file, row.names = FALSE)  # ✅ Save updated table

      showNotification("Table clininfo also updated successfully!", type = "message")
      
      newlines <- isolate(data.frame(rv$df))
      file_num$ndf <- unique(newlines)
      print(file_num$ndf)
      results_path$path <- file.path(results_folders, unique(newlines$Sample_id))
      print(results_path$path)
      str(results_path$path)
      
      # ✅ Ensure Directories Exist
      # Create directories (loop through each path)
      for (p in results_path$path) {
        dir.create(p, recursive = TRUE, showWarnings = FALSE)
        fastq_dest <- file.path(p, "fq")
        dir.create(fastq_dest, recursive = TRUE, showWarnings = FALSE) # 
      }
      
      # ✅ Copy FASTQ Files
      for (i in seq_len(nrow(file_num$ndf))) {
        if (file.exists(paste0(rawfastq_file,file_num$ndf$Rawfastq_id[i]))) {
          file.copy(paste0(rawfastq_file,file_num$ndf$Rawfastq_id[i]), file.path(results_path$path[i], "fq"))
        } else {
          showNotification(paste0("File missing: ", file_num$ndf$Rawfastq_id[i]), type = "error")
        }
      }
      
      # ✅ Run Analysis Asynchronously
      isolate({
        pgress$each_step <- paste0("\nStarting metagenomic analysis pipeline at ", Sys.time(),"\n")
      })
      for (p in results_path$path) {
        log_text("") # Clear Log before new run 
        future({
          run_pngsanalysis(p, file.path(p, "fq"), 72)
        }) %...>% (function(result) {
          isolate({
            pgress$each_step <- paste0(pgress$each_step, result)
            output$each_step <- renderText({pgress$each_step})
          })
        })
      }
    })
    
    # Render UI reactively
    #output$each_step <- renderText({pgress$each_step})
    output$each_step <- renderText({
      invalidateLater(1000) # Re-run this observe every 1 second
      if(input$exc >0 ){
        newlines <- isolate(data.frame(rv$df))
        file_num$ndf <- unique(newlines)
        if(dim(file_num$ndf)[1]>0){
          paste0(dim(file_num$ndf)[1]," will be process after sequenced",pgress$each_step)
        } else paste0("Waiting for data analysis!")
      } else paste0("Waiting for data analysis!")
    })
    

  })
}
