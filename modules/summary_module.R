# Summary Dataset 模块

# UI 部分
mod_summary_ui <- function(id) {
  ns <- NS(id)  # 用于命名空间
  tabPanel(h4('Summary Dataset'), value = 'startsets',
           includeMarkdown(paste0(getwd(),"/www/title.Rmd") ),
           div(style = "height:2px;border-width:0;color:gray;background-color:gray", hr()),
           # Single Animated Chart
           column(6, plotlyOutput(ns("location_stat"), width = "auto", height = "300px") %>% withSpinner(type = 5)),
           column(3, plotlyOutput(ns("sampletype_stat"), width = "auto", height = "300px") %>% withSpinner(type = 5)),
           column(3, plotlyOutput(ns("month_stat"), width = "auto", height = "300px") %>% withSpinner(type = 5)),
           column(12, div(style = "height:2px;border-width:0;color:gray;background-color:gray", hr()),
                  div(style = "color:#337ab7;text-align:left;", titlePanel("All Sequenced samples:"))
           ),
           DT::dataTableOutput(ns('datatable')),
           br(),
           #includeHTML("./www/html/footer_v4.html")
           includeHTML(paste0(getwd(),"/www/html/footer_v4.html"))
  )
}

# Server 逻辑
mod_summary_server <- function(id, datafr2, summary_sub_datafr2) {
  moduleServer(id, function(input, output, session) {
    
    output$location_stat <- renderPlotly({
      plot_ly(summary_sub_datafr2, x = ~Value, y = ~reorder(Location, Value), frame = ~Month_Season, 
              type = "bar", orientation = "h", color = ~Mark, 
              colors = c("Month" = "skyblue", "Season" = "orange")) %>%
        layout(title = "Sample size comparison (Monthly & Seasonal)", 
               xaxis = list(title = "Sample size"), 
               yaxis = list(title = "City"), 
               showlegend = TRUE)
    })
    
    output$sampletype_stat <- renderPlotly({
      dat_sampletype <- as.data.frame(table(datafr2$Sample_type))
      plot_ly(dat_sampletype, labels = ~Var1, values = ~Freq, type = 'pie') %>%
        layout(title = 'Percentage of each sample type',xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
        yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
    })
    
    output$month_stat <- renderPlotly({
      dat_testday <- as.Date(datafr2$Test_day, format = "%d-%m-%Y")
      dat_m <- as.data.frame(table(month(dat_testday)))
      plot_ly(dat_m, x = ~Var1, y = ~Freq, type = 'bar') %>% 
        layout(title= 'Sample size of each month in 2024', yaxis = list(title = 'Number of samples'),xaxis = list(title = 'Months'), barmode = 'stack')
    })
    
    output$datatable <- DT::renderDataTable ({
      datatable(datafr2,selection = 'none',escape = FALSE,rownames = FALSE, filter = 'top', options = list(searchHighlight = TRUE,
                               columnDefs = list(list(width = '45%', targets =17),list(className = 'dt-center', targets =0:2)),pageLength = 10)) %>%
        formatStyle(1, valueColumns=1,color = '#337ab7', cursor = 'pointer')		
    })
    
    observeEvent(input$datatable_cell_clicked,{
      info = input$datatable_cell_clicked$value
      #print(info)
      if(length(datafr1[which(datafr1$Sample_id %in% info), ])!=0){
        # Specify the file path
        file_path <- paste0(getwd(),"/results/",datafr1[which(datafr1$Sample_id %in% info), "Chip_id"],"/06finalreports/",info,".docx")  
        #print(file_path)
        if(file.exists(file_path)) {
          # File exists; open it with microsoft word
          shell.exec(file_path)
          print("File exists and has been opened successfully.")
          
        } else { if(length(datafr1[which(datafr1$Sample_id %in% info), "Chip_id"])!=0){
          showModal(modalDialog(
            title = "Report not issued yet！",
            paste("Please copy this chip ID:",datafr1[which(datafr1$Sample_id %in% info), "Chip_id"],"to the TEST REPORT panel, and generate the testing report for the selected sample"),
            easyClose = TRUE,
            footer = NULL))
        }
          
        }
      }
      
    })
    
    
  })
}
