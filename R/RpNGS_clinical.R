# V1

#' Use RpNGS through shiny interface

RpNGS <- function(){
  options(warn=-1)
  #  options(shiny.maxRequestSize=30*1024^2)
  library(sf)
  library(shiny)
  library(shinyBS) #bsButton
  library(shinycssloaders) #withSpinner
  library(shinyFiles)
  library("plotly") # plotlyOutput
  library(leaflet) #leafletOutput
  #library(feather)
  library(data.table)
  library(DT)
  library(lubridate)
  library(readr)
  library(tidyverse)
  library(shinyjs)
  library(shinydashboard)
  library(shinyWidgets)
  library(Gviz)
  library(Rsamtools)
  library(GenomicAlignments)
  library(officer)
  library(flextable)
  library(rvest) # For parsing HTML 10.4
  library(ShortRead)
  library(scales)
  library(ps)
  library(future)
  library(promises)
  
  # code loacte on server 
  #    source("/home/dell/dataanalysis/projects/RpNGS_lifegen/RpNGS/R/pngs_se_analysis_n1.R", local = TRUE)
  
  Input <- structure(list(Chip_id = c("FT10005323"), Sample_id = c("lifegen20240830"), Extracted_NAC = c(
    "2.74"), Adaptor = c("126"),library_NAC = c("58.6"), Rawfastq_id = c(
      "FT10005323_L10_126.fq.gz")), class = c(
        "spec_tbl_df", "tbl_df", "tbl",
        "data.frame"
      ), row.names = c(NA, -1L))
  
  
  ui <- shiny::navbarPage(h4("pNGS v0.1.0"),id= "rpngs",selected="startsets",
                          theme=shinythemes::shinytheme("flatly"),# inTabsetm
                          
                          tabPanel(h4('Help'), value='help'
                                   ,shiny::tags$style(
                                     '.navbar-nav {width: 90%;}
				.navbar-nav :first-child{float:right}'
                                   )
                                   ,column(1				
                                   )
                                   ,column(10
                                           ,div(style="color:#337ab7; text-align:center;", titlePanel(h2("pNGS Operation Details"))), br()
                                           # Add custom CSS
                                           ,shiny::tags$style(HTML("
                                               .markdown-container {
                                              width: 100%;
                                              height: 900px;
                                              overflow: auto;
                                              border: 1px solid #ccc;
                                            }
                                               ")),
                                           shiny::tags$div(class = "markdown-container",
                                                           includeMarkdown("RpNGS/www/html/help.Rmd")  # Replace with your markdown file
                                           )
                                           
                                           #,includeMarkdown("RpNGS/www/html/help.Rmd")
                                           #,shiny::tags$iframe(style="height:900px; width:100%", src="https://sust365-my.sharepoint.cn/:b:/g/personal/4439_sust_cn/EU3FZwzUSTtBggbUckliO5UBMM8o9XMmSzeqqbAtPvAdRQ?e=S9dHoQ")
                                   )
                                   ,column(1
                                           
                                   )
                                   
                          ),
                          tabPanel(h4('Summary dataset'), value='startsets'
                                   ,includeMarkdown("RpNGS/www/title.Rmd")
                                   ,div(style="height:2px;border-width:0;color:gray;background-color:gray", hr())
                                   ,column(4, offset=0, style='padding:0px 20px 0px 0px;'
                                           ,shiny::tags$br()
                                           ,leafletOutput("location_stat", width = "auto", height = "300px") %>% withSpinner(type=5), br()
                                   )
                                   ,column(4, offset=0, style='padding:10px 0px 0px 10px;'
                                           ,br()
                                           , plotlyOutput("sampletype_stat", width = "auto", height = "300px") %>% withSpinner(type=5), br()
                                   )
                                   ,column(4, offset=0, style='padding:10px 20px 0px 0px;'
                                           ,shiny::tags$br()
                                           ,plotlyOutput("month_stat", width = "auto", height = "300px") %>% withSpinner(type=5), br()
                                   )
                                   ,column(12, offset=0, style='padding:0px;'
                                           ,div(style="height:2px;border-width:0;color:gray;background-color:gray", hr())
                                           ,div(style="color:#337ab7;text-align:left;", titlePanel("ALL Sequenced samples:"))
                                   )
                                   ,DT::dataTableOutput('datatable'), br()
                                   ,includeHTML("RpNGS/www/html/footer_v4.html")
                                   
                          ),
                          tabPanel(h4('Analyze dataset'), value='out_seqinfo'
                                   ,column(2, offset = 0
                                           ,style='padding:0px;'
                                           #,tags$br()
                                           ,div(style="display: inline-block;vertical-align:right; ",h3('Experimential data'))
                                           ,div(style="display: inline-block;vertical-align:right;", bsButton("hint_sampleinput", label = "", icon = icon("circle-info"),style="", size = "medium"))
                                           ,bsPopover(id = "hint_sampleinput", title = "",
                                                      content = paste0("To modify or delete a record in the table, first click on the row within the table, then fill in the updated information here, and click the Modify or Delete button."),
                                                      placement = "right", trigger = "hover", options = list(container = "body")
                                           )
                                           # inputs
                                           ,textInput("chip_id", "Flowcell ID")
                                           ,textInput("sample_id", "Sample ID")
                                           ,numericInput("nucli_conce", "Nucleic acid concentration(ng/μl)",0)
                                           ,numericInput("adaptor", "Adaptor",0)
                                           ,numericInput("dna_conce", "library concentration (ng/μl)",0)
                                           ,textInput("seq_id","Sequencing file")
                                           ,actionButton("add", "INPUT")
                                           ,actionButton("edit", "EDIT")
                                           ,actionButton("deleteRows", "DELETE")
                                           #,width = 2
                                   )
                                   ,column(10, offset = 0
                                           #,style='padding:0px;'
                                           ,shiny::tags$br()
                                           ,h4("NGS info")
                                           ,dataTableOutput("TBL1")
                                           
                                   )
                                   ,column(12, offset=0, style='padding:0px;'
                                           ,div(style="height:2px;border-width:0;color:gray;background-color:gray", hr())
                                           ,br(),br()
                                   )
                                   
                                   ,column(4, offset = 0
                                           ,style='padding:0px;'
                                           #,tags$br()
                                           ,actionButton("exc", "Confirm and Analyze",icon('magnifying-glass'),style="color: #fff; box-shadow: 1px 1px 1px 1px #888; background-color: #337ab7; border-color: #2e6da4"),br()
                                           
                                           ,br()
                                           ,tabsetPanel( id="exc_analysi",
                                                         tabPanel("Processing log",value="log",br()
                                                                  # ,verbatimTextOutput("each_step")
                                                                  #,actionButton("show_log","Show log file", style="color: #fff; background-color: #337ab7; border-color: #2e6da4"), br()
                                                                  ,div(style="width:800px;", verbatimTextOutput("filenum"))
                                                                  ,shiny::tags$style(type='text/css', '#filenum {color: black; font-size:15px;}')
                                                                  ,br()
                                                                  ,div(style="width:800px;", verbatimTextOutput("each_step"))
                                                                  ,shiny::tags$style(type='text/css', '#each_step {color: black; font-size:15px;}')
                                                         )
                                                         
                                           )	
                                   )
                                   ,column(1,offset = 0
                                   )
                                   ,column(5
                                           ,p(strong("File processing progress"), class = "text-center")
                                           ,h5("[0]：Transfer Data")
                                           ,progressBar(id = "pb0", value = 0, display_pct = T,status = "danger")
                                           ,h5("[1]：Quality Control")
                                           ,progressBar(id = "pb1", value = 0, display_pct = T,status = "primary")
                                           #,tags$style(".progress-bar-custom {background-color: #25c484;}")
                                           ,h5("[2]：Host reads removal")
                                           ,progressBar(id = "pb2", value = 0, display_pct = T,status = "info")
                                           ,h5("[3]：Classification")
                                           ,progressBar(id = "pb3", value = 0, display_pct = T,status = "warning")
                                           ,h5("[4]：Statical analysis")
                                           ,progressBar(id = "pb4", value = 0, display_pct = T,status = "success")
                                           
                                   )
                                   
                                   ,column(2, offset = 0
                                           
                                   )
                                   
                          ),
                          tabPanel(h4('Test report'), value='clinicalreport'
                                   ,column(2, offset=0.5, 
                                           wellPanel(
                                             div(style="display: inline-block;vertical-align:right;", textInput('reprot_chipid',h4('Flowcell ID:'),width='300px' ))
                                             ,div(style="display: inline-block;vertical-align:right;", actionButton("chazhao", "Search",style="vertical-align:right;color: #fff; box-shadow: 1px 1px 1px 1px #888; background-color: #800000; border-color: #800000"))  
                                             ,br()
                                             ,br()
                                             ,h4("Report status:")
                                             ,br()
                                             ,div(DTOutput("sampleid_list"),style = "height:600px; overflow-y: scroll;")
                                             #,div(id = 'myDiv',DTOutput("sampleid_list"),style = "height:500px; overflow-y: scroll;overflow-x: scroll;")
                                             
                                           )
                                   )
                                   ,column(10, offset=0.5,
                                           #wellPanel(
                                           conditionalPanel("input.sampleid_list_rows_selected>0",
                                                            div(style="display: inline-block;vertical-align:right; ",h3("Clinical info:"))
                                                            ,div(style="display: inline-block;vertical-align:right;", bsButton("search_gene_bs", label = "", icon = icon("info-circle"),style="", size = "small"))
                                                            ,bsPopover(id = "search_gene_bs", title = "",
                                                                       content = paste0("Double-click to fill in, press CTRL+Enter to confirm entry."),
                                                                       placement = "right", trigger = "hover", options = list(container = "body")
                                                            )
                                                            ,DTOutput("sampleclinical_td")
                                                            ,br()
                                                            
                                                            ,br()
                                                            # #-------------------
                                                            # ,dashboardBody(
                                                            #   box(title = "Data Path", status = "primary",height = "595" ,solidHeader = T,
                                                            #       plotOutput("trace_plot")),
                                                            #   box( title = "Case Analyses Details", status = "primary", height =
                                                            #          "595",width = "6",solidHeader = T,
                                                            #        column(width = 12,
                                                            #               DT::dataTableOutput("trace_table"),style = "height:500px; overflow-y: scroll;overflow-x: scroll;"
                                                            #        )
                                                            #   ))
                                                            # #-------------------
                                                            
                                                            # #-------------------
                                                            
                                                            ,h3("Pathogen info:")
                                                            #,splitLayout(
                                                            ,column(width = 3,offset=-0.5
                                                                    ,wellPanel(
                                                                      h3("Pathogens:"),
                                                                      div(DT::dataTableOutput("pathogen_selected"),style = "height:500px; overflow-y: scroll;"))
                                                            )
                                                            ,column(width = 9,offset=-0.5
                                                                    ,wellPanel(
                                                                      h3("Raw microbes list:")
                                                                      ,div(DT::dataTableOutput("pathogen_td") %>% withSpinner(type=5),style = "height:500px; overflow-y: scroll;overflow-x: scroll;"))
                                                            )
                                                            
                                                            ,column(12, offset=0, style='padding:0px;'
                                                                    ,div(style="height:2px;border-width:0;color:gray;background-color:gray", hr())
                                                                    ,div(style="color:#337ab7;", titlePanel(""))
                                                            )
                                                            ,downloadButton("download_doc", "Report generation",style="color: #fff; box-shadow: 1px 1px 1px 1px #888; background-color: #337ab7; border-color: #2e6da4")
                                                            ,br()
                                                            
                                                            
                                                            
                                                            
                                                            
                                           )
                                           
                                   )
                                   ,conditionalPanel("input.pathogen_td_rows_selected !=0 && input.sampleid_list_rows_selected >0",
                                                     column(width = 2,offset=-0.5,
                                                            h4("Alignment analysis"),
                                                            verbatimTextOutput("statsOutput")  # Display the calculated statistics
                                                     )                   
                                                     
                                                     
                                                     ,column(width = 10,offset=-0.5
                                                             
                                                             # ,wellPanel(
                                                             ,h4("Alignment view:")
                                                             ,plotOutput("alignmentPlot") %>% withSpinner(type=5)
                                                             #,div(DT::dataTableOutput("pathogen_td") %>% withSpinner(type=5),style = "height:500px; overflow-y: scroll;overflow-x: scroll;")
                                                             #)
                                                     )
                                                     
                                   )
                                   
                                   
                                   
                          )
                          
                          
                          
                          
  )
  
  
  
  server <- function(input, output, session){
    
    ####======== main data for pNGS datatable and search in section 1 ========####
    # load data for summary
    #datafr2 <- as.data.frame(fread(paste0("RpNGS/data/RpNGS_testdatatable.csv")),stringsAsFactors=F)
    
    datafr1 <- as.data.frame(fread(paste0("RpNGS/data/response/RpNGS_sequencinginfo.csv")),stringsAsFactors=F)
    datafr2 <- as.data.frame(fread(paste0("RpNGS/data/response/RpNGS_clinicalinfo.csv")),stringsAsFactors=F) #FOR SUMMARY SAMPLE NUMBER, HAVE TO FILL A CITY OF EACH SAMPLE
    sub_datafr2 <- datafr2[which(!is.na(datafr2$Test_day)),c("Location","Test_day")]
    
    # Add month and season columns
    sub_datafr2 <- sub_datafr2 %>%
      mutate(
        month = format(Test_day, "%m"), # Extract month
        year =paste0("20",format(Test_day,"%y")),
        season = case_when(        # Define seasons
          month %in% c("12", "01", "02") ~ "Winter",
          month %in% c("03", "04", "05") ~ "Spring",
          month %in% c("06", "07", "08") ~ "Summer",
          month %in% c("09", "10", "11") ~ "Fall",
          TRUE ~ NA_character_
        )
      )
    # Extract the year
    sub_datafr2 <- sub_datafr2[which(sub_datafr2$year %in% str_extract(date(), "\\d{4}") ),] 
    
    
    # Summarize counts for each month and season
    summary_sub_datafr2 <- sub_datafr2 %>%
      mutate(month_col = paste0("M", month)) %>% # Add a month column for pivoting
      pivot_longer(cols = c(month_col, season), names_to = "category", values_to = "value") %>%
      group_by(Location, value) %>%
      summarise(count = n(), .groups = "drop") %>%
      pivot_wider(names_from = value, values_from = count, values_fill = 0)
    
    #define custom function to add columns to data frame if they do not exist
    add_cols <- function(df, cols) {
      add <- cols[!cols %in% names(df)]
      if(length(add) !=0 ) df[add] <- 0
      return(df)
    }
    
    #add season columns if they don't already exist
    summary_sub_datafr2 <- add_cols(summary_sub_datafr2, c('Winter', 'Spring', 'Summer','Fall','M01','M02','M03','M04','M05','M06','M07','M08','M09','M10','M11','M12'))
    
    
    # Arrange the columns in the desired order
    result <- summary_sub_datafr2 %>%
      select(Location, starts_with("M"), Winter, Spring, Summer, Fall)
    
    #result$name2 <- paste0(result$Location,"市")
    
    
    ## Calculate the test amount of months and seasons
    if(format(Sys.Date(), "%m") == "1"){
      
    }
    # Add the number of orders for each district and county in Shaanxi for each month (unit: pieces).
    shaanxi_data <- tibble::tribble(
      ~adcode, ~name2, 
      610100, "Xi'an",
      610200, "Tongchuan",
      610300, "Baoji",
      610400, "Xianyang",
      610500, "Weinan",
      610600, "Yanan",
      610700, "Hanzhong",
      610800, "Yulin",
      610900, "Ankang",
      611000, "Shangluo",
    )
    
    shaanxi_data <- merge(shaanxi_data,result,by.x="name2",by.y="Location")
    
    shaanxi_data <- shaanxi_data[,c("adcode",paste0("M",format(Sys.Date(),"%m")),"Winter", "Spring", "Summer", "Fall")]
    colnames(shaanxi_data)<- c("adcode","M","Winter", "Spring", "Summer", "Fall")
    
    
    # load map data 
    c_map <- read_sf('RpNGS/image/Shaanxi.json')
    # sf convert to dataframe
    c_map <- as.data.frame(c_map)
    # Merge the observation data and map data according to administrative division codes.
    shaanxi <- merge(c_map, shaanxi_data, by = "adcode")
    # dataframe convert to sf
    shaanxi <- sf::st_as_sf(shaanxi, sf_column_name = "geometry")
    # Set a segmented discrete color palette
    # pal <- colorBin("Spectral", bins = pretty(shaanxi$pop), reverse = TRUE)
    pal <- colorNumeric("Spectral", domain = NULL)
    ####======== GENERAL REVIEW: location_status ========####
    
    output$location_stat <- renderLeaflet({
      # Plot the data on the map.
      leaflet(shaanxi) |>
        addTiles(
          # From leafletCN::amap()
          urlTemplate = "http://webrd02.is.autonavi.com/appmaptile?lang=zh_cn&size=1&scale=1&style=8&x={x}&y={y}&z={z}",
          options = tileOptions(tileSize = 256, minZoom = 3, maxZoom = 17),
          attribution = "&copy; <a href=\"http://amap.com\">amp.com</a >"
        ) |>
        addPolygons(
          stroke = F, # Not display the boundaries of each district and county
          weight = 1, # Set the boundary line width
          fillOpacity = 0.7, # Set the opacity of the polygon fill.
          fillColor = ~ pal(M),
          label = lapply(paste0(
            "Location：", "<b>", shaanxi$name, "</b>", "<br/>",
            "Monthly testing volume：", shaanxi$M, "<br/>",
            "Testing volume in Winter：", shaanxi$Winter,"<br/>",
            "Testing volume in Spring:", shaanxi$Spring,"<br/>",
            "Testing volume in Summer：", shaanxi$Summer,"<br/>",
            "Testing volume in Fall：", shaanxi$Fall
          ), htmltools::HTML),
          
        ) |>
        addLegend(
          position = "bottomright", title = "Test volume",
          pal = pal, values = ~M, opacity = 1.0,
          labFormat = labelFormat(
            suffix = "个",
            transform = function(x) round(x, digits = 2)
          )
        ) |>
        addScaleBar(position = "bottomleft")
      
    })
    ####======== GENERAL REVIEW: sampletype_stat ========####
    
    output$sampletype_stat <- renderPlotly({
      dat_sampletype <- as.data.frame(table(datafr2$Sample_type))
      plot_ly(dat_sampletype, labels = ~Var1, values = ~Freq, type = 'pie',
              textposition = 'inside',
              textinfo = 'label+text',
              insidetextfont = list(color = '#FFFFFF'),
              hoverinfo = 'percent',
              #text = ~paste(Var1),
              showlegend = FALSE) %>%
        layout(title = 'Percentage of each sample type',
               xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
               yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
    })
    
    ####======== GENERAL REVIEW: sampletype_stat ========####
    output$month_stat <- renderPlotly({
      dat_testday <- as.Date(datafr2$Test_day, format = "%d-%m-%Y")
      dat_m <- as.data.frame(table(month(dat_testday)))
      
      plot_ly(dat_m, x = ~Var1, y = ~Freq, type = 'bar') %>%
        
        layout(title= 'Sample size of each month in 2024', yaxis = list(title = 'Number of samples'),xaxis = list(title = 'Months'), barmode = 'stack')
      
    })
    
    
    
    ####======== Section 4: main data table ========####
    output$datatable <- DT::renderDataTable ({
      datatable(datafr2,
                selection = 'none',
                escape = FALSE,
                rownames = FALSE, filter = 'top', 
                options = list(searchHighlight = TRUE,
                               columnDefs = list(
                                 list(width = '45%', targets =17),
                                 list(className = 'dt-center', targets =0:2
                                 )
                               ),
                               pageLength = 20
                )) %>%formatStyle(1, valueColumns=1,color = '#337ab7', cursor = 'pointer')		
    })
    
    observeEvent(input$datatable_cell_clicked,{
      info = input$datatable_cell_clicked$value
      print(info)
      if(length(datafr1[which(datafr1$Sample_id %in% info), ])!=0){
        # Specify the file path
        file_path <- paste0("D:/OneDrive - 陕西科技大学/dataanalysis/RpNGS/RpNGS/data/results/",datafr1[which(datafr1$Sample_id %in% info), "Chip_id"],"/06finalreports/",info,".docx")  # Replace with the actual file path
        print(1)
        if(file.exists(file_path)) {
          # File exists; open it with microsoft word
          shell.exec(file_path)
          print("File exists and has been opened successfully.")
          
        } else { if(length(datafr1[which(datafr1$Sample_id %in% info), "Chip_id"])!=0){
          showModal(modalDialog(
            title = "Report not issued yet！",
            paste("Please generate the testing report for this sample on the clinical report page and retrieve the chip ID:",datafr1[which(datafr1$Sample_id %in% info), "Chip_id"]),
            easyClose = TRUE,
            footer = NULL))
        }
          
        }
      }
      
    })
    
    ################################################ 2. 测序数据处理 ################################################
    # 
    # source("./RpNGS/R/RpNGS_seqdataprocessing.R", local = TRUE)
    
    
    #-------------------------------2.1 录入数据--------------------------------
    # Init with some example data
    #data <- reactiveVal(Input)
    #Input <- as.data.frame(fread(paste0("RpNGS/data/batch_RpNGS_sequencinginfo.csv")),stringsAsFactors=F)
    
    
    #Input <- list(fread(paste0("RpNGS/data/RpNGS_sequencinginfo.csv")),stringsAsFactors=F)
    rv <- reactiveValues(df = Input, row_selected = NULL) 
    
    observeEvent(input$add,{
      # start with current data
      if(input$sample_id!="" && !is.null(input$seq_id) && input$add>0){
        rv$df <- rv$df %>%
          add_row(
            Chip_id = isolate(input$chip_id),
            Sample_id = isolate(input$sample_id),
            Extracted_NAC = isolate(input$nucli_conce),
            Adaptor = isolate(input$adaptor),
            library_NAC = isolate(input$dna_conce),
            Rawfastq_id = isolate(input$seq_id)
            
          )#  %>%
        # update data value
        #data()
      }
      
      
    }
    )
    observeEvent(input$deleteRows,{
      
      if (!is.null(input$TBL1_rows_selected)) {
        #data(data()[-as.numeric(input$TBL1_rows_selected),])
        rv$df <- rv$df[-as.numeric(input$TBL1_rows_selected), ]
      }
    })
    
    observeEvent(input$edit,{
      
      if (!is.null(input$TBL1_rows_selected)) {
        cols_to_edit <- c('chip_id', 'sample_id', 'nucli_conce', 'adaptor', 'dna_conce','seq_id')
        colnms <- c('Chip_id', 'Sample_id', 'Extracted_NAC', 'Adaptor', 'library_NAC','Rawfastq_id')
        "remember the row selected"
        rv$row_selected <- input$TBL1_rows_selected
        
        walk2(cols_to_edit, colnms, ~{rv$df[input$TBL1_rows_selected, ..2] <<- input[[..1]]}) 
        
      }
      
    })
    
    output$TBL1 <- renderDataTable(
      datatable(rv$df, selection = 'single', escape = F, filter = 'none'
                ,options = list(dom = 'frtip',pageLength = if(nrow(rv$df)>=10) {10} else {nrow(rv$df)} ,columnDefs = list(list(width = '17%', targets = 0)))
                
      ) %>% formatStyle(0,  fontWeight = 'bold')
    )
    
    
    
    #------------------------------------2.2 执行分析---------------------------
    # 2.2.1 单击将记录保存至原有文件
    seqinfo <- reactiveValues()
    results_path <- reactiveValues()
    file_num <- reactiveValues()
    results_pathwindows <- reactiveValues()
    
    observe({
      if(input$exc > 0){
        newlines <- isolate(data.frame(rv$df))
        file_num$ndf <- unique(newlines)
        # Generate results path
        results_path$path <- paste0("./RpNGS/data/results/",unique(newlines[1]))
        results_pathwindows$path <- paste0("C:/Result/OutputFq/",unique(newlines[1]),"/L01/")
        # update sequencinginfo.csv 
        if(file.exists('./RpNGS/data/response/RpNGS_sequencinginfo.csv')) {
          seqinfo$df <- read.csv('./RpNGS/data/response/RpNGS_sequencinginfo.csv')
          isolate(seqinfo$df <- unique(rbind(seqinfo$df,newlines)))
          write_excel_csv(seqinfo$df, "./RpNGS/data/response/RpNGS_sequencinginfo.csv")
          #write_excel_csv(file_num$ndf, "./RpNGS/data/response/temp/RpNGS_batch.csv")
        }else write_excel_csv(file_num$ndf, "./RpNGS/data/response/RpNGS_sequencinginfo.csv")
        
        # update clinicalinfo.csv 
        if(file.exists('./RpNGS/data/response/RpNGS_clinicalinfo.csv')) {
          clinicalinfo <- read.csv('./RpNGS/data/response/RpNGS_clinicalinfo.csv')
          sampleclinical_df <- data.frame("Sample_id" = unique(newlines[2]),"Name"= "","Gender"="",
                                          "Age"="","Sampling_day"="","Test_day"="","Tel"="","Infections"="","Requesting_physician"="",
                                          "Sample_type"="","Sample_status"="","Requesting_apartment"="","Requesting_hospital"="","Location"="",
                                          "Symptoms"="","Concerning_pathogens"="","Pretreatment"="","Pathogens"="") 
          clinicalinfo <- unique(rbind(clinicalinfo,sampleclinical_df))
          write_excel_csv(clinicalinfo, "./RpNGS/data/response/RpNGS_clinicalinfo.csv")
        }else write_excel_csv(sampleclinical_df, "./RpNGS/data/response/RpNGS_clinicalinfo.csv")
        
        
        # make a fastq folder and copy fastq files from sequencer 
        for ( sj in 1:17) {
          if(!dir.exists(paste0(results_path$path,"/","fq"))){
            #Sys.sleep(18000)
            system(paste0("mkdir ",results_path$path))
            system(paste0("mkdir ",results_path$path,"/fq"))
          } else { 
            for(i in 1:dim(file_num$ndf)[1]){
              print(i)
              print(file_num$ndf[i,6])
              system(paste0("scp NWU_M@192.168.14.50:",results_pathwindows$path,file_num$ndf[i,6]," ",paste0(results_path$path,"/","fq") ))
            }
            
            details <- file.info(list.files(path=paste0(results_path$path,"/fq"),pattern = ".fq.gz$", recursive = F,full.names=T))
            files <- rownames(details[with(details, order(as.POSIXct(mtime))), ])
            # filesfq <- gsub("^.*/", "", sub(".fq.gz","",files))
            filesfq <- gsub("^.*/", "", files)
            print(filesfq)
            if(length(filesfq) != dim(file_num$ndf)[1] ){
              print(length(filesfq) != length(file_num$ndf[6]))
              
              Sys.sleep(1800)
            } else {
              print(length(filesfq))
              print(length(file_num$ndf[6]))
              # Run a long task asynchronously
              future({
                # Simulate long-running task
                run_pngsanalysis(results_path$path,paste0(results_path$path,"/","fq"),72)
                "Data analysis completed！"
              }) %...>%  # Handle the result when it finishes
                (function(result) {
                  output$each_step <- renderText(result)
                })
              
              break
            }
          }
          
          
        }
        
        
        
      }
    })
    
    
    # 2.2.2 处理记录
    #（1）预执行记录条数
    output$filenum <- renderText({
      if(input$exc > 0){
        if(dim(file_num$ndf)[1]>1){
          paste0(dim(file_num$ndf)[1]," will be process after sequencing finished.")
        } else if (dim(file_num$ndf)[1]==1){
          paste0("Waiting for processing data analysis!")
        }
        
      }else paste0("Waiting for processing data analysis！")
    }) 
    
    
    
    #2.2.3 处理进度条
    
    fileraw <- reactiveValues()
    fileqc <- reactiveValues()
    filefiltered <- reactiveValues()
    fileclassified <- reactiveValues()
    fileannotated <- reactiveValues()
    
    observe({
      if(input$exc > 0){
        invalidateLater(1000)
        if(file.exists(paste0(results_path$path,"/fq"))){
          details <- file.info(list.files(path=paste0(results_path$path,"/fq"),pattern = ".fq.gz$", recursive = F,full.names=T))
          files <- rownames(details[with(details, order(as.POSIXct(mtime))), ])
          fileraw$fn <- gsub("^.*/", "", sub(".fq.gz","",files))
        }
        if(file.exists(paste0(results_path$path,"/01fastp"))){
          details <- file.info(list.files(path=paste0(results_path$path,"/01fastp"),pattern = ".fq.gz$", recursive = F,full.names=T))
          files <- rownames(details[with(details, order(as.POSIXct(mtime))), ])
          fileqc$fn <- gsub("^.*/", "", sub(".fq.gz","",files))
        }
        if(file.exists(paste0(results_path$path,"/02decomt"))){
          details <- file.info(list.files(path=paste0(results_path$path,"/02decomt"),pattern = ".fq.gz$", recursive = F,full.names=T))
          files <- rownames(details[with(details, order(as.POSIXct(mtime))), ])
          filefiltered$fn <- gsub("^.*/", "", sub(".fq.gz","",files))
        }
        if(file.exists(paste0(results_path$path,"/03kraken"))){
          details <- file.info(list.files(path=paste0(results_path$path,"/03kraken"),pattern = ".tsv$", recursive = F,full.names=T))
          files <- rownames(details[with(details, order(as.POSIXct(mtime))), ])
          fileclassified$fn <- gsub("^.*/", "", sub(".tsv","",files))
        }
        if(file.exists(paste0(results_path$path,"/04rawresults"))){
          details <- file.info(list.files(path=paste0(results_path$path,"/04rawresults"),pattern = ".csv$", recursive = F,full.names=T))
          files <- rownames(details[with(details, order(as.POSIXct(mtime))), ])
          fileannotated$fn <- gsub("^.*/", "", sub(".csv","",files))
        }
        
      }
    })
    
    
    output$each_step <- renderText({
      #paste0("This dataset has ",length(fileqc$fn)," sample runs and might take longer to create the correlation plot.")
      paste0("")
    })
    
    observe({
      updateProgressBar(id = "pb0", value = length(fileraw$fn), total =dim(file_num$ndf)[1]-1 ,title = paste("Transfer Data", trunc(length(fileraw$fn)/10)))
      updateProgressBar(id = "pb1", value = length(fileqc$fn), total =dim(file_num$ndf)[1]-1 ,title = paste("Quality Control", trunc(length(fileqc$fn)/10)))
      updateProgressBar(id = "pb2", value = length(filefiltered$fn), total =dim(file_num$ndf)[1]-1 ,title = paste("Host reads removal", trunc(length(filefiltered$fn)/10)))
      updateProgressBar(id = "pb3", value = length(fileclassified$fn), total =dim(file_num$ndf)[1]-1 ,title = paste("Classification", trunc(length(fileclassified$fn)/10)))
      updateProgressBar(id = "pb4", value = length(fileannotated$fn), total =dim(file_num$ndf)[1]-1 ,title = paste("Statical analysis", trunc(length(fileannotated$fn)/10)))
      
    })
    
    ################################################ 3. 临床报告 ################################################
    ############### 临床报告 ##############
    
    
    #---------------------------------------------------------------------------
    # search for samples id based on chip id 
    
    baogao_seqinfo <- reactiveValues()
    #clinicalinfo <- reactiveValues()
    #hide("sampleid_list")
    #hide("myDiv")
    observeEvent(input$chazhao,{
      reprot_chipid = trimws(input$reprot_chipid, which ="both")
      if(input$reprot_chipid !=""){
        sequencinginfo <- as.data.frame(fread(paste0("RpNGS/data/response/RpNGS_sequencinginfo.csv")),stringsAsFactors=F)
        chipid_samplesid_data<- sequencinginfo[which(sequencinginfo$Chip_id %in% reprot_chipid),c("Chip_id","Sample_id")]
        if(nrow(chipid_samplesid_data)==0){
          
          showModal(modalDialog(
            title = "Note！",
            paste("Unkown Chip ID"),
            easyClose = TRUE,
            footer = NULL))
          baogao_seqinfo$chipsampleid <- NULL
        } else {
          baogao_details <- file.info(list.files(path=paste0("./RpNGS/data/results/",reprot_chipid,"/06finalreports"),pattern = ".docx$", recursive = F,full.names=T))
          baogao_files <- rownames(baogao_details[with(baogao_details, order(as.POSIXct(mtime))), ])
          baogao_filesword <- gsub("^.*/", "", sub(".docx","",baogao_files))
          chipid_samplesid_data <- cbind(chipid_samplesid_data,"Check required")
          chipid_samplesid_data[chipid_samplesid_data$Sample_id %in% baogao_filesword,'"Check required"'] <- "Reported"
          colnames(chipid_samplesid_data) <- c("Chip_id","Sample_id","Report_status")
          baogao_seqinfo$chipsampleid <- unique(chipid_samplesid_data)
          
          
        }
        
      } else{
        
        showModal(modalDialog(
          title = "Note!",
          paste("Chip ID is required"),
          easyClose = TRUE,
          footer = NULL))
        hide("sampleid_list")
        baogao_seqinfo$chipsampleid <- NULL
      }
    })
    
    
    # USE DT TABLE 
    
    output$sampleid_list <- renderDataTable({
      baogao_seqinfo$chipsampleid
    },options = list(searchHighlight = TRUE, paging = FALSE,
                     scrollX = TRUE), rownames = FALSE,selection ="single")
    
    
    
    observeEvent(input$chazhao,{
      if(length(reactiveValuesToList(baogao_seqinfo)) !=0 ||input$reprot_chipid =="" ){
        hide("sampleid_list")
      } else
        show("sampleid_list" )
      
    })
    
    
    
    #---------------------------------------------------------------------------------   
    # show editable table of basic information about one selected sample
    output$sampleclinical_td <- renderDataTable({
      clinicalinfo <- as.data.frame(fread(paste0("RpNGS/data/response/RpNGS_clinicalinfo.csv")),stringsAsFactors=F)
      s = input$sampleid_list_rows_selected
      onelineclinical_df <-clinicalinfo[which(clinicalinfo$Sample_id %in% baogao_seqinfo$chipsampleid[s,2]),]
      
    },
    options = list(pageLength = 1, lengthChange = FALSE,dom = 't',  initComplete = JS(
      "function(settings, json) {",
      "$(this.api().table().header()).css({'background-color': '#696969', 'color': '#fff'});",
      "}")), 
    rownames = FALSE,
    selection ="none",
    class = 'cell-border stripe',
    editable = list(
      target = 'row', disable = list(columns = c(0))
      #target = 'cell', disable = list(columns = c(0))
    ))
    
    # save the edited table 
    observeEvent(input$sampleclinical_td_cell_edit, {
      clinicalinfo <- as.data.frame(fread(paste0("RpNGS/data/response/RpNGS_clinicalinfo.csv")),stringsAsFactors=F)
      s = input$sampleid_list_rows_selected
      onelineclinical_df <-clinicalinfo[which(clinicalinfo$Sample_id %in% baogao_seqinfo$chipsampleid[s,2]),]
      onelineclinical_df$Tel<- as.numeric(onelineclinical_df$Tel)
      onelineclinical_df$Sampleing_day <- as.Date(onelineclinical_df$Sampling_day)
      onelineclinical_df$Test_day <- as.Date(onelineclinical_df$Test_day)
      new_onelineclinical_df <<- editData(onelineclinical_df, input$sampleclinical_td_cell_edit, 'sampleclinical_td',rownames = FALSE)
      #
      
      clinicalinfo[which(clinicalinfo$Sample_id %in% new_onelineclinical_df$Sample_id),] <- new_onelineclinical_df
      write_excel_csv(clinicalinfo, "./RpNGS/data/response/RpNGS_clinicalinfo.csv")
    })
    #---------------------------------------------------------------------------
    # 2.4 pathogen table of each sample
    pathogen_human <- as.data.frame(fread(paste0("RpNGS/data/response/human_pathogens.csv")),stringsAsFactors=F) # all pre-known pathogens
    
    pathogeninfo <- reactive({ s = input$sampleid_list_rows_selected
    
    samples_id <- baogao_seqinfo$chipsampleid[s,2]
    if(file.exists(paste0("./RpNGS/data/results/",trimws(input$reprot_chipid, which ="both"),"/04rawresults/",samples_id,".csv"))){
      pathogeninfo <- as.data.frame(fread(paste0("./RpNGS/data/results/",trimws(input$reprot_chipid, which ="both"),"/04rawresults/",samples_id,".csv")),stringsAsFactors=F)
      
      return(pathogeninfo)
    } else {return(NULL)}
    
    })
    
    output$pathogen_td <- renderDataTable({
      pathogeninfo <- pathogeninfo()
      datatable(pathogeninfo,
                options = list(searchHighlight = TRUE,paging = FALSE),
                rownames = FALSE,
                selection ="multiple",
                # selection =list(row=c(1,3,5)),   
                class = 'cell-border stripe')  %>%formatStyle(1, valueColumns=1,color = '#337ab7', cursor = 'pointer') %>% formatRound(columns = c('Perct','ss_esti'),
                                                                                                                                       digits = 2) %>%
        formatStyle(7,
                    backgroundColor = styleEqual(pathogen_human[,1], "lightgreen"))
    })
    #%>% formatStyle(1,  fontWeight = 'bold') 
    
    #--------------------------------------------------
    # 2.5 select clinical pathogens from raw pathogen list
    selected_pathogen <- reactive({
      sel <- input$pathogen_td_rows_selected
      if(length(pathogeninfo())){
        pathogeninfo <- pathogeninfo()
        selected_pathogen <- pathogeninfo[sel,]
        return(selected_pathogen)
      }
    })
    
    
    output$pathogen_selected = DT::renderDataTable({
      selected_pathogen <- selected_pathogen()
      td_selected_pathogen <- selected_pathogen[,c("Reads","species","Kingdom")]
      
    },options = list(searchHighlight = TRUE,paging = FALSE,dom = 't',columnDefs = list(list(width = '17%', targets = 0))),
    rownames = FALSE,
    selection ="none", server = FALSE)
    
    #--------------------------------------------------
    # 2.6 Dynamic Plot of alignment from Selected Rows 
    # Load BAM file and calculate statistics
    # Reactive value to store BAM details
    bamDetails <- reactive({
      selected_rows <- input$pathogen_td_rows_selected
      if (length(selected_rows) == 0) {
        return(NULL)  # No rows selected, no plot
      } else {
        # Use the last selected row
        last_selected <- selected_rows[length(selected_rows)]
        pathogeninfo <- pathogeninfo()
        
        selected_data <- pathogeninfo[last_selected, "species"]
        # Load BAM file paths
        s = input$sampleid_list_rows_selected
        samples_id <- baogao_seqinfo$chipsampleid[s,2]
        bamPath <- paste0("D:/OneDrive - 陕西科技大学/dataanalysis/RpNGS/RpNGS/data/results/",trimws(input$reprot_chipid, which ="both"),"/05alignment/",samples_id,"_",selected_data,".bam")
        baiPath <- paste0("D:/OneDrive - 陕西科技大学/dataanalysis/RpNGS/RpNGS/data/results/",trimws(input$reprot_chipid, which ="both"),"/05alignment/",samples_id,"_",selected_data,".bam.bai")
        if (file.exists(bamPath) && file.exists(baiPath)){
          # Ensure BAM is indexed
          indexBam(bamPath, baiPath)
          
          # Extract BAM header to find the genome
          bamHeader <- scanBamHeader(bamPath)[[1]]$targets
          genomeName <- names(bamHeader)[1]
          genomeLength <- bamHeader[[1]]
          
          # Calculate statistics using GenomicAlignments
          bam <- readGAlignments(bamPath)
          cov <- coverage(bam)  # Returns an RleList
          
          # Extract coverage for the genome
          genomeCov <- as.numeric(cov[[genomeName]])  # Convert Rle to numeric for the contig
          
          # Calculate statistics
          avgCoverage <- mean(genomeCov)  # Average coverage
          coveredBases <- sum(genomeCov > 0)  # Count bases with coverage > 0
          percentCovered <- (coveredBases / genomeLength) * 100  # Percent covered
          totalReads <- length(bam)  # Total reads
          
          # Store details and stats
          bamDetails<-list(
            path = bamPath,
            species = selected_data,
            genomeName = genomeName,
            genomeLength = genomeLength,
            avgCoverage = avgCoverage,
            coveredBases = coveredBases,
            percentCovered = percentCovered,
            totalReads = totalReads
          )
          return(bamDetails)
        } else {return(NULL)}
      }
    })
    
    
    # Display statistics
    output$statsOutput <- renderPrint({
      if (length(bamDetails())){
        stats <- bamDetails()
        #if (length(stats)>0){
        cat("Species: ",stats$species,"\n")
        cat("Genome Name: ", stats$genomeName, "\n")
        cat("Genome Length: ", stats$genomeLength, "bp\n")
        cat("Average Coverage: ", round(stats$avgCoverage, 2), "x\n")
        cat("Covered Bases: ", stats$coveredBases, "bp\n")
        cat("Percent Covered: ", round(stats$percentCovered, 2), "%\n")
        cat("Total Reads: ", stats$totalReads, "\n")
        
      }
      
    })
    
    
    # Generate alignment plot
    output$alignmentPlot <- renderPlot({
      if (length(bamDetails())){
        stats <- bamDetails()
        # Create an AlignmentsTrack
        alignmentTrack <- AlignmentsTrack(
          stats$path, 
          isPaired = FALSE, 
          genome = "bacteria",options(ucscChromosomeNames=FALSE)
        )
        # Render the alignment plot
        plotTracks(
          list(alignmentTrack),
          chromosome = stats$genomeName,
          from = 1,
          to = stats$genomeLength,
          main = paste("Bowtie2 Alignments for", stats$species),
          background.title = "lightblue",
          fontcolor.title = "black"
        )
      }
      
    })
    
    
    #--------------------------------------------------
    # 2.7 To export a reactive table from Shiny to a pre-existing Word template
    
    # Generate the Word document on download
    output$download_doc <- downloadHandler(
      
      filename = function() {
        s = input$sampleid_list_rows_selected
        samples_id <- baogao_seqinfo$chipsampleid[s,2]
        paste0(samples_id, ".docx")
      },
      content = function(file) {
        create_word_document(selected_pathogen(),pathogeninfo(),pathogen_human,file)
      }
    )
    
    
    
    # Write a Function to Fill in the Template and Add the Table
    #create_word_document <- function(pathogen_data,alignment_data, output_file) {
    create_word_document <- function(pathogen_data,pathogeninfo,pathogen_human,output_file) {
      
      # Load the Word template
      #doc <- read_docx("template.docx")
      #doc <- read_docx("D:/OneDrive - 陕西科技大学/dataanalysis/RpNGS_lifegen/RpNGS/data/word_template/wordtemplate.docx")
      doc <- read_docx("D:/OneDrive - 陕西科技大学/dataanalysis/RpNGS/RpNGS/data/word_template/clinicalreporttemplate.docx")
      # Extract basic information of each sample
      clinicalinfo <- as.data.frame(fread(paste0("RpNGS/data/response/RpNGS_clinicalinfo.csv")),stringsAsFactors=F)
      s = input$sampleid_list_rows_selected
      onelineclinical_df <-clinicalinfo[which(clinicalinfo$Sample_id %in% baogao_seqinfo$chipsampleid[s,2]),]
      # add basic information 
      doc <- headers_replace_all_text(doc, "YBBH", as.character(onelineclinical_df[1,1]))
      doc <- body_replace_all_text(doc, "xingming", as.character(onelineclinical_df[1,2]))
      doc <- body_replace_all_text(doc, "xingbie", as.character(onelineclinical_df[1,3]))
      doc <- body_replace_all_text(doc, "nianling", as.character(onelineclinical_df[1,4]))
      doc <- body_replace_all_text(doc, "CYRQ", as.character(onelineclinical_df[1,5]))
      doc <- body_replace_all_text(doc, "jieyangriqi", as.character(onelineclinical_df[1,6]))
      doc <- body_replace_all_text(doc, "dianhua", as.character(onelineclinical_df[1,7]))
      doc <- body_replace_all_text(doc, "songjianyisheng", as.character(onelineclinical_df[1,8]))
      doc <- body_replace_all_text(doc, "BBLX", as.character(onelineclinical_df[1,9]))
      doc <- body_replace_all_text(doc, "yangbenzhuangtai", as.character(onelineclinical_df[1,10]))
      doc <- body_replace_all_text(doc, "songjiankeshi", as.character(onelineclinical_df[1,11]))
      doc <- body_replace_all_text(doc, "songjianyiyuan", as.character(onelineclinical_df[1,12]))
      doc <- body_replace_all_text(doc, "zhusu", as.character(onelineclinical_df[1,13]))
      doc <- body_replace_all_text(doc, "bingyuanleixing", as.character(onelineclinical_df[1,14]))
      doc <- body_replace_all_text(doc, "kangganran", as.character(onelineclinical_df[1,15]))
      
      
      
      
      #----------------------------------------------------------
      # 一、人工筛选后的表格
      # Create a flextable object from the reactive data
      selectedpathogen_data <- pathogen_data[,c("菌属中文名","菌种中文名","Gram+/-","species","Reads","Kingdom")] 
      
      colnames(selectedpathogen_data) <- c("菌属中文名","菌种中文名","类型","拉丁名","序列数","Kingdom")
      selectedpathogen_data$覆盖度 <- ""
      print(ncol(selectedpathogen_data))
      # ---
      # 1. avgcoverage of selelction pathogen
      
      for (var1 in selectedpathogen_data$拉丁名) {
        bamPath <- paste0("D:/OneDrive - 陕西科技大学/dataanalysis/RpNGS/RpNGS/data/results/",trimws(input$reprot_chipid, which ="both"),"/05alignment/",baogao_seqinfo$chipsampleid[s,2],"_",var1,".bam")
        baiPath <- paste0("D:/OneDrive - 陕西科技大学/dataanalysis/RpNGS/RpNGS/data/results/",trimws(input$reprot_chipid, which ="both"),"/05alignment/",baogao_seqinfo$chipsampleid[s,2],"_",var1,".bam.bai")
        if (file.exists(bamPath) && file.exists(baiPath)){
          # Ensure BAM is indexed
          indexBam(bamPath, baiPath)
          
          # Extract BAM header to find the genome
          bamHeader <- scanBamHeader(bamPath)[[1]]$targets
          genomeName <- names(bamHeader)[1]
          genomeLength <- bamHeader[[1]]
          
          # Calculate statistics using GenomicAlignments
          bam <- readGAlignments(bamPath)
          cov <- coverage(bam)  # Returns an RleList
          
          # Extract coverage for the genome
          genomeCov <- as.numeric(cov[[genomeName]])  # Convert Rle to numeric for the contig
          
          # Calculate statistics
          avgCoverage <- mean(genomeCov)  # Average coverage
          coveredBases <- sum(genomeCov > 0)  # Count bases with coverage > 0
          #percentCovered <- (coveredBases / genomeLength) * 100  # Percent covered
          percentCovered <- percent((coveredBases / genomeLength),accuracy = 0.01)   # Percent covered
          totalReads <- length(bam)  # Total reads
          
          # Store details and stats
          selectedpathogen_data[selectedpathogen_data$拉丁名 ==var1,"覆盖度"] <- percentCovered
        } 
      }
      
      
      
      
      
      #----------------------------------------------------------
      # 二、检测结果综述
      
      # Aggregate the data by group
      aggregated_data <- selectedpathogen_data %>%
        group_by(Kingdom) %>%
        summarise(
          merged_data = paste( 拉丁名, 序列数,"reads", collapse = "；"),
          .groups = 'drop'
        )
      print(aggregated_data)
      
      # repalce kingdom to chinese
      aggregated_data$Kingdom[aggregated_data$Kingdom %in% "Bacteria"] <- "Bacteria" 
      aggregated_data$Kingdom[aggregated_data$Kingdom %in% "Eukaryota"] <- "Eukaryota" 
      aggregated_data$Kingdom[aggregated_data$Kingdom %in% "Viruses"] <- "Viruses" 
      
      # # Extract only the 'merged_data' column
      # export_data <- aggregated_data %>%
      #   select(merged_data)
      
      # Aggregate into a single character vector
      aggregated_vector <- paste(aggregated_data$merged_data, collapse = "；")
      
      # Update RpNGS_clinicalinfo.csv
      if(dim(selectedpathogen_data)[1]!=0){
        clinicalinfo[which(clinicalinfo$Sample_id %in% baogao_seqinfo$chipsampleid[s,2]),"Pathogens"] <- aggregated_vector
      } else {
        clinicalinfo[which(clinicalinfo$Sample_id %in% baogao_seqinfo$chipsampleid[s,2]),"Pathogens"] <- "No Pathogens is Detected"
      }
      write_excel_csv(clinicalinfo, "./RpNGS/data/response/RpNGS_clinicalinfo.csv")
      
      # # Create a flextable without the header
      # ft_z <- flextable(export_data) %>%
      #   delete_part(part = "header")  # Remove the header row
      ft_z <- flextable(aggregated_data, col_keys = c("Summary"))%>%
        mk_par(j = "Summary", value = as_paragraph(Kingdom, "：",merged_data)) %>%
        delete_part(part = "header")  # Remove the header row
      
      
      # Set background color for body content to light orange
      ft_z <- bg(ft_z, part = "body", bg = "#FEF1E6") # Light orange color
      
      
      # Remove the vertical outer border by setting it to NULL
      # Set outer border properties
      ft_z <- border_outer(ft_z, border = fp_border(color = "#BFBFBF", width = 0.5))
      
      # Define no border for the vertical outer borders
      no_border <- fp_border(width = 0, color = "transparent")
      
      # Remove the vertical outer borders
      ft_z <- border(ft_z, part = "all", border.left = no_border, border.right = no_border)
      
      # # Set inner horizontal borders
      # ft_z <- border_inner_h(ft_z, border = fp_border(color = "#BFBFBF", width = 0.5))
      # # 
      # # # Set inner vertical borders
      # ft <- border_inner_v(ft, border = fp_border(color = "white", width = 0.5))
      
      
      # Adjust the table for better fit within the Word document
      #ft <- autofit(ft)
      #ft <- set_table_properties(ft,width = 1, layout = "autofit" )
      # Adjust the flextable width to 16.38 cm
      ft_z <- autofit(ft_z) # Fit content first
      ft_z <- flextable::width(ft_z, width = rep(16.38/2.54/(ncol(selectedpathogen_data)-6), (ncol(selectedpathogen_data)-6))) # Divide total width across columns
      # ncol(selectedpathogen_data) = 7
      
      
      # Set text alignment to centre for all columns
      ft_z <- align(ft_z, align = "left", part = "all") 
      
      # ft_z <- flextable(aggregated_data)%>%
      #   merge_v(j = ~Kingdom) %>% align(align = "left", part = "all")
      
      #mk_par(j = "综述", value = as_paragraph(菌种中文名, as_i(拉丁名))) 
      # set_header_labels(种名称 = "种名称")
      
      
      #----------------------------------------------------------
      # 三、检出结果列表
      
      
      #----------------------------------------------------------
      # 三、检出结果列表
      
      
      # ----------
      # 3.1 细菌检出列表
      # 3.1.1 细菌
      selected_bacteria <- selectedpathogen_data[which(selectedpathogen_data$Kingdom %in% "Bacteria"),c("菌属中文名","菌种中文名","类型","拉丁名","序列数","覆盖度")]
      selected_bacteria_b <- selected_bacteria[which(selected_bacteria$菌属中文名 != "分支杆菌属"),]
      colnames(selected_bacteria_b) <- c("Genus_Chinese","Species_Chinese","Types","Species_latin","Reads","Coverage")
      print(selected_bacteria_b)
      #print(dim(selected_bacteria_b))
      ft_b <- flextable(selected_bacteria_b, col_keys = c("Species","Types","Reads","Coverage"))%>%
        mk_par(j = "Species", value = as_paragraph(Species_Chinese, as_i(Species_latin))) %>%
        set_header_labels(Species = "Species")
      
      # value = as_paragraph(as_b(菌种中文名), as_i(拉丁名))
      # Set background color for header to white
      ft_b <- bg(ft_b, part = "header", bg = "white")
      
      # Set background color for body content to light orange
      ft_b <- bg(ft_b, part = "body", bg = "#FEF1E6") # Light orange color
      
      
      # Remove the vertical outer border by setting it to NULL
      # Set outer border properties
      ft_b <- border_outer(ft_b, border = fp_border(color = "#BFBFBF", width = 0.5))
      
      # Define no border for the vertical outer borders
      no_border <- fp_border(width = 0, color = "transparent")
      
      # Remove the vertical outer borders
      ft_b <- border(ft_b, part = "all", border.left = no_border, border.right = no_border)
      
      # # Set inner horizontal borders
      ft_b <- border_inner_h(ft_b, border = fp_border(color = "#BFBFBF", width = 0.5))
      # # 
      # # # Set inner vertical borders
      # ft <- border_inner_v(ft, border = fp_border(color = "white", width = 0.5))
      
      
      # Adjust the table for better fit within the Word document
      #ft <- autofit(ft)
      #ft <- set_table_properties(ft,width = 1, layout = "autofit" )
      # Adjust the flextable width to 16.38 cm
      ft_b <- autofit(ft_b) # Fit content first
      ft_b <- flextable::width(ft_b, width = rep(16.38/2.54/(ncol(selectedpathogen_data)-3), (ncol(selectedpathogen_data)-3))) # Divide total width across columns
      # ncol(selectedpathogen_data)-1 =5, ft_b also have 5 columns
      
      # Set text alignment to centre for all columns
      ft_b <- align(ft_b, align = "center", part = "all") 
      #---
      # 3.1.2结核
      selected_bacteria_mycobacterium <- selected_bacteria[which(selected_bacteria$菌属中文名 %in% "分支杆菌属"),]
      print(selected_bacteria_mycobacterium)
      colnames(selected_bacteria_mycobacterium) <- c("Genus_Chinese","Species_Chinese","Types","Species_latin","Reads","Coverage")
      # ft_m <- flextable(selected_bacteria_mycobacterium, col_keys = c("菌属中文名","种名称","类型","序列数","覆盖度"))%>%
      #   mk_par(j = "种名称", value = as_paragraph(as_b(菌种中文名), as_i(拉丁名))) %>%
      #   set_header_labels(种名称 = "种名称")
      
      #print(dim(selected_bacteria_b))
      ft_m <- flextable(selected_bacteria_mycobacterium, col_keys = c("Species","Types","Reads","Coverage"))%>%
        mk_par(j = "Species", value = as_paragraph(Species_Chinese, as_i(Species_latin))) %>%
        set_header_labels(Species = "Species")
      
      
      # Set background color for header to white
      ft_m <- bg(ft_m, part = "header", bg = "white")
      
      # Set background color for body content to light orange
      ft_m <- bg(ft_m, part = "body", bg = "#FEF1E6") # Light orange color
      
      
      # Remove the vertical outer border by setting it to NULL
      # Set outer border properties
      ft_m <- border_outer(ft_m, border = fp_border(color = "#BFBFBF", width = 0.5))
      
      # Define no border for the vertical outer borders
      no_border <- fp_border(width = 0, color = "transparent")
      
      # Remove the vertical outer borders
      ft_m <- border(ft_m, part = "all", border.left = no_border, border.right = no_border)
      
      # # Set inner horizontal borders
      ft_m <- border_inner_h(ft_m, border = fp_border(color = "#BFBFBF", width = 0.5))
      # # 
      # # # Set inner vertical borders
      # ft <- border_inner_v(ft, border = fp_border(color = "white", width = 0.5))
      
      # Adjust the flextable width to 16.38 cm
      ft_m <- autofit(ft_m) # Fit content first
      ft_m <- flextable::width(ft_m, width = rep(16.38/2.54/(ncol(selectedpathogen_data)-3), (ncol(selectedpathogen_data)-3))) # Divide total width across columns
      
      # Set text alignment to centre for all columns
      ft_m <- align(ft_m, align = "center", part = "all") 
      
      # ----------
      # 3.2 真核生物检出列表
      # 3.2.1 真菌
      selected_euk <- selectedpathogen_data[which(selectedpathogen_data$Kingdom %in% "Eukaryota"),c("菌属中文名","菌种中文名","类型","拉丁名","序列数","覆盖度")]
      selected_euk_euk <- selected_euk[which(selected_euk$类型 != "par"),]
      colnames(selected_euk_euk) <- c("Genus_Chinese","Species_Chinese","Types","Species_latin","Reads","Coverage")
      # ft_f <- flextable(selected_euk_euk, col_keys = c("菌属中文名","种名称","类型","序列数","覆盖度"))%>%
      #   mk_par(j = "种名称", value = as_paragraph(菌种中文名, as_i(拉丁名))) %>%
      #   set_header_labels(种名称 = "种名称")
      
      #print(dim(selected_bacteria_b))
      ft_f <- flextable(selected_euk_euk, col_keys = c("Species","Types","Reads","Coverage"))%>%
        mk_par(j = "Species", value = as_paragraph(Species_Chinese, as_i(Species_latin))) %>%
        set_header_labels(Species = "Species")
      
      
      # Set background color for header to white
      ft_f <- bg(ft_f, part = "header", bg = "white")
      
      # Set background color for body content to light orange
      ft_f <- bg(ft_f, part = "body", bg = "#FEF1E6") # Light orange color
      
      
      # Remove the vertical outer border by setting it to NULL
      # Set outer border properties
      ft_f <- border_outer(ft_f, border = fp_border(color = "#BFBFBF", width = 0.5))
      
      # Define no border for the vertical outer borders
      no_border <- fp_border(width = 0, color = "transparent")
      
      # Remove the vertical outer borders
      ft_f <- border(ft_f, part = "all", border.left = no_border, border.right = no_border)
      
      # # Set inner horizontal borders
      ft_f <- border_inner_h(ft_f, border = fp_border(color = "#BFBFBF", width = 0.5))
      # # 
      # # # Set inner vertical borders
      # ft <- border_inner_v(ft, border = fp_border(color = "white", width = 0.5))
      
      
      # Adjust the table for better fit within the Word document
      #ft <- autofit(ft)
      #ft <- set_table_properties(ft,width = 1, layout = "autofit" )
      # Adjust the flextable width to 16.38 cm
      ft_f <- autofit(ft_f) # Fit content first
      ft_f <- flextable::width(ft_f, width = rep(16.38/2.54/(ncol(selectedpathogen_data)-3), (ncol(selectedpathogen_data)-3))) # Divide total width across columns
      
      # Set text alignment to centre for all columns
      ft_f <- align(ft_f, align = "center", part = "all") 
      #---
      # 3.2.2 寄生虫
      selected_euk_par <- selected_euk[which(selected_euk$类型 %in% "par"),]
      
      colnames(selected_euk_par) <- c("Genus_Chinese","Species_Chinese","Types","Species_latin","Reads","Coverage")
      # ft_f <- flextable(selected_euk_euk, col_keys = c("菌属中文名","种名称","类型","序列数","覆盖度"))%>%
      #   mk_par(j = "种名称", value = as_paragraph(菌种中文名, as_i(拉丁名))) %>%
      #   set_header_labels(种名称 = "种名称")
      
      #print(dim(selected_bacteria_b))
      ft_p <- flextable(selected_euk_par, col_keys = c("Species","Types","Reads","Coverage"))%>%
        mk_par(j = "Species", value = as_paragraph(Species_Chinese, as_i(Species_latin))) %>%
        set_header_labels(Species = "Species")
      # 
      # ft_p <- flextable(selected_euk_par, col_keys = c("菌属中文名","种名称","类型","序列数","覆盖度"))%>%
      #   mk_par(j = "种名称", value = as_paragraph(菌种中文名, as_i(拉丁名))) %>%
      #   set_header_labels(种名称 = "种名称")
      # Set background color for header to white
      ft_p <- bg(ft_p, part = "header", bg = "white")
      
      # Set background color for body content to light orange
      ft_p <- bg(ft_p, part = "body", bg = "#FEF1E6") # Light orange color
      
      
      # Remove the vertical outer border by setting it to NULL
      # Set outer border properties
      ft_p <- border_outer(ft_p, border = fp_border(color = "#BFBFBF", width = 0.5))
      
      # Define no border for the vertical outer borders
      no_border <- fp_border(width = 0, color = "transparent")
      
      # Remove the vertical outer borders
      ft_p <- border(ft_p, part = "all", border.left = no_border, border.right = no_border)
      
      # # Set inner horizontal borders
      ft_p <- border_inner_h(ft_p, border = fp_border(color = "#BFBFBF", width = 0.5))
      # # 
      # # # Set inner vertical borders
      # ft <- border_inner_v(ft, border = fp_border(color = "white", width = 0.5))
      
      
      # Adjust the table for better fit within the Word document
      #ft <- autofit(ft)
      #ft <- set_table_properties(ft,width = 1, layout = "autofit" )
      # Adjust the flextable width to 16.38 cm
      ft_p <- autofit(ft_p) # Fit content first
      ft_p <- flextable::width(ft_p, width = rep(16.38/2.54/(ncol(selectedpathogen_data)-3), (ncol(selectedpathogen_data)-3))) # Divide total width across columns
      
      # Set text alignment to centre for all columns
      ft_p <- align(ft_p, align = "center", part = "all") 
      
      
      # ----------
      # 3.3 病毒检出列表
      selected_vir <- selectedpathogen_data[which(selectedpathogen_data$Kingdom %in% "Viruses"),c("菌属中文名","菌种中文名","类型","拉丁名","序列数","覆盖度")]
      colnames(selected_vir) <- c("Genus_Chinese","Species_Chinese","Types","Species_latin","Reads","Coverage")
      # ft_f <- flextable(selected_euk_euk, col_keys = c("菌属中文名","种名称","类型","序列数","覆盖度"))%>%
      #   mk_par(j = "种名称", value = as_paragraph(菌种中文名, as_i(拉丁名))) %>%
      #   set_header_labels(种名称 = "种名称")
      
      #print(dim(selected_bacteria_b))
      ft_v <- flextable(selected_vir, col_keys = c("Species","Types","Reads","Coverage"))%>%
        mk_par(j = "Species", value = as_paragraph(Species_Chinese, as_i(Species_latin))) %>%
        set_header_labels(Species = "Species")
      
      # 
      # ft_v <- flextable(selected_vir, col_keys = c("菌属中文名","种名称","类型","序列数","覆盖度"))%>%
      #   mk_par(j = "种名称", value = as_paragraph(菌种中文名, as_i(拉丁名))) %>%
      #   set_header_labels(种名称 = "种名称")
      
      # Set background color for header to white
      ft_v <- bg(ft_v, part = "header", bg = "white")
      
      # Set background color for body content to light orange
      ft_v <- bg(ft_v, part = "body", bg = "#FEF1E6") # Light orange color
      
      
      # Remove the vertical outer border by setting it to NULL
      # Set outer border properties
      ft_v <- border_outer(ft_v, border = fp_border(color = "#BFBFBF", width = 0.5))
      
      # Define no border for the vertical outer borders
      no_border <- fp_border(width = 0, color = "transparent")
      
      # Remove the vertical outer borders
      ft_v <- border(ft_v, part = "all", border.left = no_border, border.right = no_border)
      
      # # Set inner horizontal borders
      ft_v <- border_inner_h(ft_v, border = fp_border(color = "#BFBFBF", width = 0.5))
      # # 
      # # # Set inner vertical borders
      # ft <- border_inner_v(ft, border = fp_border(color = "white", width = 0.5))
      
      
      # Adjust the table for better fit within the Word document
      #ft <- autofit(ft)
      #ft <- set_table_properties(ft,width = 1, layout = "autofit" )
      # Adjust the flextable width to 16.38 cm
      ft_v <- autofit(ft_v) # Fit content first
      ft_v <- flextable::width(ft_v, width = rep(16.38/2.54/(ncol(selectedpathogen_data)-3), (ncol(selectedpathogen_data)-3))) # Divide total width across columns
      
      # Set text alignment to centre for all columns
      ft_v <- align(ft_v, align = "center", part = "all") 
      
      # ----------
      # 3.4 疑似检出列表
      #print(pathogeninfo[which(pathogeninfo[,"species"] %in% pathogen_human[,1]),])
      d_pathogens <- pathogeninfo[which(pathogeninfo[,"species"] %in% pathogen_human[,1]),]
      d_pathogens <- d_pathogens[-which(d_pathogens[,"species"] %in% pathogen_data[,"species"]),]
      d_pathogens <- d_pathogens[,c("菌属中文名","菌种中文名","Gram+/-","species","Reads")] 
      
     # colnames(d_pathogens) <- c("菌属中文名","菌种中文名","类型","拉丁名","序列数","名词解释")
      # d_pathogens$种名称 <- paste(d_pathogens$菌种中文名,d_pathogens$species)
      # d_pathogens <- d_pathogens[,c("种名称","Gram+/-","Reads","名词解释")]  
      # colnames(d_pathogens) <- c("种名称","类型","序列数","名词解释")
      colnames(d_pathogens) <- c("Genus_Chinese","Species_Chinese","Types","Species_latin","Reads")
      
      ft_y <- flextable(d_pathogens, col_keys = c("Species","Types","Reads"))%>%
        mk_par(j = "Species", value = as_paragraph(Species_Chinese, as_i(Species_latin))) %>%
        set_header_labels(Species = "Species")
      
      # ft_y <- flextable(d_pathogens, col_keys = c("种名称","类型","序列数","名词解释"))%>%
      #   mk_par(j = "种名称", value = as_paragraph(菌种中文名, as_i(拉丁名))) %>%
      #   set_header_labels(种名称 = "种名称")
      
      #ft_y <- flextable(d_pathogens)
      # Set background color for header to white
      ft_y <- bg(ft_y, part = "header", bg = "white")
      
      # Set background color for body content to light orange
      ft_y <- bg(ft_y, part = "body", bg = "#FEF1E6") # Light orange color
      
      
      # Remove the vertical outer border by setting it to NULL
      # Set outer border properties
      ft_y <- border_outer(ft_y, border = fp_border(color = "#BFBFBF", width = 0.5))
      
      # Define no border for the vertical outer borders
      no_border <- fp_border(width = 0, color = "transparent")
      
      # Remove the vertical outer borders
      ft_y <- border(ft_y, part = "all", border.left = no_border, border.right = no_border)
      
      # # Set inner horizontal borders
      ft_y <- border_inner_h(ft_y, border = fp_border(color = "#BFBFBF", width = 0.5))
      # # 
      # # # Set inner vertical borders
      # ft <- border_inner_v(ft, border = fp_border(color = "white", width = 0.5))
      
      
      # Adjust the table for better fit within the Word document
      #ft <- autofit(ft)
      #ft <- set_table_properties(ft,width = 1, layout = "autofit" )
      # Adjust the flextable width to 16.38 cm
      ft_y <- autofit(ft_y) # Fit content first
      # ft_y <- flextable::width(ft_y, width = rep(16.38/2.54/(ncol(selectedpathogen_data)-3), (ncol(selectedpathogen_data)-3))) # Divide total width across columns
      ft_y <- flextable::width(ft_y,j=1, width = 4.45) # separately set column widths
      ft_y <- flextable::width(ft_y,j=2, width = 1)
      ft_y <- flextable::width(ft_y,j=3, width = 1)
      #ft_y <- flextable::width(ft_y,j=4, width = 3.16) 
      
      # Set text alignment to centre for all columns
      ft_y <- align(ft_y, align = "center", part = "all") 
      
      
      
      
      
      #----------------------------------------------------------
      # 4.  
      # To export the cell values of a flextable as numbered content into a pre-existing Word template,
      # Iterate over each row of the flextable's data
      
      # This work but \n is not work; can not change to a new line
      # Build numbered list content as a single string
      # list_content <- paste0(seq_along(pathogen_data[,'菌种中文名']), "). ",pathogen_data[,'菌种中文名'],":",pathogen_data$名词解释, collapse = "\n")
      # # Replace placeholder with the numbered list
      # doc <- doc %>%
      #   body_replace_all_text("JSSM", list_content)
      
      
      # This code can add first item in the right position, but the rest are still append in the final of word. 
      # # Replace placeholder with the first line of the list (to initialize the position)
      # doc <- doc %>%
      #   body_replace_all_text("JSSM", paste0("1. ", pathogen_data[1,'菌种中文名'],":",pathogen_data[1,"名词解释"]))
      # 
      # # Append the rest of the list line by line
      # for (i in 2:length(pathogen_data)) {
      #   doc <- doc %>% body_add_par(value = paste0(i, ". ", pathogen_data[i,'菌种中文名'],":",pathogen_data[i,"名词解释"]), pos = "after")
      # }
      
      
      #----------------------------------------------------------
      # 四、测序质量
      
      #---
      # 4.1 加载实验质控
      sequencinginfo <- as.data.frame(fread(paste0("RpNGS/data/response/RpNGS_sequencinginfo.csv")),stringsAsFactors=F)
      onesequencinginfo <- sequencinginfo[which(sequencinginfo$样本编号 %in% baogao_seqinfo$chipsampleid[s,2]),]
      doc <- body_replace_all_text(doc, "HSND", as.character(onesequencinginfo[1,3]))
      doc <- body_replace_all_text(doc, "WKND", as.character(onesequencinginfo[1,5]))
      
      #---
      # 4.2 数据质控
      #---
      
      #--
      # fastp file 
      #html_file <- paste0('D:/OneDrive - 陕西科技大学/dataanalysis/RpNGS/RpNGS/data/results/','FT10005323/qc/lifegen20240830_8.html')
      html_file<- paste0("D:/OneDrive - 陕西科技大学/dataanalysis/RpNGS/RpNGS/data/results/",trimws(input$reprot_chipid, which ="both"),"/01fastp/",baogao_seqinfo$chipsampleid[s,2],".html")
      
      # Extract Filtered Reads (total reads after filtering)
      
      # Read the HTML file
      html_content <- read_html(html_file)
      # using css info to extract table. 
      quality_table <- html_content%>%
        html_nodes("div#after_filtering_summary table") %>%
        html_table(fill=TRUE)
      # convert tibble list to data.frame 
      quality_table <- as.data.frame(do.call(rbind, quality_table))
      print(quality_table)
      # Extract Q30 bases and total reads
      q30_bases <- quality_table[quality_table$X1 == "Q30 bases:", "X2"]
      print(q30_bases)
      # Extract the percentage using regex
      q30_bases <- sub(".*\\((\\d+\\.\\d+)%\\).*", "\\1", q30_bases)
      
      # Convert to numeric, round to 2 decimal places, and append "%"
      q_bases <- paste0(round(as.numeric(q30_bases), 2), "%")
      
      
      total_reads <- quality_table[quality_table$X1 == "total reads:", "X2"]
      doc <- body_replace_all_text(doc, "ZXLS", as.character(total_reads))
      doc <- body_replace_all_text(doc, "QSLB", as.character(q_bases))
      # 
      
      
      #--
      #Human reads filtered file 
      fastq_file <- paste0("D:/OneDrive - 陕西科技大学/dataanalysis/RpNGS/RpNGS/data/results/",trimws(input$reprot_chipid, which ="both"),"/02decomt/",baogao_seqinfo$chipsampleid[s,2],".fastq.gz")
      num_reads <- countFastq(fastq_file)
      doc <- body_replace_all_text(doc, "FRYXL", as.character(num_reads$scores))
      
      
      
      
      
      # Add the table at the bookmark
      if(dim(aggregated_data)[1]!=0){
        doc <- body_replace_flextable_at_bkm(doc, ft_z, bookmark = "ft_z")
      }
      
      if(dim(selected_bacteria_b)[1]!=0){
        doc <- body_replace_flextable_at_bkm(doc, ft_b, bookmark = "ft_b")
      }
      if(dim(selected_bacteria_mycobacterium)[1]!=0){
        doc <- body_replace_flextable_at_bkm(doc, ft_m, bookmark = "ft_m")
      }
      if(dim(selected_euk_euk)[1]!=0){
        doc <- body_replace_flextable_at_bkm(doc, ft_f, bookmark = "ft_f")
      }
      if(dim(selected_euk_par)[1]!=0){
        doc <- body_replace_flextable_at_bkm(doc, ft_p, bookmark = "ft_p")
      }
      if(dim(selected_vir)[1]!=0){
        doc <- body_replace_flextable_at_bkm(doc, ft_v, bookmark = "ft_v")
      }
      if(dim(d_pathogens)[1]!=0){
        doc <- body_replace_flextable_at_bkm(doc, ft_y, bookmark = "ft_y")
      }
      
      
      
      
      
      
      # Save the new Word document
      print(doc, target = output_file)
    }
    #==============================
    
    
    
  }
  shiny::runApp(shiny::shinyApp(ui, server), quiet=FALSE, launch.browser=TRUE)
  #shinyApp(ui, server)
}

RpNGS()


