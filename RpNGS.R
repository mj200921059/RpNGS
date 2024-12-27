

# Load packages 
source('./R/packages.R')
# Perset
source('./R/perset.R')
# File transfer function.
#source("./R/seqer_fastq2server.R")
# Activate data process pipeline
source('./R/rpngs_se_analysis.R')
# export pathogen, sequencing info and clinical info to word as a clinical detection report.
source('./R/export2word.R')

# Define UI for application 
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
                                                         includeMarkdown(help_file)  # Replace with your markdown file
                                         )
                                         

                                 )
                                 ,column(1
                                         
                                 )
                                 
                        ),
                        tabPanel(h4('Summary dataset'), value='startsets'
                                 ,includeMarkdown(title_file)
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
                                 ,includeHTML(footer_file)
                                 
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

# Define server
server <- function(input, output, session){
  
  ################################################ Main data for pNGS datatable and search in section 1 ################################################
  
  datafr1 <- as.data.frame(fread(seqinfo_file),stringsAsFactors=F)
  datafr2 <- as.data.frame(fread(clininfo_file),stringsAsFactors=F) #FOR SUMMARY SAMPLE NUMBER, HAVE TO FILL A CITY OF EACH SAMPLE
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
  c_map <- read_sf(map_file)
  # sf convert to dataframe
  c_map <- as.data.frame(c_map)
  # Merge the observation data and map data according to administrative division codes.
  shaanxi <- merge(c_map, shaanxi_data, by = "adcode")
  # dataframe convert to sf
  shaanxi <- sf::st_as_sf(shaanxi, sf_column_name = "geometry")
  # Set a segmented discrete color palette
  # pal <- colorBin("Spectral", bins = pretty(shaanxi$pop), reverse = TRUE)
  pal <- colorNumeric("Spectral", domain = NULL)
  
  
  
  ####========location_status ========####
  
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
  ####======== sampletype_stat ========####
  
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
  
  
  
  ####======= main data table ========####
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
      file_path <- paste0(path4results,datafr1[which(datafr1$Sample_id %in% info), "Chip_id"],"/06finalreports/",info,".docx")  # Replace with the actual file path
      print(1)
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
  
  ################################################ 2. Analysis dataset in section 2################################################
  
  ####========2.1 Input experimental datasets ========####
  
  # Provide an example for NGS input info.
  rv <- reactiveValues(df = Input, row_selected = NULL) 
  
  # For input data
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
          
        )}
  })
  # For delete existing rows in TBL1 Table
  observeEvent(input$deleteRows,{
    
    if (!is.null(input$TBL1_rows_selected)) {
      #data(data()[-as.numeric(input$TBL1_rows_selected),])
      rv$df <- rv$df[-as.numeric(input$TBL1_rows_selected), ]
    }
  })
  
  # For edit row in TBL1 Table 
  observeEvent(input$edit,{
    
    if (!is.null(input$TBL1_rows_selected)) {
      cols_to_edit <- c('chip_id', 'sample_id', 'nucli_conce', 'adaptor', 'dna_conce','seq_id')
      colnms <- c('Chip_id', 'Sample_id', 'Extracted_NAC', 'Adaptor', 'library_NAC','Rawfastq_id')
      "remember the row selected"
      rv$row_selected <- input$TBL1_rows_selected
      
      walk2(cols_to_edit, colnms, ~{rv$df[input$TBL1_rows_selected, ..2] <<- input[[..1]]}) 
      
    }
    
  })
  
  # Show Table 
  output$TBL1 <- renderDataTable(
    datatable(rv$df, selection = 'single', escape = F, filter = 'none'
              ,options = list(dom = 'frtip',pageLength = if(nrow(rv$df)>=10) {10} else {nrow(rv$df)} ,columnDefs = list(list(width = '17%', targets = 0)))
              
    ) %>% formatStyle(0,  fontWeight = 'bold')
  )
  
  ####========2.2 Process data analysis ========####
  
  seqinfo <- reactiveValues()
  results_path <- reactiveValues()
  file_num <- reactiveValues()
  results_pathwindows <- reactiveValues()
  
  # Update the RpNGS_sequencinginfo.csv file and run a data anlaysis.
  observe({
    if(input$exc > 0){
      newlines <- isolate(data.frame(rv$df))
      newlines <- newlines[-1,] #remove example info 
      file_num$ndf <- unique(newlines)
      # SET results path
      results_path$path <- paste0(path4results,unique(newlines[1]))
      results_pathwindows$path <- paste0(path4sequencer,unique(newlines[1]),"/L01/")
      
      # update sequencinginfo.csv 
      if(file.exists(seqinfo_file)) {
        seqinfo$df <- read.csv(seqinfo_file)
        isolate(seqinfo$df <- unique(rbind(seqinfo$df,newlines)))
        write_excel_csv(seqinfo$df, seqinfo_file)
      }else write_excel_csv(file_num$ndf, seqinfo_file)
      
      # update clinicalinfo.csv 
      if(file.exists(clininfo_file)) {
        clinicalinfo <- read.csv(clininfo_file)
        sampleclinical_df <- data.frame("Sample_id" = unique(newlines[2]),"Name"= "","Gender"="",
                                        "Age"="","Sampling_day"="","Test_day"="","Tel"="","Infections"="","Requesting_physician"="",
                                        "Sample_type"="","Sample_status"="","Requesting_apartment"="","Requesting_hospital"="","Location"="",
                                        "Symptoms"="","Concerning_pathogens"="","Pretreatment"="","Pathogens"="") 
        clinicalinfo <- unique(rbind(clinicalinfo,sampleclinical_df))
        write_excel_csv(clinicalinfo, clininfo_file)
      }else write_excel_csv(sampleclinical_df, clininfo_file)
      
      # 
      # make a fastq folder and copy fastq files from sequencer 
      for ( sj in 1:17) {
        if(!dir.exists(paste0(results_path$path,"/","fq"))){
          system(paste0("mkdir ",results_path$path))
          system(paste0("mkdir ",results_path$path,"/fq"))
        } else { 
          for(i in 1:dim(file_num$ndf)[1]){
            print(i)
            print(file_num$ndf[i,6])
            system(paste0("scp ",seqer_ip,":",results_pathwindows$path,file_num$ndf[i,6]," ",paste0(results_path$path,"/","fq") ))
          }
          
          details <- file.info(list.files(path=paste0(results_path$path,"/fq"),pattern = ".fq.gz$", recursive = F,full.names=T))
          files <- rownames(details[with(details, order(as.POSIXct(mtime))), ])
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
              # Long-running task
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
  
  ####========2.2.2 LOG info========####

  #（1）Sample Number 
  output$filenum <- renderText({
    if(input$exc > 0){
      if(dim(file_num$ndf)[1]>1){
        paste0(dim(file_num$ndf)[1]," will be process after sequencing finished.")
      } else if (dim(file_num$ndf)[1]==1){
        paste0("Waiting for processing data analysis!")
      }
      
    }else paste0("Waiting for processing data analysis！")
  }) 
  
  
  ####========2.2.3 Progress bar ========####

  
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
    updateProgressBar(id = "pb4", value = length(fileannotated$fn), total =dim(file_num$ndf)[1]-1 ,title = paste("Validation analysis", trunc(length(fileannotated$fn)/10)))
    
  })
  
  ################################################ 3. Test report ################################################

  
  ####========3.1  search for samples id based on chip id  ========####  

  baogao_seqinfo <- reactiveValues()
  observeEvent(input$chazhao,{
    reprot_chipid = trimws(input$reprot_chipid, which ="both")
    if(input$reprot_chipid !=""){
      sequencinginfo <- as.data.frame(fread(seqinfo_file),stringsAsFactors=F)
      chipid_samplesid_data<- sequencinginfo[which(sequencinginfo$Chip_id %in% reprot_chipid),c("Chip_id","Sample_id")]
      if(nrow(chipid_samplesid_data)==0){
        
        showModal(modalDialog(
          title = "Note！",
          paste("Unkown Chip ID"),
          easyClose = TRUE,
          footer = NULL))
        baogao_seqinfo$chipsampleid <- NULL
      } else {
        baogao_details <- file.info(list.files(path=paste0(path4results,reprot_chipid,"/06finalreports"),pattern = ".docx$", recursive = F,full.names=T))
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
    clinicalinfo <- as.data.frame(fread(clininfo_file),stringsAsFactors=F)
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
    clinicalinfo <- as.data.frame(fread(clininfo_file),stringsAsFactors=F)
    s = input$sampleid_list_rows_selected
    onelineclinical_df <-clinicalinfo[which(clinicalinfo$Sample_id %in% baogao_seqinfo$chipsampleid[s,2]),]
    onelineclinical_df$Tel<- as.numeric(onelineclinical_df$Tel)
    onelineclinical_df$Sampleing_day <- as.Date(onelineclinical_df$Sampling_day)
    onelineclinical_df$Test_day <- as.Date(onelineclinical_df$Test_day)
    new_onelineclinical_df <<- editData(onelineclinical_df, input$sampleclinical_td_cell_edit, 'sampleclinical_td',rownames = FALSE)
    #
    
    clinicalinfo[which(clinicalinfo$Sample_id %in% new_onelineclinical_df$Sample_id),] <- new_onelineclinical_df
    write_excel_csv(clinicalinfo, clininfo_file)
  })
  #---------------------------------------------------------------------------
  # 2.4 pathogen table of each sample
  pathogen_human <- as.data.frame(fread(paste0(path4datasets,"human_pathogens.csv")),stringsAsFactors=F) # all pre-known pathogens
  
  pathogeninfo <- reactive({ s = input$sampleid_list_rows_selected
  
  samples_id <- baogao_seqinfo$chipsampleid[s,2]
  if(file.exists(paste0(path4results,trimws(input$reprot_chipid, which ="both"),"/04rawresults/",samples_id,".csv"))){
    pathogeninfo <- as.data.frame(fread(paste0(path4results,trimws(input$reprot_chipid, which ="both"),"/04rawresults/",samples_id,".csv")),stringsAsFactors=F)
    
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
      bamPath <- paste0(path4results,trimws(input$reprot_chipid, which ="both"),"/05alignment/",samples_id,"_",selected_data,".bam")
      baiPath <- paste0(path4results,trimws(input$reprot_chipid, which ="both"),"/05alignment/",samples_id,"_",selected_data,".bam.bai")
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

}

# Run the application 
# shinyApp(ui = ui, server = server)
shiny::runApp(shiny::shinyApp(ui, server), quiet=FALSE, launch.browser=TRUE)
