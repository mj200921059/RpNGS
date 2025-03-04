# RpNGS v0.1.1
# 2025.02.03

#source("config.R")


#source('./R/packages.R')
# source('./R/perset.R')

# source('./R/export2word.R')
#---------------
# R functions 
source("./R/data_prep.R")  # 载入数据处理脚本
source('./R/run_pngsanalysis.R')

# Load modules 
source('./modules/help_module.R')
source('./modules/summary_module.R')
source('./modules/analysis_module.R')
#source('./modules/report_module.R')
# Define UI
ui <- shiny::navbarPage(
  h4("pNGS v0.1.1"), id = "rpngs", selected = "startsets",
  theme = shinythemes::shinytheme("flatly"),
  
  mod_help_ui("help"),  # 调用 Help 模块
  mod_summary_ui("summary"),  # 调用 Summary Dataset 模块
  mod_analysis_ui("analysis"), # 调用 analysis 模块
  #mod_report_ui("report") # 调用 report 模块
  
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
                                              h3("Raw microbes list:"),
                                              shiny::tags$style(HTML("
                                                table.dataTable tbody tr {
                                                              height:40px !important; 
                                                                 }
                                                table.dataTable td{
                                                  overflow:hidden;
                                                  test-overflow:ellipsis;
                                                  white-space:nowrap;
                                                  max-width:150px;
                                                }
                                              "))
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

# 定义服务器逻辑
server <- function(input, output, session) {
  mod_help_server("help")  # 绑定 Help 模块
  mod_summary_server("summary", datafr2, summary_sub_datafr2)  # 绑定 Summary Dataset 模块
  mod_analysis_server("analysis",sample_data,datafr1,datafr2) # 绑定 anslysis 模块
  #mod_report_server("report",clininfo_file,pathogen_human,results_folders,seqinfo_file) 
  #mod_report_server("report",report_data) 
  
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
        # Create file paths and check for existence
        chipid_samplesid_data <- chipid_samplesid_data %>%
          mutate(
            report_path = file.path(results_folders, Sample_id, "06finalreports", paste0(Sample_id,".docx")),
            Status = ifelse(file.exists(report_path), "Reported", "Check Required")
          ) %>%
          select(-report_path)  # Remove the temporary path column
        # View the updated data frame
        baogao_seqinfo$chipsampleid <- unique(chipid_samplesid_data)
        print(baogao_seqinfo$chipsampleid)
        
        
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
    target = 'row', disable = list(columns = c(0,17))
    #target = 'cell', disable = list(columns = c(0))
  ))
  
  # save the edited table
  observeEvent(input$sampleclinical_td_cell_edit, {
    clinicalinfo <- as.data.frame(fread(clininfo_file),stringsAsFactors=F)
    s = input$sampleid_list_rows_selected
    onelineclinical_df <-clinicalinfo[which(clinicalinfo$Sample_id %in% baogao_seqinfo$chipsampleid[s,2]),]
    onelineclinical_df$Tel<- as.numeric(onelineclinical_df$Tel)
    onelineclinical_df$Sampling_day <- as.Date(onelineclinical_df$Sampling_day)
    onelineclinical_df$Test_day <- as.Date(onelineclinical_df$Test_day)
    new_onelineclinical_df <<- editData(onelineclinical_df, input$sampleclinical_td_cell_edit, 'sampleclinical_td',rownames = FALSE)
    #
    
    clinicalinfo[which(clinicalinfo$Sample_id %in% new_onelineclinical_df$Sample_id),] <- new_onelineclinical_df
    write_excel_csv(clinicalinfo, clininfo_file)
  })
  #---------------------------------------------------------------------------
  # 2.4 pathogen table of each sample
  pathogen_human <- as.data.frame(fread(pathogen_human),stringsAsFactors=F) # all pre-known pathogens
  
  pathogeninfo <- reactive({ s = input$sampleid_list_rows_selected
  
  samples_id <- baogao_seqinfo$chipsampleid[s,2]
  print(samples_id)
  if(file.exists(paste0(results_folders,"/",samples_id,"/04rawresults/",samples_id,".csv"))){
    pathogeninfo <- as.data.frame(fread(paste0(results_folders,"/",samples_id,"/04rawresults/",samples_id,".csv")),stringsAsFactors=F)
    #pathogeninfo$zscore <-""
    return(pathogeninfo)
  } else {return(NULL)}
  
  })

  output$pathogen_td <- renderDataTable({
    pathogeninfo <- pathogeninfo()
    datatable(pathogeninfo,
              options = list(searchHighlight = TRUE,paging = FALSE,
                             rowCallback =JS(
                               "function(row,data,index){",
                               # " $('td',row).attr('title',data);", # This code hover all content of one row. bleow code only show content of specific column  
                               " $('td:eq(13),td:eq(7)',row).attr('title',function(){ return $(this).text(); });", 
                               "}")
                             ),
              rownames = FALSE,
              selection ="multiple",
              # selection =list(row=c(1,3,5)),
              class = 'cell-border stripe')  %>%formatStyle(1, valueColumns=1,color = '#337ab7', cursor = 'pointer') %>% formatRound(columns = c('Perct','zscore'),
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
      bamPath <- paste0(results_folders,"/",samples_id,"/05alignment/",samples_id,"_",selected_data,".bam")
      baiPath <- paste0(results_folders,"/",samples_id,"/05alignment/",samples_id,"_",selected_data,".bam.bai")
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
      # set a region insteads of the entire genome
      from_pos <- 1 # set a starting position
      to_pos <- min(stats$genomeLength,100000) # Limit to 100kb region to reduce load
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
        from = from_pos,
        to = to_pos,
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

#shinyApp(ui, server)
shiny::runApp(shiny::shinyApp(ui, server), quiet=FALSE, launch.browser=TRUE)
