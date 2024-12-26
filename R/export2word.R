

# Write a Function to Fill in the Template and Add the Table
create_word_document <- function(pathogen_data,pathogeninfo,pathogen_human,output_file) {
  
  # Load the Word template
  
  doc <- read_docx(paste0(path4template,"clinicalreporttemplate.docx"))
  # Extract basic information of each sample
  clinicalinfo <- as.data.frame(fread(clininfo_file),stringsAsFactors=F)
  s = input$sampleid_list_rows_selected
  onelineclinical_df <-clinicalinfo[which(clinicalinfo$Sample_id %in% baogao_seqinfo$chipsampleid[s,2]),]
  
  #---------------------------------
  # Section 1. add basic information 
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
  # Pathogens by manually selected from raw microbes list.
  # Create a flextable object from the reactive data
  selectedpathogen_data <- pathogen_data[,c("菌属中文名","菌种中文名","Gram+/-","species","Reads","Kingdom")] 
  
  colnames(selectedpathogen_data) <- c("菌属中文名","菌种中文名","类型","拉丁名","序列数","Kingdom")
  selectedpathogen_data$覆盖度 <- ""
  print(ncol(selectedpathogen_data))
  # ---
  # avgcoverage of selelction pathogen
  
  for (var1 in selectedpathogen_data$拉丁名) {
    bamPath <- paste0(path4results,trimws(input$reprot_chipid, which ="both"),"/05alignment/",baogao_seqinfo$chipsampleid[s,2],"_",var1,".bam")
    baiPath <- paste0(path4results,trimws(input$reprot_chipid, which ="both"),"/05alignment/",baogao_seqinfo$chipsampleid[s,2],"_",var1,".bam.bai")
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

      percentCovered <- percent((coveredBases / genomeLength),accuracy = 0.01)   # Percent covered
      totalReads <- length(bam)  # Total reads
      
      # Store details and stats
      selectedpathogen_data[selectedpathogen_data$拉丁名 ==var1,"覆盖度"] <- percentCovered
    } 
  }
  
  #----------------------------------------------------------
  # Section 2、Summary pathogen test results
  
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
  

  # Adjust the table for better fit within the Word document

  # Adjust the flextable width to 16.38 cm
  ft_z <- autofit(ft_z) # Fit content first
  ft_z <- flextable::width(ft_z, width = rep(16.38/2.54/(ncol(selectedpathogen_data)-6), (ncol(selectedpathogen_data)-6))) 
  
  
  # Set text alignment to centre for all columns
  ft_z <- align(ft_z, align = "left", part = "all") 
  

  #----------------------------------------------------------
  # Section 3 List of test result.
  
  #-------------------
  # 3.1 Bacteria list
  # 3.1.1 Bacteria
  selected_bacteria <- selectedpathogen_data[which(selectedpathogen_data$Kingdom %in% "Bacteria"),c("菌属中文名","菌种中文名","类型","拉丁名","序列数","覆盖度")]
  selected_bacteria_b <- selected_bacteria[which(selected_bacteria$菌属中文名 != "分支杆菌属"),]
  colnames(selected_bacteria_b) <- c("Genus_Chinese","Species_Chinese","Types","Species_latin","Reads","Coverage")
  print(selected_bacteria_b)
  #print(dim(selected_bacteria_b))
  ft_b <- flextable(selected_bacteria_b, col_keys = c("Species","Types","Reads","Coverage"))%>%
    mk_par(j = "Species", value = as_paragraph(Species_Chinese, as_i(Species_latin))) %>%
    set_header_labels(Species = "Species")
  

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

  # Adjust the flextable width to 16.38 cm
  ft_b <- autofit(ft_b) # Fit content first
  ft_b <- flextable::width(ft_b, width = rep(16.38/2.54/(ncol(selectedpathogen_data)-3), (ncol(selectedpathogen_data)-3))) # Divide total width across columns
  # ncol(selectedpathogen_data)-1 =5, ft_b also have 5 columns
  
  # Set text alignment to centre for all columns
  ft_b <- align(ft_b, align = "center", part = "all") 
  #---
  # 3.1.2 Mycobacterium
  selected_bacteria_mycobacterium <- selected_bacteria[which(selected_bacteria$菌属中文名 %in% "分支杆菌属"),]
  print(selected_bacteria_mycobacterium)
  colnames(selected_bacteria_mycobacterium) <- c("Genus_Chinese","Species_Chinese","Types","Species_latin","Reads","Coverage")

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
  # 3.2 Eukaryota
  # 3.2.1 Fungi
  selected_euk <- selectedpathogen_data[which(selectedpathogen_data$Kingdom %in% "Eukaryota"),c("菌属中文名","菌种中文名","类型","拉丁名","序列数","覆盖度")]
  selected_euk_euk <- selected_euk[which(selected_euk$类型 != "par"),]
  colnames(selected_euk_euk) <- c("Genus_Chinese","Species_Chinese","Types","Species_latin","Reads","Coverage")

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
  # 3.2.2 Par
  selected_euk_par <- selected_euk[which(selected_euk$类型 %in% "par"),]
  
  colnames(selected_euk_par) <- c("Genus_Chinese","Species_Chinese","Types","Species_latin","Reads","Coverage")

  
  #print(dim(selected_bacteria_b))
  ft_p <- flextable(selected_euk_par, col_keys = c("Species","Types","Reads","Coverage"))%>%
    mk_par(j = "Species", value = as_paragraph(Species_Chinese, as_i(Species_latin))) %>%
    set_header_labels(Species = "Species")

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
  # 3.3 Viruses list
  selected_vir <- selectedpathogen_data[which(selectedpathogen_data$Kingdom %in% "Viruses"),c("菌属中文名","菌种中文名","类型","拉丁名","序列数","覆盖度")]
  colnames(selected_vir) <- c("Genus_Chinese","Species_Chinese","Types","Species_latin","Reads","Coverage")

  
  #print(dim(selected_bacteria_b))
  ft_v <- flextable(selected_vir, col_keys = c("Species","Types","Reads","Coverage"))%>%
    mk_par(j = "Species", value = as_paragraph(Species_Chinese, as_i(Species_latin))) %>%
    set_header_labels(Species = "Species")
  
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
  # 3.4 Potential pathogen list

  d_pathogens <- pathogeninfo[which(pathogeninfo[,"species"] %in% pathogen_human[,1]),]
  d_pathogens <- d_pathogens[-which(d_pathogens[,"species"] %in% pathogen_data[,"species"]),]
  d_pathogens <- d_pathogens[,c("菌属中文名","菌种中文名","Gram+/-","species","Reads")] 

  colnames(d_pathogens) <- c("Genus_Chinese","Species_Chinese","Types","Species_latin","Reads")
  
  ft_y <- flextable(d_pathogens, col_keys = c("Species","Types","Reads"))%>%
    mk_par(j = "Species", value = as_paragraph(Species_Chinese, as_i(Species_latin))) %>%
    set_header_labels(Species = "Species")

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
  # Section 4、Sequencing quality
  
  #---
  # 4.1 Quality control
  sequencinginfo <- as.data.frame(fread(paste0("RpNGS/data/response/RpNGS_sequencinginfo.csv")),stringsAsFactors=F)
  onesequencinginfo <- sequencinginfo[which(sequencinginfo$样本编号 %in% baogao_seqinfo$chipsampleid[s,2]),]
  doc <- body_replace_all_text(doc, "HSND", as.character(onesequencinginfo[1,3]))
  doc <- body_replace_all_text(doc, "WKND", as.character(onesequencinginfo[1,5]))
  
  #---
  # 4.2 Data control
  #---
  
  #--
  # fastp file 
  html_file<- paste0(path4results,trimws(input$reprot_chipid, which ="both"),"/01fastp/",baogao_seqinfo$chipsampleid[s,2],".html")
  
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
  fastq_file <- paste0(path4results,trimws(input$reprot_chipid, which ="both"),"/02decomt/",baogao_seqinfo$chipsampleid[s,2],".fastq.gz")
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