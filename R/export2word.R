replace_text <- function(doc, placeholder, value) {
  body_replace_all_text(doc, placeholder, as.character(value))
}

create_word_document <- function(pathogen_data, pathogeninfo, pathogen_human, output_file) {
  doc <- read_docx("clinicalreporttemplate.docx")
  
  # 读取临床信息
  clinicalinfo <- as.data.frame(fread("clininfo_file.csv"))
  s <- input$sampleid_list_rows_selected
  sample_id <- baogao_seqinfo$chipsampleid[s,2]
  onelineclinical_df <- clinicalinfo[clinicalinfo$Sample_id == sample_id, ]
  
  # 替换临床信息
  fields <- c("YBBH", "xingming", "xingbie", "nianling", "CYRQ", "jieyangriqi", "dianhua")
  values <- onelineclinical_df[1, 1:7]
  
  for (i in seq_along(fields)) {
    doc <- replace_text(doc, fields[i], values[i])
  }
  
  # 计算 BAM 覆盖度（优化计算）
  get_coverage <- function(var1) {
    bamPath <- glue("{path4results}/{input$reprot_chipid}/05alignment/{sample_id}_{var1}.bam")
    baiPath <- glue("{bamPath}.bai")
    
    if (file.exists(bamPath) && file.exists(baiPath)) {
      bam <- readGAlignments(bamPath)
      coverage_vec <- as.numeric(coverage(bam)[[1]])
      percentCovered <- sum(coverage_vec > 0) / length(coverage_vec) * 100
      return(glue("{round(percentCovered, 2)}%"))
    } else {
      return("N/A")
    }
  }
  
  selectedpathogen_data$覆盖度 <- sapply(selectedpathogen_data$拉丁名, get_coverage)
  
  # 生成 Word 文档
  doc <- body_add_flextable(doc, flextable(selectedpathogen_data))
  print(doc, target = output_file)
}