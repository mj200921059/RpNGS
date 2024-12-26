


# 0. sequencer@ip address 
seqer_ip <- "NWU_M@192.168.14.50" # SHOULD CHANGE BASED ON THE SSH CONNECTION 

# 1. Set path
path4datasets <- paste0(getwd(),"/data/response/")
path4results <- paste0(getwd(),"/data/results/")
path4template <- paste0(getwd(),"/data/word_template/")
path4sequencer <- "C:/Result/OutputFq/" # The path of fastq files in BGI G99 sequencer 

#****
seqinfo_file <- paste0(path4datasets,"RpNGS_sequencinginfo.csv")
clininfo_file <- paste0(path4datasets,"RpNGS_clinicalinfo.csv")
help_file <- paste0(getwd(),"/www/html/help.Rmd")
title_file <-  paste0(getwd(),"/www/title.Rmd") 
footer_file <- paste0(getwd(),"/www/html/footer_v4.html") 
map_file <-  paste0(getwd(),"/image/Shaanxi.json") 

#-------------------------------------------------------------------------------
# 2. Input example data in analysis sets panel

Input <- structure(list(Chip_id = c("FT10005323"), Sample_id = c("lifegen20240830"), Extracted_NAC = c(
  "2.74"), Adaptor = c("126"),library_NAC = c("58.6"), Rawfastq_id = c(
    "FT10005323_L10_126.fq.gz")), class = c(
      "spec_tbl_df", "tbl_df", "tbl",
      "data.frame"
    ), row.names = c(NA, -1L))




