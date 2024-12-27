#' mircrobes identification based on fastq file analysis 
#' User should create a conda env and install fastp, bowtie2, kraken2,bracken.
#' 
#' 
#' 
#' @param distdir directory where all the results will be saved.
#' @param fastq_idr directory of the fastq files.
#' @param n_thread number of cores to use.
#' 
#' @return csv files in 04rawresults folder.
# 
#' @example 
#' 
#'  

# run_pngsanalysis("/home/dell/dataanalysis/projects/web/fastq/20240326","/home/dell/dataanalysis/projects/mngs/20240326/fastq",48)




run_pngsanalysis <- function(destdir,fastq_dir,n_thread){
  library(dplyr)
  library(data.table)
  # ----------------------------------------------------------------------------
  # STEP 1 REMOVING AND FILTERING
  cat(paste("Running Preprocessing ...",Sys.time(),"\n",sep = ""))
  fastq_files = list.files(fastq_dir,pattern = ".fq.gz$",full.names = T)
  #fastq_files = list.files(fastq_dir,pattern = ".fastq.gz$",full.names = T)
  setwd(destdir)
  if(!dir.exists("01fastp")){
    system(paste0("mkdir ",destdir,"/01fastp"))
  }

  for (n in fastq_files) {
    fn <- paste0(destdir,"/01fastp/",last(unlist(strsplit(n, split='/', fixed=TRUE))))
    comn <- last(unlist(strsplit(n, split='/', fixed=TRUE)))
    jn <- paste0(destdir,"/01fastp/",unlist(strsplit(comn, split='.', fixed=TRUE))[1],".json")
    hn <- paste0(destdir,"/01fastp/",unlist(strsplit(comn, split='.', fixed=TRUE))[1],".html")
    system(paste0("conda run -n fastp fastp -i ",n," -o ",fn," -w ",n_thread," -5"," -3"," -W 4"," -p ","--length_required 36"," -q 15 -u 40"," -y "," -Y 30"," -D"," -j ",jn," -h ",hn))
  }
  # ----------------------------------------------------------------------------
  # STEP 2 DECOMT
  cat(paste("Running reads cleaning ...",Sys.time(),"\n",sep = ""))
  fastp_dir <-  paste0(destdir,"/01fastp")
  fastq_files = list.files(fastp_dir,pattern = ".fq.gz$",full.names = T)

  if(!dir.exists("02decomt")){
    system(paste0("mkdir ",destdir,"/02decomt"))
  }

  for (n in fastq_files) {
    comn <- last(unlist(strsplit(n, split='/', fixed=TRUE)))
    fn <- paste0(destdir,"/02decomt/",comn)
    sn <- paste0(destdir,"/02decomt/",unlist(strsplit(comn, split='.', fixed=TRUE))[1],".sam")
    system(paste0("conda run -n fastp bowtie2 -p ",n_thread," -q -x /home/dell/dataanalysis/pipelines/databases/bowtie2_index/hg37dec_v0.1 -U ",n ," --very-sensitive --un-gz ",fn," -S ",sn))
    system(paste0('rm ',destdir,"/02decomt/*sam"))
  }
  # ----------------------------------------------------------------------------
  # STEP 3 CLASSIFICATION
  cat(paste("Running reads classification ...",Sys.time(),"\n",sep = ""))
  fastp_dir <-  paste0(destdir,"/02decomt")
  fastq_files = list.files(fastp_dir,pattern = ".fq.gz$",full.names = T)
  if(!dir.exists("03kraken")){
    system(paste0("mkdir ",destdir,"/03kraken"))
  }

  for (n in fastq_files) {
    comn <- last(unlist(strsplit(n, split='/', fixed=TRUE)))
    on <- paste0(destdir,"/03kraken/",unlist(strsplit(comn, split='.', fixed=TRUE))[1],"_k2output.txt")
    rn <- paste0(destdir,"/03kraken/",unlist(strsplit(comn, split='.', fixed=TRUE))[1],"_k2report.txt")
    bn <- paste0(destdir,"/03kraken/",unlist(strsplit(comn, split='.', fixed=TRUE))[1],"_brackentaxa.tsv")
    #system(paste0("conda run -n fastp kraken2 --db /home/dell/dataanalysis/pipelines/databases/karken_database_20220607 --threads ",n_thread," --use-names --minimum-hit-groups 3 --report-minimizer-data ",n," --output ",on," --report ",rn))
    #system(paste0("conda run -n fastp kraken2 --db /home/dell/dataanalysis/pipelines/databases/h_bavfp_k2db/pngsk2db20240220 --threads ",n_thread," --use-names --minimum-hit-groups 3 --report-minimizer-data ",n," --output ",on," --report ",rn))
    system(paste0("conda run -n fastp kraken2 --db /home/dell/dataanalysis/pipelines/databases/h_bavfp_k2db/pngsk2db20240418 --threads ",n_thread," --use-names --minimum-hit-groups 3 --report-minimizer-data ",n," --output ",on," --report ",rn))
    #system(paste0("conda run -n fastp bracken -d /home/dell/dataanalysis/pipelines/databases/karken_database_20220607 "," -i ",rn," -o ",bn," -r 50 -l S"))
    #system(paste0("conda run -n fastp bracken -d /home/dell/dataanalysis/pipelines/databases/h_bavfp_k2db/pngsk2db20240220 "," -i ",rn," -o ",bn," -r 50 -l S"))
    system(paste0("conda run -n fastp bracken -d /home/dell/dataanalysis/pipelines/databases/h_bavfp_k2db/pngsk2db20240418 "," -i ",rn," -o ",bn," -r 50 -l S"))
  }
  # ----------------------------------------------------------------------------
  # STEP 4 Raw results with annotation 
  cat(paste("Running Annotation ...",Sys.time(),"\n",sep = ""))
  taxo_dir <-  paste0(destdir,"/03kraken")
  taxo_files = list.files(taxo_dir,pattern = ".tsv$",full.names = T)

  if(!dir.exists("04rawresults")){
    system(paste0("mkdir ",destdir,"/04rawresults"))
  }
  load('~/dataanalysis/projects/lifegen/pre_datasets/microbes_taxonomy.rda')
  microbes_taxonomy <- microbes_taxonomy[,-1]
#-----------------------  
  # GET FILE NAMES
  folder_files <- list.files(path = "./03kraken",pattern = "tsv" ,full.names = T, recursive = FALSE)
  
  filenames = ""
  for (name in folder_files) {
    filenames <- c(filenames, unlist(strsplit(name, split='/', fixed=TRUE))[3])}
  filenames <- filenames[-1]
  
  
  sample_names = ""
  for (nam in filenames) {
    sample_names <- c(sample_names, unlist(strsplit(nam, split='_brackentaxa', fixed=TRUE))[1])}
  sample_names <- sample_names[-1]
  sample_names <- unique(sample_names)
  
  # Create files with a loop
  for (x in sample_names) {
    # 01 load report data for DNA virues 
    taxa_report <- paste("./03kraken/",x,"_k2report",".txt", sep="")
    taxa_report_abun <- fread(taxa_report) 
    #taxa_report_abun <- read.table(file = x, header = T, quote = "", sep = "\t")
    colnames(taxa_report_abun) <- c("Perct","Reads","Specific_reads","Kmers","Specific Kmers","Level",'taxonomy_id','Taxa_names')
    
    taxa_report_abun <- taxa_report_abun[which(taxa_report_abun$Level %in% c("S","S1","S2","S3","S4")),]
    taxa_report_abun <- taxa_report_abun[!taxa_report_abun$Taxa_names %in% "Homo sapiens",]
    taxa_report_abun <- merge(taxa_report_abun,microbes_taxonomy,by='taxonomy_id')
    taxa_report_abun <- taxa_report_abun[taxa_report_abun$Kingdom %in% "Viruses",]
    #taxa_report_abun <- taxa_report_abun %>% group_by(species) %>% slice_min(V6, n=1) %>% ungroup() 
    taxa_report_abun <- taxa_report_abun %>% group_by(Kingdom) %>% mutate(Perct = Specific_reads/sum(Specific_reads)) %>% ungroup() 
    taxa_report_abun <- taxa_report_abun[,-5:-6]
    
    # 02 Load bracken data for microbes without viruses
    taxa_bracken <- paste("./03kraken/",x,"_k2report","_bracken_species.txt",sep="")
    taxa_bracken_abun <- fread(taxa_bracken) 
    #taxa_report_abun <- read.table(file = x, header = T, quote = "", sep = "\t")
    taxa_bracken_abun <- taxa_bracken_abun[which(taxa_bracken_abun$V4 %in% c("S","S1","S2","S3","S4")),]
    taxa_bracken_abun <- taxa_bracken_abun[!taxa_bracken_abun$V6 %in% "Homo sapiens",]
    taxa_bracken_abun <- merge(taxa_bracken_abun,microbes_taxonomy ,by.x = "V5", by.y = "taxonomy_id")
    taxa_bracken_abun <- taxa_bracken_abun[taxa_bracken_abun$Kingdom != "Viruses",]
    #taxa_bracken_abun <- taxa_bracken_abun %>% group_by(species) %>% slice_min(V4, n=1) %>% ungroup() 
    taxa_bracken_abun <- taxa_bracken_abun %>%
      group_by(Kingdom) %>%
      mutate(V1 = V3/sum(V3)) %>% ungroup() 
    colnames(taxa_bracken_abun)[1:6] <- c('taxonomy_id',"Perct","Reads","Specific_reads","Level",'Taxa_names')
    
    # 03 merge two tables 
    taxa_bracken_abun <- rbind(taxa_bracken_abun,taxa_report_abun)
    rn<- unlist(strsplit(x, split='_k2report', fixed=TRUE))[1]
    write.csv(taxa_bracken_abun,paste("./04rawresults/",rn,".csv", sep=""))
  }
  
}





