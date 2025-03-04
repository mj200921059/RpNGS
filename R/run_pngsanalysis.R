run_pngsanalysis <- function(destdir, fastq_dir, n_thread) {
  # library(dplyr)
  # library(data.table)

  log_text(paste0(log_text(),"\nFor sample file ",p,"\nMetagenomic analysis start at ", Sys.time())) # Update Log
  # 读取所有 fastq 文件
  fastq_files <- list.files(fastq_dir, pattern = "\\.fq.gz$", full.names = TRUE)

  # Step 1: 质量控制 (fastp)
  run_fastp <- function(file, destdir, n_thread) {
    create_dir(file.path(destdir, "01fastp"))
    output_file <- file.path(destdir, "01fastp", basename(file))
    json_file <- sub("\\.fq.gz$", ".json", output_file)
    html_file <- sub("\\.fq.gz$", ".html", output_file)
    #output$fpp <- renderText(paste0("File Processing Progress:",p))
    cmd <- sprintf(
      "/home/majun/miniconda3/condabin/conda run -n fastp fastp -i %s -o %s -w %d -q 15 -l 36 --json %s --html %s",
      file, output_file, n_thread, json_file, html_file
    ) # reminder: using which conda to find the path of conda in the system and replace the conda path
    
    tryCatch({
      oput <-system(cmd, intern = TRUE) # Capture output
      log_text(paste0(log_text(),paste(oput, collapse = "\n"))) # Append outpu to Log
    } , error = function(e) {
      #message("Error in fastp: ", e$message)
      log_text(paste0(log_text(),"\nError:", e$message, "\n"))
    })
    
  }
  
  #message("Running Fastp for quality control...")
  log_text(paste0(log_text(),"\nRunning Fastp for quality control...")) # Update Log
  # run QC fun
  lapply(fastq_files, run_fastp, destdir, n_thread)
  # run_fastp("./results/S1/fq/FT10002025_L0_12.fq.gz","./results/S1",6)

  

  
  
  # Step 2: 去宿主序列 (Bowtie2)
  run_bowtie2 <- function(file, destdir, n_thread) {
    qc_file <- file.path(destdir, "01fastp", basename(file))
    create_dir(file.path(destdir, "02decomt"))
    output_file <- file.path(destdir, "02decomt", basename(file))
    sam_file <- sub("\\.fq.gz$", ".sam", output_file)

    cmd <- sprintf(
      "/home/majun/miniconda3/condabin/conda run -n fastp bowtie2 -p %d -x /home/majun/data_analysis/pipelines/databases/bowtie2/GRCh37/GRCh37 -U %s --very-sensitive --un-gz %s -S %s",
      n_thread, qc_file, output_file, sam_file
    )

    tryCatch(system(cmd, intern = TRUE), error = function(e) {
      message("Error in bowtie2: ", e$message)
    })

    unlink(sam_file) # 删除 SAM 文件
  }

  #message("Running Bowtie2 for host read removal...")
  log_text(paste0(log_text(),"\nRunning Bowtie2 for host read removal...")) # Update Log
  lapply(fastq_files, run_bowtie2, destdir, n_thread)
  
  # Step 3: 物种分类 (Kraken2 & Bracken)
  run_kraken <- function(file, destdir, n_thread) {
    create_dir(file.path(destdir, "03kraken"))
    removed_file <- file.path(destdir, "02decomt", basename(file))
    base_name <- sub("\\.fq.gz$", "", basename(file))
    output_file <- file.path(destdir, "03kraken", paste0(base_name, "_k2output.txt"))
    report_file <- file.path(destdir, "03kraken", paste0(base_name, "_k2report.txt"))
    bracken_file <- file.path(destdir, "03kraken", paste0(base_name, "_brackentaxa.tsv"))

    cmd_kraken <- sprintf(
      "/home/majun/miniconda3/condabin/conda run -n fastp kraken2 --db /home/dell/dataanalysis/pipelines/databases/h_bavfp_k2db/pngsk2db20240418 --threads %d --use-names --minimum-hit-groups 3 --report-minimizer-data %s --output %s --report %s",
      n_thread, removed_file, output_file, report_file
    )

    cmd_bracken <- sprintf(
      "/home/majun/miniconda3/condabin/conda run -n fastp bracken -d /home/dell/dataanalysis/pipelines/databases/h_bavfp_k2db/pngsk2db20240418 -i %s -o %s -r 50 -l S",
      report_file, bracken_file
    )

    tryCatch(system(cmd_kraken, intern = TRUE), error = function(e) {
      message("Error in Kraken2: ", e$message)
    })
    tryCatch(system(cmd_bracken, intern = TRUE), error = function(e) {
      message("Error in Bracken: ", e$message)
    })
  }

  #message("Running Kraken2 for taxonomy classification...")
  log_text(paste0(log_text(),"\nRunning Kraken2 for taxonomy classification...")) # Update Log
  lapply(fastq_files, run_kraken, destdir, n_thread)
  # 
  # Step 4: Raw results with annotation
  run_annot <- function(file, destdir) {
      create_dir(file.path(destdir, "04rawresults"))
      base_name <- sub("\\.fq.gz$", "", basename(file))
      output_file <- file.path(destdir, "04rawresults", paste0(base_name, "_k2report.csv"))
      report_file <- file.path(destdir, "03kraken", paste0(base_name, "_k2report.txt"))
      bracken_file <- file.path(destdir, "03kraken", paste0(base_name, "_brackentaxa.tsv"))

      # 01 load report data for DNA virues
      taxa_report_abun <- fread(report_file)
      colnames(taxa_report_abun) <- c("Perct","Reads","Specific_reads","Kmers","Specific Kmers","Level",'taxonomy_id','Taxa_names')
      taxa_report_abun <- taxa_report_abun[which(taxa_report_abun$Level %in% c("S","S1","S2","S3","S4")),]
      taxa_report_abun <- taxa_report_abun[!taxa_report_abun$Taxa_names %in% "Homo sapiens",]
      taxa_report_abun <- merge(taxa_report_abun,microbes_taxonomy,by='taxonomy_id')
      taxa_report_abun <- taxa_report_abun[taxa_report_abun$Kingdom %in% "Viruses",]
      taxa_report_abun <- taxa_report_abun %>% group_by(Kingdom) %>% mutate(Perct = Specific_reads/sum(Specific_reads)) %>% ungroup()
      taxa_report_abun <- taxa_report_abun[,-5:-6]

      # 02 Load bracken data for microbes without viruses
      taxa_bracken_abun <- fread(bracken_file)
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
      write.csv(taxa_bracken_abun,output_file,row.names = FALSE)



    }

  #message("Running Kraken2 for taxonomy annotation...")
  log_text(paste0(log_text(),"\nLoading microbe.db for taxonomy annotation...")) # Update Log
  lapply(fastq_files, run_annot, destdir)
  
  #message("Metagenomic analysis completed at ", Sys.time())
  log_text(paste0(log_text(),"\nMetagenomic analysis completed at ", Sys.time(), "\n")) # Update Log
 
  return(log_text())
}
