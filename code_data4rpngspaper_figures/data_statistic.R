# For regenerate the Pathogens identification pipeline evaluation figure

### NOTE: Given the data size is large. We upload the tar file into github. 
# So uncompressed the idcz_Sample_Taxo_ Reports.tar file into "./idcz_Sample_Taxo_Reports/" before run this code. 

# ---------------------------
# Calculate P, R, and F1 of IDseq nt, nr and combo nt+nr. and RpNGS

# 1.1 WITHOUT Z-SCORE CORRECT
# Load the package
library(readxl)
library(stringr)
library(readr)
#------------
# Common data 
# Read the Excel file
tb1 <- read_excel("table_S1.xlsx", sheet = 1)
tb1 <- tb1[-24:-25,]
#tb1_dna <- tb1[-which(tb1$Type %in% "RNA virus"),]
#---
tb2 <- read_excel("table_S2.xlsx", sheet = 1)
# Extract "DNA" and "RNA"
tb2$type <- str_extract(tb2$Experiment_title, "(DNA|RNA)")

tb1idseq <- read_excel("table_S1idseq.xlsx", sheet = 1)
tb1idseq <- tb1idseq[-24:-25,]
#----------
# IDseq

# Directory where the CSV files are stored
input_dir <- "./idcz_Sample_Taxo_Reports/rpngs_18874"  


# Create an empty data table  to store the data frames
all_statdt <- data.frame()

# Loop over each row of the tibble
for (i in 1:nrow(tb2)) {
  # Construct a pattern to search for the file (e.g., HRR1209919_*_taxon_report.csv)
  #pattern <- paste0("^", tb2$Accession[i], "_[0-9]+_taxon_report\\.csv$")
  pattern <- paste0("^", tb2$Accession[i], "_.*_taxon_report\\.csv$")
  # Get the list of files matching the pattern
  matching_files <- list.files(input_dir, pattern = pattern, full.names = TRUE)
  
  # If a matching file is found, read it
  if (length(matching_files) > 0) {
    # Read the first matched file (assuming one file per Accession)
    report <- read_csv(matching_files[1])
    report <- report[which(report$tax_level %in% 1),]
    
    if (tb2$type[i] %in% "DNA") {
      tb1_groundtruth_dnatax <- tb1idseq[-which(tb1idseq$Type %in% "RNA virus"),c("Microbes",tb2$Sample_Name[i])]
      colnames(tb1_groundtruth_dnatax) <- c("Microbes","value")
      tb1_groundtruth_dnatax <- tb1_groundtruth_dnatax[which(tb1_groundtruth_dnatax$value != "\\"),]
      report_nt <- report[which(report$nt_rpm >0),]

      TP <- length(intersect(report_nt$name, tb1_groundtruth_dnatax$Microbes))
      FP <- length(report_nt$tax_level) - length(intersect(report_nt$name, tb1_groundtruth_dnatax$Microbes))
      
      FN <- length(tb1_groundtruth_dnatax$Microbes) - length(intersect(report_nt$name, tb1_groundtruth_dnatax$Microbes))
      P <- TP/(TP+FP)
      R <- TP/(TP+FN)
      F1=2/((1/R) + (1/P))
      rpm_type <- "NT"
      statdt_dna <- data.frame(tb2$Sample_Name[i],TP,FP,FN,P,R,F1,tb2$type[i],tb2$Accession[i],rpm_type)
      # Store the data frame in the list
      all_statdt<- rbind(all_statdt,statdt_dna )
      
      report_nr <- report[which(report$nr_rpm >0),]
      
      TP <- length(intersect(report_nr$name, tb1_groundtruth_dnatax$Microbes))
      FP <- length(report_nr$tax_level) - length(intersect(report_nr$name, tb1_groundtruth_dnatax$Microbes))
      
      FN <- length(tb1_groundtruth_dnatax$Microbes) - length(intersect(report_nr$name, tb1_groundtruth_dnatax$Microbes))
      P <- TP/(TP+FP)
      R <- TP/(TP+FN)
      F1=2/((1/R) + (1/P))
      rpm_type <- "NR"
      statdt_dna <- data.frame(tb2$Sample_Name[i],TP,FP,FN,P,R,F1,tb2$type[i],tb2$Accession[i],rpm_type)
      # Store the data frame in the list
      all_statdt<- rbind(all_statdt,statdt_dna )
      
      
      
    } else {
      tb1_groundtruth_rnatax <- tb1idseq[which(tb1idseq$Type %in% "RNA virus"),c("Microbes",tb2$Sample_Name[i])]
      colnames(tb1_groundtruth_rnatax) <- c("Microbes","value")
      tb1_groundtruth_rnatax <- tb1_groundtruth_rnatax[which(tb1_groundtruth_rnatax$value != "\\"),]
      report_nt <- report[which(report$nt_rpm >0),]
      
      TP <- length(intersect(report_nt$name, tb1_groundtruth_rnatax$Microbes))
      FP <- length(report_nt$tax_level) - length(intersect(report_nt$name, tb1_groundtruth_rnatax$Microbes))
      
      FN <- length(tb1_groundtruth_rnatax$Microbes) - length(intersect(report_nt$name, tb1_groundtruth_rnatax$Microbes))
      P <- TP/(TP+FP)
      R <- TP/(TP+FN)
      F1=2/((1/R) + (1/P))
      rpm_type <- "NT"
      statdt_rna <- data.frame(tb2$Sample_Name[i],TP,FP,FN,P,R,F1,tb2$type[i],tb2$Accession[i],rpm_type)
      all_statdt<- rbind(all_statdt,statdt_rna )
      
      report_nr <- report[which(report$nr_rpm >0),]
      
      TP <- length(intersect(report_nr$name, tb1_groundtruth_rnatax$Microbes))
      FP <- length(report_nr$tax_level) - length(intersect(report_nr$name, tb1_groundtruth_rnatax$Microbes))
      
      FN <- length(tb1_groundtruth_rnatax$Microbes) - length(intersect(report_nr$name, tb1_groundtruth_rnatax$Microbes))
      P <- TP/(TP+FP)
      R <- TP/(TP+FN)
      F1=2/((1/R) + (1/P))
      rpm_type <- "NR"
      statdt_rna <- data.frame(tb2$Sample_Name[i],TP,FP,FN,P,R,F1,tb2$type[i],tb2$Accession[i],rpm_type)
      all_statdt<- rbind(all_statdt,statdt_rna )
      
      
    }


  } else {
    message(paste("No matching file for:", tb2$Accession[i]))
  }
}

# 
all_statdt$platform <- "IDseq" 
save(all_statdt, file = "idseq_all_statdt.rda")

#----------
# RpNGS

# Directory where the CSV files are stored
input_dir <- "./pngs/04rawresults"


# Create an empty data table  to store the data frames
all_statdt_pngs <- data.frame()

# Loop over each row of the tibble
for (i in 1:nrow(tb2)) {
  # Construct a pattern to search for the file (e.g., HRR1209919_*_taxon_report.csv)
  #pattern <- paste0("^", tb2$Accession[i], "_.*_taxon_report\\.csv$")
  pattern <- paste0(tb2$Accession[i], ".csv")
  # Get the list of files matching the pattern
  matching_files <- list.files(input_dir, pattern = pattern, full.names = TRUE)
  
  # If a matching file is found, read it
  if (length(matching_files) > 0) {
    # Read the first matched file (assuming one file per Accession)
    report <- read_csv(matching_files[1])
    #report <- report[which(report$tax_level %in% 1),]
    
    if (tb2$type[i] %in% "DNA") {
      tb1_groundtruth_dnatax <- tb1[-which(tb1$Type %in% "RNA virus"),c("Microbes",tb2$Sample_Name[i])]
      colnames(tb1_groundtruth_dnatax) <- c("Microbes","value")
      tb1_groundtruth_dnatax <- tb1_groundtruth_dnatax[which(tb1_groundtruth_dnatax$value != "\\"),]

      TP <- length(unique(intersect(report$Taxa_names, tb1_groundtruth_dnatax$Microbes)))
      FP <- length(unique(report$Taxa_names)) - length(unique(intersect(report$Taxa_names, tb1_groundtruth_dnatax$Microbes)))
      
      FN <- length(tb1_groundtruth_dnatax$Microbes) - length(unique(intersect(report$Taxa_names, tb1_groundtruth_dnatax$Microbes)))
      P <- TP/(TP+FP)
      R <- TP/(TP+FN)
      F1=2/((1/R) + (1/P))
      rpm_type <- "kmer"
      statdt_dna <- data.frame(tb2$Sample_Name[i],TP,FP,FN,P,R,F1,tb2$type[i],tb2$Accession[i],rpm_type)
      # Store the data frame in the list
      all_statdt_pngs<- rbind(all_statdt_pngs,statdt_dna )
      
    } else {
      tb1_groundtruth_rnatax <- tb1[which(tb1$Type %in% "RNA virus"),c("Microbes",tb2$Sample_Name[i])]
      colnames(tb1_groundtruth_rnatax) <- c("Microbes","value")
      tb1_groundtruth_rnatax <- tb1_groundtruth_rnatax[which(tb1_groundtruth_rnatax$value != "\\"),]

      TP <- length(unique(intersect(report$Taxa_names, tb1_groundtruth_rnatax$Microbes)))
      FP <- length(unique(report$Taxa_names)) - length(unique(intersect(report$Taxa_names, tb1_groundtruth_rnatax$Microbes)))
      
      FN <- length(tb1_groundtruth_rnatax$Microbes) - length(unique(intersect(report$Taxa_names, tb1_groundtruth_rnatax$Microbes)))
      P <- TP/(TP+FP)
      R <- TP/(TP+FN)
      F1=2/((1/R) + (1/P))
      rpm_type <- "kmer"
      statdt_rna <- data.frame(tb2$Sample_Name[i],TP,FP,FN,P,R,F1,tb2$type[i],tb2$Accession[i],rpm_type)
      all_statdt_pngs<- rbind(all_statdt_pngs,statdt_rna )
      
    }
    
    
  } else {
    message(paste("No matching file for:", tb2$Accession[i]))
  }
}

all_statdt_pngs$platform <- "RpNGS"
save(all_statdt_pngs, file = "rpngs_all_statdt.rda")




#==============================
# merge two tables and make figures
load("rpngs_all_statdt.rda")
load("idseq_all_statdt.rda")

all_statdt <- rbind(all_statdt_pngs,all_statdt)
write.csv(all_statdt, "allstatdt.csv")

#=====================
df <- read.csv("allstatdt.csv")

# Load necessary libraries
library(ggplot2)
library(dplyr)

### F-scores Plot (using F1 column) ###
df_F1 <- df %>%
  select(tb2.Sample_Name.i., tb2.type.i., rpm_type, platform, F1) %>%
  rename(score = F1)

plot_F1 <- ggplot(df_F1, aes(x = tb2.Sample_Name.i., y = score)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  facet_grid(tb2.type.i. ~ rpm_type + platform) +
  labs(
    title = "F-scores",
    x = "Sample Name",
    y = "F-scores"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

### Recall Plot (using R column) ###
df_R <- df %>%
  select(tb2.Sample_Name.i., tb2.type.i., rpm_type, platform, R) %>%
  rename(score = R)

plot_R <- ggplot(df_R, aes(x = tb2.Sample_Name.i., y = score)) +
  geom_bar(stat = "identity", fill = "darkgreen") +
  facet_grid(tb2.type.i. ~ rpm_type + platform) +
  labs(
    title = "Recall",
    x = "Sample Name",
    y = "Recall"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

### Precision Plot (using P column) ###
df_P <- df %>%
  select(tb2.Sample_Name.i., tb2.type.i., rpm_type, platform, P) %>%
  rename(score = P)

plot_P <- ggplot(df_P, aes(x = tb2.Sample_Name.i., y = score)) +
  geom_bar(stat = "identity", fill = "firebrick") +
  facet_grid(tb2.type.i. ~ rpm_type + platform) +
  labs(
    title = "Precision",
    x = "Sample Name",
    y = "Precision"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# To view the plots separately:
pdf(file = "fscore.pdf",width = 12, height = 8,)
plot(plot_F1)
dev.off()

pdf(file = "r.pdf",width = 12, height = 8,)
plot(plot_R)
dev.off()

pdf(file = "p.pdf",width = 12, height = 8,)
plot(plot_P)
dev.off()




#-------------------------------------------------------------------------------
# 1.2 WITH Z-SCORE CORRECT

# *****1.2.1 load idseq csv file to combine all abundance tables into one data frame *****#
# Directory where the CSV files are stored
input_dir <- "./idcz_Sample_Taxo_Reports/rpngs_18874"
# Create an empty list to store the data frames
all_reports_dna <- list()
all_reports_rna <- list()
# Loop over each row of the tibble
for (i in 1:nrow(tb2)) {
  # Construct a pattern to search for the file (e.g., HRR1209919_*_taxon_report.csv)
  pattern <- paste0("^", tb2$Accession[i], "_.*_taxon_report\\.csv$")
  # Get the list of files matching the pattern
  matching_files <- list.files(input_dir, pattern = pattern, full.names = TRUE)
  
  # If a matching file is found, read it
  if (length(matching_files) > 0 ) {
    # Read the first matched file (assuming one file per Accession)
    report <- read_csv(matching_files[1])
    report <- report[which(report$tax_level %in% 1),]
    if(tb2$type[i] %in% "DNA"){
      all_reports_dna[[tb2$Accession[i]]] <- report
    } else {
      all_reports_rna[[tb2$Accession[i]]] <- report
    }
    

  } else {
    message(paste("No matching file for:", tb2$Accession[i]))
  }
}


# 
save(all_reports_dna ,file = "all_reports_dna.rda")
save(all_reports_rna ,file = "all_reports_rna.rda")


#------------------------------
# 1.2.2  Calculate z scores FOR IDseq DNA SAMPLES 
load("all_reports_dna.rda")
all_reports<- all_reports_dna
combined_df <- bind_rows(
  lapply(names(all_reports), function(sample) {
    df <- all_reports[[sample]]
    df$Sample <- sample
    return(df)
  })
)
#######For NT
# remove nt_rpm is na WHICH is only detected by nr_rpm 
# For nt_rpm detected microbies
nt_combined_df <- combined_df[which(!is.na(combined_df$nt_rpm)),]  



# Add pseudo-count to avoid zero variance
 pseudo_count <- 0.1
# nt_combined_df <- nt_combined_df %>%
#   mutate(
#     nt_rpm = ifelse(is.na(nt_rpm), pseudo_count, nt_rpm + pseudo_count),
#     nr_rpm = ifelse(is.na(nr_rpm), pseudo_count, nr_rpm + pseudo_count)
#   )

# Recalculate control statistics
nt_control_stats <- nt_combined_df %>%
  filter(Sample %in% c("HRR1209919", "HRR1209929")) %>%
  group_by(name) %>%
  summarize(
    nt_mean = mean(nt_rpm, na.rm = TRUE),
    nt_sd = sd(nt_rpm, na.rm = TRUE)
  ) %>%
  # Handle zero SD by setting a minimum value
  #mutate(nt_sd = ifelse(is.na(nt_sd) | nt_sd == 0, pseudo_count, nt_sd))
  mutate(nt_sd = ifelse(is.na(nt_sd), pseudo_count, nt_sd))

# nr_control_stats <- combined_df %>%
#   filter(Sample %in% c("HRR1209919", "HRR1209929")) %>%
#   group_by(name) %>%
#   summarize(
#     nr_mean = mean(nr_rpm, na.rm = TRUE),
#     nr_sd = sd(nr_rpm, na.rm = TRUE)
#   ) %>%
#   mutate(nr_sd = ifelse(is.na(nr_sd) | nr_sd == 0, pseudo_count, nr_sd))

# Recalculate Z-scores
nt_z_score_df <- nt_combined_df %>%
  filter(!Sample %in% c("HRR1209919", "HRR1209929")) %>%
  left_join(nt_control_stats, by = "name") %>%
  #left_join(nr_control_stats, by = "name") %>%
  mutate(
    nt_z_score = (nt_rpm - nt_mean) / nt_sd,
#    nr_z_score = (nr_rpm - nr_mean) / nr_sd,
    # Handle taxa absent in all controls
    nt_z_score = ifelse(is.na(nt_z_score) & !is.na(nt_rpm), Inf, nt_z_score),
#    nr_z_score = ifelse(is.na(nr_z_score) & !is.na(nr_rpm), Inf, nr_z_score),
    # Flag unique taxa
    #is_unique = ifelse(is.na(nt_mean) & is.na(nr_mean), TRUE, FALSE)
  )

idseq_nt_z_score_dna <- nt_z_score_df[nt_z_score_df$nt_z_score >1 ,]
save(idseq_nt_z_score_dna , file = "idseq_nt_z_score_dna.rda")
# # Add z-score and unique flag columns back to the original all_reports list
# for (sample in unique(z_score_df$Sample)) {
#   sample_z <- z_score_df %>% filter(Sample == sample) %>%
#     select(name, nt_z_score, nr_z_score, is_unique)
#   
#   all_reports[[sample]] <- all_reports[[sample]] %>%
#     left_join(sample_z, by = "name")
# }

####### For NR ########
# remove nt_rpm is na WHICH is only detected by nr_rpm 
# For nt_rpm detected microbies
nr_combined_df <- combined_df[which(!is.na(combined_df$nr_rpm)),]  



# Add pseudo-count to avoid zero variance
pseudo_count <- 0.1
# nt_combined_df <- nt_combined_df %>%
#   mutate(
#     nt_rpm = ifelse(is.na(nt_rpm), pseudo_count, nt_rpm + pseudo_count),
#     nr_rpm = ifelse(is.na(nr_rpm), pseudo_count, nr_rpm + pseudo_count)
#   )

# Recalculate control statistics
nr_control_stats <- nr_combined_df %>%
  filter(Sample %in% c("HRR1209919", "HRR1209929")) %>%
  group_by(name) %>%
  summarize(
    nr_mean = mean(nr_rpm, na.rm = TRUE),
    nr_sd = sd(nr_rpm, na.rm = TRUE)
  ) %>%
  # Handle NA SD by setting a minimum value
  mutate(nr_sd = ifelse(is.na(nr_sd), pseudo_count, nr_sd))

# Recalculate Z-scores
nr_z_score_df <- nr_combined_df %>%
  filter(!Sample %in% c("HRR1209919", "HRR1209929")) %>%
  #left_join(nt_control_stats, by = "name") %>%
  left_join(nr_control_stats, by = "name") %>%
  mutate(
    #nt_z_score = (nt_rpm - nt_mean) / nt_sd,
    nr_z_score = (nr_rpm - nr_mean) / nr_sd,
    # Handle taxa absent in all controls
    #nt_z_score = ifelse(is.na(nt_z_score) & !is.na(nt_rpm), Inf, nt_z_score),
    nr_z_score = ifelse(is.na(nr_z_score) & !is.na(nr_rpm), Inf, nr_z_score),
    # Flag unique taxa
    #is_unique = ifelse(is.na(nt_mean) & is.na(nr_mean), TRUE, FALSE)
  )

idseq_nr_z_score_dna <- nr_z_score_df[nr_z_score_df$nr_z_score >1 ,]
save(idseq_nr_z_score_dna , file = "idseq_nr_z_score_dna.rda")




#------------------------------
# 1.2.3 Calculate z scores FOR IDseq RNA SAMPLES
load("all_reports_rna.rda")
all_reports<- all_reports_rna
combined_df <- bind_rows(
  lapply(names(all_reports), function(sample) {
    df <- all_reports[[sample]]
    df$Sample <- sample
    return(df)
  })
)
#######For NT
# remove nt_rpm is na WHICH is only detected by nr_rpm 
# For nt_rpm detected microbies
nt_combined_df <- combined_df[which(!is.na(combined_df$nt_rpm)),]  



# Add pseudo-count to avoid zero variance
pseudo_count <- 0.1
# Recalculate control statistics
nt_control_stats <- nt_combined_df %>%
  filter(Sample %in% c("HRR1209920", "HRR1209930")) %>%
  group_by(name) %>%
  summarize(
    nt_mean = mean(nt_rpm, na.rm = TRUE),
    nt_sd = sd(nt_rpm, na.rm = TRUE)
  ) %>%
  mutate(nt_sd = ifelse(is.na(nt_sd), pseudo_count, nt_sd))

# Recalculate Z-scores
nt_z_score_df <- nt_combined_df %>%
  filter(!Sample %in% c("HRR1209920", "HRR1209930")) %>%
  left_join(nt_control_stats, by = "name") %>%
  #left_join(nr_control_stats, by = "name") %>%
  mutate(
    nt_z_score = (nt_rpm - nt_mean) / nt_sd,
    # Handle taxa absent in all controls
    nt_z_score = ifelse(is.na(nt_z_score) & !is.na(nt_rpm), Inf, nt_z_score),

  )

idseq_nt_z_score_rna <- nt_z_score_df[nt_z_score_df$nt_z_score >1 ,]
save(idseq_nt_z_score_rna , file = "idseq_nt_z_score_rna.rda")


####### For NR ########
# remove nt_rpm is na WHICH is only detected by nr_rpm 
# For nt_rpm detected microbies
nr_combined_df <- combined_df[which(!is.na(combined_df$nr_rpm)),]  



# Add pseudo-count to avoid zero variance
pseudo_count <- 0.1
# Recalculate control statistics
nr_control_stats <- nr_combined_df %>%
  filter(Sample %in% c("HRR1209920", "HRR1209930")) %>%
  group_by(name) %>%
  summarize(
    nr_mean = mean(nr_rpm, na.rm = TRUE),
    nr_sd = sd(nr_rpm, na.rm = TRUE)
  ) %>%
  # Handle NA SD by setting a minimum value
  mutate(nr_sd = ifelse(is.na(nr_sd), pseudo_count, nr_sd))

# Recalculate Z-scores
nr_z_score_df <- nr_combined_df %>%
  filter(!Sample %in% c("HRR1209920", "HRR1209930")) %>%
  #left_join(nt_control_stats, by = "name") %>%
  left_join(nr_control_stats, by = "name") %>%
  mutate(
    #nt_z_score = (nt_rpm - nt_mean) / nt_sd,
    nr_z_score = (nr_rpm - nr_mean) / nr_sd,
    # Handle taxa absent in all controls
    #nt_z_score = ifelse(is.na(nt_z_score) & !is.na(nt_rpm), Inf, nt_z_score),
    nr_z_score = ifelse(is.na(nr_z_score) & !is.na(nr_rpm), Inf, nr_z_score),
    # Flag unique taxa
    #is_unique = ifelse(is.na(nt_mean) & is.na(nr_mean), TRUE, FALSE)
  )

idseq_nr_z_score_rna <- nr_z_score_df[nr_z_score_df$nr_z_score >1 ,]
save(idseq_nr_z_score_rna , file = "idseq_nr_z_score_rna.rda")

#-----------------------------
# 1.2.4 For RpNGS

#  load csv file to combine all abundance tables into one data frame

# Directory where the CSV files are stored
input_dir <- "./pngs/04rawresults"
# Create an empty list to store the data frames
all_reports_dna <- list()
all_reports_rna <- list()
# Loop over each row of the tibble
for (i in 1:nrow(tb2)) {
  # Construct a pattern to search for the file (e.g., HRR1209919_*_taxon_report.csv)
  pattern <- paste0( tb2$Accession[i], ".csv")
  # Get the list of files matching the pattern
  matching_files <- list.files(input_dir, pattern = pattern, full.names = TRUE)
  
  # If a matching file is found, read it
  if (length(matching_files) > 0 ) {
    # Read the first matched file (assuming one file per Accession)
    report <- read_csv(matching_files[1])
    if(tb2$type[i] %in% "DNA"){
      all_reports_dna[[tb2$Accession[i]]] <- report
    } else {
      all_reports_rna[[tb2$Accession[i]]] <- report
    }
    
    
  } else {
    message(paste("No matching file for:", tb2$Accession[i]))
  }
}


# 
save(all_reports_dna ,file = "all_reports_pngsdna.rda")
save(all_reports_rna ,file = "all_reports_pngsrna.rda")


#------------------------------
# 1.2.5  Calculate z scores FOR RpNGS DNA SAMPLES 
load("all_reports_pngsdna.rda")
all_reports<- all_reports_dna
combined_df <- bind_rows(
  lapply(names(all_reports), function(sample) {
    df <- all_reports[[sample]]
    df$Sample <- sample
    return(df)
  })
)

# Add pseudo-count to avoid zero variance
pseudo_count <- 0.1
# Recalculate control statistics
control_stats <- combined_df %>%
  filter(Sample %in% c("HRR1209919", "HRR1209929")) %>%
  group_by(Taxa_names) %>%
  summarize(
    nt_mean = mean(Reads, na.rm = TRUE),
    nt_sd = sd(Reads, na.rm = TRUE)
  ) %>%
  # Handle zero SD by setting a minimum value
  mutate(nt_sd = ifelse(is.na(nt_sd), pseudo_count, nt_sd))

# Recalculate Z-scores
z_score_df <- combined_df %>%
  filter(!Sample %in% c("HRR1209919", "HRR1209929")) %>%
  left_join(control_stats, by = "Taxa_names") %>%
  mutate(
    z_score = (Reads - nt_mean) / nt_sd,
    # Handle taxa absent in all controls
    z_score = ifelse(is.na(z_score), Inf, z_score),
  )

pngs_z_score_dna <- z_score_df[z_score_df$z_score >1 ,]
save(pngs_z_score_dna , file = "pngs_z_score_dna.rda")






#------------------------------
# 1.2.6 Calculate z scores FOR RpNGS RNA SAMPLES
load("all_reports_pngsrna.rda")
all_reports<- all_reports_rna
combined_df <- bind_rows(
  lapply(names(all_reports), function(sample) {
    df <- all_reports[[sample]]
    df$Sample <- sample
    return(df)
  })
)

# Add pseudo-count to avoid zero variance
pseudo_count <- 0.1
# Recalculate control statistics
control_stats <- combined_df %>%
  filter(Sample %in% c("HRR1209920", "HRR1209930")) %>%
  group_by(Taxa_names) %>%
  summarize(
    nt_mean = mean(Reads, na.rm = TRUE),
    nt_sd = sd(Reads, na.rm = TRUE)
  ) %>%
  # Handle zero SD by setting a minimum value
  mutate(nt_sd = ifelse(is.na(nt_sd), pseudo_count, nt_sd))

# Recalculate Z-scores
z_score_df <- combined_df %>%
  filter(!Sample %in% c("HRR1209920", "HRR1209930")) %>%
  left_join(control_stats, by = "Taxa_names") %>%
  mutate(
    z_score = (Reads - nt_mean) / nt_sd,
    # Handle taxa absent in all controls
    z_score = ifelse(is.na(z_score), Inf, z_score),
  )

pngs_z_score_rna <- z_score_df[z_score_df$z_score >1 ,]
save(pngs_z_score_rna , file = "pngs_z_score_rna.rda")

#==========================================================================
# 1.2.7 MAKE P, R, F1 figures 

# Common data 
# Read the Excel file
tb1 <- read_excel("table_S1.xlsx", sheet = 1)
tb1 <- tb1[-24:-25,]
#tb1_dna <- tb1[-which(tb1$Type %in% "RNA virus"),]
#---
tb2 <- read_excel("table_S2.xlsx", sheet = 1)
# Extract "DNA" and "RNA"
tb2$type <- str_extract(tb2$Experiment_title, "(DNA|RNA)")

tb1idseq <- read_excel("table_S1idseq.xlsx", sheet = 1)
tb1idseq <- tb1idseq[-24:-25,]



# Create an empty data table  to store the data frames
zscore_all_statdt <- data.frame()
# For Rpngs
load("pngs_z_score_rna.rda")
load("pngs_z_score_dna.rda")
# IDseq_nt
load("idseq_nt_z_score_dna.rda")
load("idseq_nt_z_score_rna.rda")

# IDseq_nr
load("idseq_nr_z_score_dna.rda")
load("idseq_nr_z_score_rna.rda")
# Loop over each row of the tibble
for (i in 1:nrow(tb2)) {
    if (tb2$type[i] %in% "DNA") {
      
      # rpngs
      tb1_groundtruth_dnatax <- tb1[-which(tb1$Type %in% "RNA virus"),c("Microbes",tb2$Sample_Name[i])]
      colnames(tb1_groundtruth_dnatax) <- c("Microbes","value")
      tb1_groundtruth_dnatax <- tb1_groundtruth_dnatax[which(tb1_groundtruth_dnatax$value != "\\"),]

      pngs_taxa <- pngs_z_score_dna[which(pngs_z_score_dna$Sample %in% tb2$Accession[i]),]
      TP <- length(unique(intersect(pngs_taxa$Taxa_names, tb1_groundtruth_dnatax$Microbes)))
      FP <- length(unique(pngs_taxa$Taxa_names)) - TP
      
      FN <- length(tb1_groundtruth_dnatax$Microbes) - TP
      P <- TP/(TP+FP)
      R <- TP/(TP+FN)
      F1=2/((1/R) + (1/P))
      rpm_type <- "kmer"
      platform <- "RpNGS"
      statdt_dna <- data.frame(tb2$Sample_Name[i],TP,FP,FN,P,R,F1,tb2$type[i],tb2$Accession[i],rpm_type,platform)
      # Store the data frame in the list
      zscore_all_statdt<- rbind(zscore_all_statdt,statdt_dna )
      
      
      # For idseq_nt
      tb1_groundtruth_dnatax <- tb1idseq[-which(tb1idseq$Type %in% "RNA virus"),c("Microbes",tb2$Sample_Name[i])]
      colnames(tb1_groundtruth_dnatax) <- c("Microbes","value")
      tb1_groundtruth_dnatax <- tb1_groundtruth_dnatax[which(tb1_groundtruth_dnatax$value != "\\"),]
      
      idseq_taxa <- idseq_nt_z_score_dna[which(idseq_nt_z_score_dna$Sample %in% tb2$Accession[i]),]
      TP <- length(unique(intersect(idseq_taxa$name, tb1_groundtruth_dnatax$Microbes)))
      FP <- length(unique(idseq_taxa$name)) - TP
      
      FN <- length(tb1_groundtruth_dnatax$Microbes) - TP
      P <- TP/(TP+FP)
      R <- TP/(TP+FN)
      F1=2/((1/R) + (1/P))
      rpm_type <- "NT"
      platform <- "IDseq"
      statdt_dna <- data.frame(tb2$Sample_Name[i],TP,FP,FN,P,R,F1,tb2$type[i],tb2$Accession[i],rpm_type,platform)
      # Store the data frame in the list
      zscore_all_statdt<- rbind(zscore_all_statdt,statdt_dna )
      
      # For idseq_nr
      idseq_taxa <- idseq_nr_z_score_dna[which(idseq_nr_z_score_dna$Sample %in% tb2$Accession[i]),]
      TP <- length(unique(intersect(idseq_taxa$name, tb1_groundtruth_dnatax$Microbes)))
      FP <- length(unique(idseq_taxa$name)) - TP
      
      FN <- length(tb1_groundtruth_dnatax$Microbes) - TP
      P <- TP/(TP+FP)
      R <- TP/(TP+FN)
      F1=2/((1/R) + (1/P))
      rpm_type <- "NR"
      platform <- "IDseq"
      statdt_dna <- data.frame(tb2$Sample_Name[i],TP,FP,FN,P,R,F1,tb2$type[i],tb2$Accession[i],rpm_type,platform)
      # Store the data frame in the list
      zscore_all_statdt<- rbind(zscore_all_statdt,statdt_dna )

      
    } else {
      tb1_groundtruth_rnatax <- tb1[which(tb1$Type %in% "RNA virus"),c("Microbes",tb2$Sample_Name[i])]
      colnames(tb1_groundtruth_rnatax) <- c("Microbes","value")
      tb1_groundtruth_rnatax <- tb1_groundtruth_rnatax[which(tb1_groundtruth_rnatax$value != "\\"),]
      
      pngs_taxa <- pngs_z_score_rna[which(pngs_z_score_rna$Sample %in% tb2$Accession[i]),]
      TP <- length(unique(intersect(pngs_taxa$Taxa_names, tb1_groundtruth_rnatax$Microbes)))
      FP <- length(unique(pngs_taxa$Taxa_names)) - TP
      
      FN <- length(tb1_groundtruth_rnatax$Microbes) - TP
      P <- TP/(TP+FP)
      R <- TP/(TP+FN)
      F1=2/((1/R) + (1/P))
      rpm_type <- "kmer"
      platform <- "RpNGS"
      statdt_dna <- data.frame(tb2$Sample_Name[i],TP,FP,FN,P,R,F1,tb2$type[i],tb2$Accession[i],rpm_type,platform)
      # Store the data frame in the list
      zscore_all_statdt<- rbind(zscore_all_statdt,statdt_dna )
      
      # For idseq_nt
      tb1_groundtruth_rnatax <- tb1idseq[which(tb1idseq$Type %in% "RNA virus"),c("Microbes",tb2$Sample_Name[i])]
      colnames(tb1_groundtruth_rnatax) <- c("Microbes","value")
      tb1_groundtruth_rnatax <- tb1_groundtruth_rnatax[which(tb1_groundtruth_rnatax$value != "\\"),]
      
      idseq_taxa <- idseq_nt_z_score_rna[which(idseq_nt_z_score_rna$Sample %in% tb2$Accession[i]),]
      TP <- length(unique(intersect(idseq_taxa$name, tb1_groundtruth_rnatax$Microbes)))
      FP <- length(unique(idseq_taxa$name)) - TP
      
      FN <- length(tb1_groundtruth_rnatax$Microbes) - TP
      P <- TP/(TP+FP)
      R <- TP/(TP+FN)
      F1=2/((1/R) + (1/P))
      rpm_type <- "NT"
      platform <- "IDseq"
      statdt_dna <- data.frame(tb2$Sample_Name[i],TP,FP,FN,P,R,F1,tb2$type[i],tb2$Accession[i],rpm_type,platform)
      # Store the data frame in the list
      zscore_all_statdt<- rbind(zscore_all_statdt,statdt_dna )
      
      # For idseq_nr
      idseq_taxa <- idseq_nr_z_score_rna[which(idseq_nr_z_score_rna$Sample %in% tb2$Accession[i]),]
      TP <- length(unique(intersect(idseq_taxa$name, tb1_groundtruth_rnatax$Microbes)))
      FP <- length(unique(idseq_taxa$name)) - TP
      
      FN <- length(tb1_groundtruth_rnatax$Microbes) - TP
      P <- TP/(TP+FP)
      R <- TP/(TP+FN)
      F1=2/((1/R) + (1/P))
      rpm_type <- "NR"
      platform <- "IDseq"
      statdt_dna <- data.frame(tb2$Sample_Name[i],TP,FP,FN,P,R,F1,tb2$type[i],tb2$Accession[i],rpm_type,platform)
      # Store the data frame in the list
      zscore_all_statdt<- rbind(zscore_all_statdt,statdt_dna )
      
      
      
    }
    
    
}
save(zscore_all_statdt, file="zscore_all_statdt.rda")
#=====================
load("zscore_all_statdt.rda")
df <- zscore_all_statdt
# Load necessary libraries
library(ggplot2)
library(dplyr)

### F-scores Plot (using F1 column) ###
df_F1 <- df %>%
  select(tb2.Sample_Name.i., tb2.type.i., rpm_type, platform, F1) %>%
  rename(score = F1)

plot_F1 <- ggplot(df_F1, aes(x = tb2.Sample_Name.i., y = score)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  facet_grid(tb2.type.i. ~ rpm_type + platform) +
  labs(
    title = "F-scores",
    x = "Sample Name",
    y = "F-scores"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

### Recall Plot (using R column) ###
df_R <- df %>%
  select(tb2.Sample_Name.i., tb2.type.i., rpm_type, platform, R) %>%
  rename(score = R)

plot_R <- ggplot(df_R, aes(x = tb2.Sample_Name.i., y = score)) +
  geom_bar(stat = "identity", fill = "darkgreen") +
  facet_grid(tb2.type.i. ~ rpm_type + platform) +
  labs(
    title = "Recall",
    x = "Sample Name",
    y = "Recall"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

### Precision Plot (using P column) ###
df_P <- df %>%
  select(tb2.Sample_Name.i., tb2.type.i., rpm_type, platform, P) %>%
  rename(score = P)

plot_P <- ggplot(df_P, aes(x = tb2.Sample_Name.i., y = score)) +
  geom_bar(stat = "identity", fill = "firebrick") +
  facet_grid(tb2.type.i. ~ rpm_type + platform) +
  labs(
    title = "Precision",
    x = "Sample Name",
    y = "Precision"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# To view the plots separately:
pdf(file = "zscore_fscore.pdf",width = 12, height = 8,)
plot(plot_F1)
dev.off()

pdf(file = "zscore_r.pdf",width = 12, height = 8,)
plot(plot_R)
dev.off()

pdf(file = "zscore_p.pdf",width = 12, height = 8,)
plot(plot_P)
dev.off()

#=========================================================
#DONE





# --------------------------------------------------------------------------------
# HETAMAP FOR RpNGS clinical samples 

library(ggplot2)
library(reshape2)
library(pheatmap)

# Define file paths
file_paths <- list(
  "./clinicalsample_se/04rawresults/CHRF0000_1.csv",
  "./clinicalsample_se/04rawresults/CHRF0002_1.csv",
  "./clinicalsample_se/04rawresults/CHRF0094_1.csv"
)

# Define subset of taxa
subset_taxa <- c(
  "Cutibacterium acnes", "Sphingomonas sp. AAP5", "Streptococcus pneumoniae",
  "Corynebacterium kefirresidentii", "Malassezia restricta", "Staphylococcus warneri",
  "Streptococcus mitis", "Lactobacillus iners", "Micrococcus luteus", "Staphylococcus epidermidis",
  "Salmonella enterica", "Malassezia globosa", "Streptococcus oralis", "Chikungunya virus",
  "Escherichia coli", "Solanum pennellii", "uncultured bacterium", "Porcine picobirnavirus",
  "Streptomyces laurentii"
)

# Read and merge data
data_list <- lapply(seq_along(file_paths), function(i) {
  df <- read.csv(file_paths[[i]])
  df$Sample <- paste0("Sample_", i)  # Assign sample names
  return(df)
})

# Combine all data into one DataFrame
df_combined <- do.call(rbind, data_list)

# Filter for selected taxa
df_filtered <- subset(df_combined, Taxa_names %in% subset_taxa)

# Reshape data for heatmap (wide format)
df_pivot <- dcast(df_filtered, Taxa_names ~ Sample, value.var = "Reads", fun.aggregate = sum)

# Set row names
rownames(df_pivot) <- df_pivot$Taxa_names
df_pivot$Taxa_names <- NULL

# Convert NA to 0
df_pivot[is.na(df_pivot)] <- 0
write.csv(df_pivot, "rpngs_clinicalheatmap.csv")
# Generate heatmap

pdf(file = "rpngs_clinicalheatmap.pdf",width = 12, height = 8,)
pheatmap(
  df_pivot,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("white", "red"))(50),
  main = "Heatmap of reads across three samples"
)
dev.off()
