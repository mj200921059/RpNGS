# data_prep.R
# 该脚本用于加载和预处理数据

# 使用pacman管理包
if (!require(pacman)) install.packages("pacman")
pacman::p_load("sf","shiny","shinyBS","shinycssloaders","shinyFiles","plotly","leaflet","data.table",
               "DT","lubridate","readr","tidyverse","shinyjs","shinydashboard","shinyWidgets","Gviz",
               "Rsamtools","GenomicAlignments","officer","flextable","rvest","ShortRead","scales",
               "ps","future","promises","dplyr","markdown","shinythemes")

# 读取 `RpNGS_sequencinginfo.csv` 和 `RpNGS_clinicalinfo.csv`
seqinfo_file <- "./data/RpNGS_sequencinginfo.csv"
clininfo_file <- "./data/RpNGS_clinicalinfo.csv"
pathogen_human <- "./data/human_pathogens.csv"
rawfastq_file <- "./mnt/"
results_folders <- "./results"
# load('~/dataanalysis/projects/lifegen/pre_datasets/microbes_taxonomy.rda')
# microbes_taxonomy <- microbes_taxonomy[,-1]

# 读取 sequencing info 数据
# Load previous data if it exists
if (file.exists(seqinfo_file)) {
  datafr1 <- read.csv(seqinfo_file, stringsAsFactors = FALSE)
} else {
  datafr1 <- data.frame(
    Chip_id = character(),
    Sample_id = character(),
    Extracted_NAC = numeric(),
    Adaptor = numeric(),
    Library_NAC = numeric(),
    Rawfastq_id = character(),
    stringsAsFactors = FALSE
  )
}
#datafr1 <- as.data.frame(fread(seqinfo_file), stringsAsFactors = FALSE)

# 读取 clinical info 数据
datafr2 <- as.data.frame(fread(clininfo_file), stringsAsFactors = FALSE)

# 处理时间格式
datafr2$Test_day <- as.Date(datafr2$Test_day, format = "%d-%m-%Y")

# 只选取包含 Test_day 的样本
sub_datafr2 <- datafr2 %>%
  filter(!is.na(Test_day)) %>%
  mutate(
    month = format(Test_day, "%m"),  # 提取月份
    year = paste0("20", format(Test_day, "%y")),  # 提取年份
    season = case_when(
      month %in% c("12", "01", "02") ~ "Winter",
      month %in% c("03", "04", "05") ~ "Spring",
      month %in% c("06", "07", "08") ~ "Summer",
      month %in% c("09", "10", "11") ~ "Fall",
      TRUE ~ NA_character_
    )
  )

# 仅保留当前年份的数据
current_year <- format(Sys.Date(), "%Y")
# summary_sub_datafr2 <- sub_datafr2 %>% filter(year == current_year)
summary_sub_datafr2 <- sub_datafr2 %>% filter(year == "2024") # For demo, I chose 2024
# 统计各个 `Location` 的 `Month` 和 `Season` 统计数据
# summary_sub_datafr2 <- summary_sub_datafr2 %>%
#   group_by(Location, season) %>%
#   summarise(count = n(), .groups = "drop") %>%
#   pivot_wider(names_from = season, values_from = count, values_fill = 0)

# Summarize counts for each month and season
summary_sub_datafr2 <- summary_sub_datafr2 %>%
  mutate(month_col = paste0("M", month)) %>% # Add a month column for pivoting
  pivot_longer(cols = c(month_col, season), names_to = "category", values_to = "value") %>%
  group_by(Location, value) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(names_from = value, values_from = count, values_fill = 0)



# 确保数据完整，补充缺失的月份
add_cols <- function(df, cols) {
  missing_cols <- setdiff(cols, names(df))
  df[missing_cols] <- 0
  return(df)
}
summary_sub_datafr2 <- add_cols(summary_sub_datafr2, c('Winter', 'Spring', 'Summer','Fall','M01','M02','M03','M04','M05','M06','M07','M08','M09','M10','M11','M12'))

# Arrange the columns in the desired order
summary_sub_datafr2 <- summary_sub_datafr2 %>%
  select(Location, starts_with("M"), Winter, Spring, Summer, Fall)

# Convert the data frame to long format, keeping 'Location' intact
summary_sub_datafr2 <- summary_sub_datafr2 %>%
  pivot_longer(cols = -Location, 
               names_to = "Month_Season", 
               values_to = "Value")


# Add fourth column for labeling as "Month" or "Season"
summary_sub_datafr2 <- summary_sub_datafr2 %>%
  mutate(Mark = case_when(
    grepl("^M", Month_Season) ~ "Month",  # Labels months based on 'M' prefix
    Month_Season %in% c("Summer", "Spring", "Fall", "Winter") ~ "Season",  # Labels seasons
    TRUE ~ "Other"   # In case of any unexpected column names
  ))

# 导出数据（如果需要）
# write.csv(summary_sub_datafr2, "summary_sub_datafr2.csv", row.names = FALSE)


#----------------------------------------------------------
# Summary pathogen based on each infection 
#----------------------------------------------------------
# Load required libraries
library(dplyr)
library(tidyr)
library(stringr)

# Read the data
#clinical_data <- read.csv("RpNGS_clinicalinfo.csv", stringsAsFactors = FALSE)

# Clean and prepare the data
pathogen_data <- datafr2 %>%
  # Remove empty rows
  filter(!is.na(Infections) & Infections != "") %>%
  # Select relevant columns
  select(Sample_id, Infections, Pathogens) %>%
  # Remove rows with no pathogen data
  filter(!is.na(Pathogens) & Pathogens != "") %>%
  # Separate pathogens into individual rows
  mutate(Pathogens = str_split(Pathogens, ";")) %>%
  unnest(Pathogens) %>%
  # Clean up pathogen names and extract read counts
  mutate(
    Pathogens = str_trim(Pathogens),
    pathogen_name = str_extract(Pathogens, "^[A-Za-z ]+"),
    read_count = as.numeric(str_extract(Pathogens, "[0-9]+"))
  ) %>%
  # Remove any remaining empty rows
  filter(!is.na(pathogen_name) & pathogen_name != "")

# Summarize pathogens by infection type
pathogen_summary <- pathogen_data %>%
  group_by(Infections, pathogen_name) %>%
  summarize(
    total_reads = sum(read_count, na.rm = TRUE),
    sample_count = n(),
    .groups = "drop"
  ) %>%
  arrange(Infections, desc(total_reads))

# View the summary
print(pathogen_summary)

# Optional: Create a more detailed view
detailed_summary <- pathogen_summary %>%
  group_by(Infections) %>%
  mutate(
    percent_of_total = round(total_reads / sum(total_reads) * 100, 1)
  ) %>%
  arrange(Infections, desc(total_reads))

#===============================================================================
# analysis module
#----------------------------------
# Ensure input_data exists before calling the analysis module
sample_data <- data.frame(
  Chip_id = character(),
  Sample_id = character(),
  Extracted_NAC = numeric(),
  Adaptor = numeric(),
  Library_NAC = numeric(),
  Rawfastq_id = character(),
  stringsAsFactors = FALSE
)

# 创建所需的目录结构
create_dir <- function(dir) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  }
}
