# Check and install required packages

required_packages <- c("sf","shiny","shinyBS","shinycssloaders","shinyFiles","plotly","leaflet","data.table",
                       "DT","lubridate","readr","tidyverse","shinyjs","shinydashboard","shinyWidgets","Gviz",
                       "Rsamtools","GenomicAlignments","officer","flextable","rvest","ShortRead","scales",
                       "ps","future","promises")

new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

# 
# PACKAGES USED IN THIS APP
# JUN MA
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