# Code 6: Descriptive expsition ----

rm(list=ls())
## Settings ----
source("00_Code/0.1 Settings.R")
source("00_Code/0.2 Packages.R")
source("00_Code/0.3 Functions.R")

# Data path 
data_inp <- "01_Data/Output/"
data_out <- "01_Data/Output/"

## Load data ----
exposure <- rio::import(paste0(data_inp, "Contamination_Climate_Data_2010_2020.RData"))
