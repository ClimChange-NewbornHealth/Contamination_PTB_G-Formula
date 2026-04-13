# Code 6: Descriptive exposition — exposure summary table and plots ----

rm(list = ls())

## Settings ----
source("00_Code/0.1 Settings.R")
source("00_Code/0.2 Packages.R")

data_inp <- "01_Data/Output/"
data_out <- "02_Output/Descriptives/"

## Load exposure (contaminación + metadatos comunales: lat, long, sup) ----
exposure <- rio::import(paste0(data_inp, "Contamination_Climate_Data_2010_2020.RData"))
glimpse(exposure)

