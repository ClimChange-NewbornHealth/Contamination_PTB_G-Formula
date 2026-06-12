# 10.0 G-Formula — construcción de objetos de intervención (ejecutar una vez) ----
#
# Genera todos los RDS de intervenciones globales en:
#   02_Output/G-Form/Interventions/
#
# Uso (desde la raíz del proyecto):
#   Rscript "00_Code/10.0 G-Form_build_interventions.R"

rm(list = ls())

## Settings ----
source("00_Code/0.1 Settings.R")
source("00_Code/0.2 Packages.R")
source("00_Code/0.3 Functions.R")
source("00_Code/10.1 G-Form_functions.R")

data_inp <- "01_Data/Output/"
data_out_g <- "02_Output/G-Form/"
dir_interventions <- file.path(data_out_g, "Interventions")

## Parámetros ----
overwrite_interventions <- FALSE
n_workers <- 10L

message("=== G-Formula Etapa 1: construcción de intervenciones ===")
message("Workers: ", n_workers)
message("Destino: ", dir_interventions)

## Carga de datos (cohorte completa) ----
message("Cargando datos...")

births <- rio::import(paste0(data_inp, "births_exposure_period_metrics_full_d30_d4_tri.RData")) |>
  dplyr::select(id, ndvi_full) |>
  dplyr::distinct(id, .keep_all = TRUE)

load(paste0(data_inp, "births_2010_2020_exposure_weeks_lagged.RData"))

data_long <- data_long |>
  dplyr::left_join(births, by = "id") |>
  dplyr::mutate(
    month_week1 = factor(month_week1),
    year_week1 = factor(year_week1),
    covid = factor(covid)
  )

message("Nacimientos: ", dplyr::n_distinct(data_long$id))

## Construcción paralela ----
tictoc::tic("build_all_interventions")
built_paths <- build_all_interventions(
  data_long = data_long,
  registry = GFORM_INTERVENTION_REGISTRY,
  output_dir = dir_interventions,
  overwrite = overwrite_interventions,
  parallel = FALSE
)
tictoc::toc()

message("\nIntervenciones generadas (", length(built_paths), "):")
for (p in built_paths) message("  - ", p)

beepr::beep(8)
