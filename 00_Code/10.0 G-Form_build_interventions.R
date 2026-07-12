# 10.0 G-Formula — construcción de objetos de intervención (ejecutar una vez) ----
#
# Genera todos los RDS de intervenciones globales en:
#   02_Output/G-Form/Interventions/
#
# Uso (desde la raíz del proyecto):
#   Rscript "00_Code/10.0 G-Form_build_interventions.R"
#   GFORM_EXEC_MODE=server Rscript "00_Code/10.0 G-Form_build_interventions.R"

rm(list = ls())

## Settings ----
source("00_Code/0.2 Packages_gform.R")
source("00_Code/0.1 Settings.R")
source("00_Code/0.3 Functions.R")
source("00_Code/10.1 G-Form_functions.R")

data_inp <- "01_Data/Output/"
data_out_g <- "02_Output/G-Form/"
dir_interventions <- file.path(data_out_g, "Interventions")

## Parámetros ----
overwrite_interventions <- FALSE

execution_mode <- tolower(Sys.getenv("GFORM_EXEC_MODE", "auto"))
if (identical(execution_mode, "auto")) {
  execution_mode <- if (gform_is_linux_server()) "server" else "local"
}

parallel_config <- gform_parallel_config()
options(gform.parallel = parallel_config)
run_parallel <- execution_mode == "server"

on.exit(gform_finalize_run(), add = TRUE)

message("=== G-Formula Etapa 1: construcción de intervenciones ===")
message("Modo: ", execution_mode)
message("Destino: ", dir_interventions)
if (run_parallel) {
  gform_setup_parallel(task = "build", config = parallel_config)
  message("Construcción paralela de ", length(GFORM_INTERVENTION_REGISTRY), " intervenciones.")
} else {
  message("Construcción secuencial.")
}

timing_log <- gform_timing_log_init()

block <- gform_time_block("Cargar datos de cohorte", {
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
  data_long
})
timing_log <- gform_timing_log_add(timing_log, block$timing)
data_long <- block$result

built <- build_all_interventions(
  data_long = data_long,
  registry = GFORM_INTERVENTION_REGISTRY,
  output_dir = dir_interventions,
  overwrite = overwrite_interventions,
  parallel = run_parallel,
  n_workers = if (run_parallel) parallel_config$n_workers_build else NULL
)
timing_log <- gform_timing_log_merge(timing_log, built$timing)
built_paths <- built$paths

timing_log$finished_at <- Sys.time()
timing_log$wall_sec <- as.numeric(
  difftime(timing_log$finished_at, timing_log$started_at, units = "secs")
)
timing_log$total_sec <- gform_timing_total_sec(timing_log)
gform_print_timing_summary(timing_log, "Etapa 1 — construcción de intervenciones")

message("\nIntervenciones generadas (", length(built_paths), "):")
for (p in built_paths) message("  - ", p)

beepr::beep(8)
