# 16.1 G-Formula por períodos multicontaminante — construcción de intervenciones ----
#
# Genera 6 RDS (pct20 × 3 contaminantes × 2 alcances) en:
#   02_Output/G-Form-Period-Multi/Interventions/
#
# Trimestre: reduce T1/T2/T3 del focal; co-contaminantes en curso natural.
# Overall:   reduce full del focal; co-contaminantes en curso natural.
#
# Uso (desde la raíz del proyecto, preferible en servidor):
#   Rscript "00_Code/16.1 G-Form_period_multi_build_interventions.R"
#   GFORM_EXEC_MODE=server Rscript "00_Code/16.1 G-Form_period_multi_build_interventions.R"

rm(list = ls())

source("00_Code/0.2 Packages_gform.R")
source("00_Code/0.1 Settings.R")
source("00_Code/0.3 Functions.R")
source("00_Code/16.0 G-Form_period_multi_functions.R")

data_inp <- "01_Data/Output"
data_out_g <- GFORM_PERIOD_MULTI_DATA_OUT
dir_interventions <- file.path(data_out_g, "Interventions")

overwrite_interventions <- gform_env_bool("GFORM_PERIOD_MULTI_OVERWRITE", default = FALSE)

execution_mode <- tolower(Sys.getenv("GFORM_EXEC_MODE", "auto"))
if (identical(execution_mode, "auto")) {
  execution_mode <- if (gform_is_linux_server()) "server" else "local"
}

parallel_config <- gform_parallel_config()
options(gform.parallel = parallel_config)
run_parallel <- execution_mode == "server"

on.exit(gform_finalize_run(), add = TRUE)

message("=== G-Formula Período Multi — Etapa 1: construcción de intervenciones ===")
message("Modo: ", execution_mode)
message("Destino: ", dir_interventions)
message("Escenarios: ", length(GFORM_PERIOD_MULTI_INTERVENTION_REGISTRY))
print_gform_period_multi_intervention_menu()

if (run_parallel) {
  gform_setup_parallel(task = "build", config = parallel_config)
  message("Construcción paralela de ", length(GFORM_PERIOD_MULTI_INTERVENTION_REGISTRY), " RDS.")
} else {
  message("Construcción secuencial.")
}

timing_log <- gform_timing_log_init()

block <- gform_time_block("Cargar cohorte con exposiciones agregadas", {
  births_df <- prepare_period_births_cohort(data_inp)
  message("Nacimientos: ", nrow(births_df))
  births_df
})
timing_log <- gform_timing_log_add(timing_log, block$timing)
births_df <- block$result

built <- build_all_period_multi_interventions(
  births_df = births_df,
  registry = GFORM_PERIOD_MULTI_INTERVENTION_REGISTRY,
  output_dir = dir_interventions,
  overwrite = overwrite_interventions,
  parallel = run_parallel,
  n_workers = if (run_parallel) parallel_config$n_workers_build else NULL
)
timing_log <- gform_timing_log_merge(timing_log, built$timing)
built_paths <- built$paths

rm(births_df)
gc(verbose = FALSE)

timing_log$finished_at <- Sys.time()
timing_log$wall_sec <- as.numeric(
  difftime(timing_log$finished_at, timing_log$started_at, units = "secs")
)
timing_log$total_sec <- gform_timing_total_sec(timing_log)
gform_print_timing_summary(timing_log, "Etapa 1 — intervenciones período multi")

message("\nIntervenciones período multi generadas (", length(built_paths), "):")
for (p in built_paths) message("  - ", p)

beepr::beep(8)
