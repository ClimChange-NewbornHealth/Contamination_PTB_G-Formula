# 11.1 Exposure–perinatal outcome models — estimación Cox ----
#
# Modelos Cox por contaminante (PM2.5, NO2, O3; kriging):
#   - período completo (tot): modelo simple
#   - trimestres gestacionales (t1+t2+t3): modelo conjunto
# Escalas: raw e IQR (variables *_iqr en births_2010_2020_exposure.RData)
#
# Uso (desde la raíz del proyecto):
#   Rscript "00_Code/11.1 Exposure_PO_Models.R"
#
# Salida:
#   02_Output/Exposure_PO/Models/Exposure_models_PO_cox.xlsx
#   02_Output/Exposure_PO/Models/Exposure_models_PO_cox.RData
#   02_Output/Exposure_PO/Models/List_models_exposure_PO.xlsx

rm(list = ls())
gc(verbose = FALSE)

source("00_Code/0.1 Settings.R")
source("00_Code/0.2 Packages.R")
source("00_Code/11.0 Exposure_PO_functions.R")

data_inp <- "01_Data/Output/"
data_out <- "02_Output/Exposure_PO/"
dir_models <- file.path(data_out, "Models")
dir.create(dir_models, recursive = TRUE, showWarnings = FALSE)

## 1. Cargar y preparar datos ----

births_exposure <- rio::import(file.path(data_inp, "births_2010_2020_exposure.RData"))

message("Verificando variables de exposición en births_2010_2020_exposure.RData...")
required_full_tri <- unlist(lapply(
  PO_EXPOSURE_POLLUTANTS,
  function(p) {
    c(
      po_exposure_var(p, "tot", "raw"),
      po_exposure_var(p, "tot", "iqr"),
      po_exposure_var(p, "t1", "raw"),
      po_exposure_var(p, "t2", "raw"),
      po_exposure_var(p, "t3", "raw"),
      po_exposure_var(p, "t1", "iqr"),
      po_exposure_var(p, "t2", "iqr"),
      po_exposure_var(p, "t3", "iqr")
    )
  }
))
missing_full_tri <- setdiff(required_full_tri, names(births_exposure))
if (length(missing_full_tri)) {
  stop(
    "Faltan variables de exposición full/trimestre: ",
    paste(missing_full_tri, collapse = ", ")
  )
}
message("  OK: full + t1/t2/t3 (raw e IQR) presentes para pm25, no2, o3.")

dependent_vars <- c(
  "birth_preterm",
  "birth_very_preterm",
  "birth_moderately_preterm",
  "birth_late_preterm",
  "lbw",
  "tlbw",
  "sga"
)

control_vars <- c(
  "sex", "age_group_mom", "educ_group_mom", "job_group_mom",
  "age_group_dad", "educ_group_dad", "job_group_dad",
  "month_week1", "year_week1", "covid", "vulnerability",
  "tad_full", "ndvi_full"
)

data_model <- births_exposure |>
  dplyr::filter(!is.na(.data$weeks), .data$weeks >= 28L) |>
  dplyr::mutate(
    tstart = 27L,
    month_week1 = factor(.data$month_week1),
    year_week1 = factor(.data$year_week1),
    covid = factor(.data$covid)
  )

rm(births_exposure)
gc(verbose = FALSE)

missing_deps <- setdiff(dependent_vars, names(data_model))
if (length(missing_deps)) {
  stop("Faltan outcomes: ", paste(missing_deps, collapse = ", "))
}
missing_controls <- setdiff(control_vars, names(data_model))
if (length(missing_controls)) {
  stop("Faltan covariables: ", paste(missing_controls, collapse = ", "))
}

available_predictors <- names(data_model)[
  grepl("_krg_(full|t1|t2|t3)", names(data_model)) |
    grepl("_krg_(full|t1|t2|t3)_iqr", names(data_model))
]

message(
  "Datos modelados: n = ", nrow(data_model),
  " | outcomes = ", length(dependent_vars),
  " | predictores disponibles = ", length(available_predictors)
)

## 2. Grilla de modelos ----

combinations <- build_po_model_grid(
  dependent_vars = dependent_vars,
  available_predictors = available_predictors
)

writexl::write_xlsx(
  combinations,
  path = file.path(dir_models, "List_models_exposure_PO.xlsx")
)
message("Grilla de modelos: ", nrow(combinations), " especificaciones Cox.")

## 3. Estimación Cox ----

run_po_model_grid <- function(combinations, fit_fun, label) {
  message("Estimando modelos ", label, " (n = ", nrow(combinations), ")...")
  if (requireNamespace("furrr", quietly = TRUE) &&
      requireNamespace("future", quietly = TRUE)) {
    n_workers <- as.integer(Sys.getenv("PO_COX_WORKERS", unset = "4"))
    if (!is.finite(n_workers) || n_workers < 1L) n_workers <- 4L
    message("Paralelización: ", n_workers, " workers (multisession).")
    future::plan(future::multisession, workers = n_workers)
    options(future.globals.maxSize = 2 * 1024^3)
    on.exit(future::plan(future::sequential), add = TRUE)
    furrr::future_map(
      seq_len(nrow(combinations)),
      function(i) {
        row <- combinations[i, , drop = FALSE]
        fit_fun(
          dependent = row$dependent[[1L]],
          predictor = row$predictor[[1L]],
          tiempo = row$tiempo[[1L]],
          contaminante = row$contaminante[[1L]],
          tipo = row$tipo[[1L]],
          model_type = row$model_type[[1L]],
          data = data_model,
          control_vars = control_vars,
          adjustment = row$adjustment[[1L]],
          exposure_scale = row$exposure_scale[[1L]]
        )
      },
      .options = furrr::furrr_options(seed = TRUE)
    )
  } else {
    lapply(seq_len(nrow(combinations)), function(i) {
      row <- combinations[i, , drop = FALSE]
      fit_fun(
        dependent = row$dependent[[1L]],
        predictor = row$predictor[[1L]],
        tiempo = row$tiempo[[1L]],
        contaminante = row$contaminante[[1L]],
        tipo = row$tipo[[1L]],
        model_type = row$model_type[[1L]],
        data = data_model,
        control_vars = control_vars,
        adjustment = row$adjustment[[1L]],
        exposure_scale = row$exposure_scale[[1L]]
      )
    })
  }
}

results_list_cox <- run_po_model_grid(
  combinations,
  fit_fun = fit_po_cox_model,
  label = "Cox"
)

results_cox <- dplyr::bind_rows(results_list_cox)

## 4. Guardar ----

save(
  results_list_cox,
  file = file.path(dir_models, "Exposure_models_PO_cox.RData")
)

writexl::write_xlsx(
  list(cox_models = results_cox),
  path = file.path(dir_models, "Exposure_models_PO_cox.xlsx")
)

message("Resultados guardados en: ", dir_models)
message("  Filas Cox: ", nrow(results_cox))
