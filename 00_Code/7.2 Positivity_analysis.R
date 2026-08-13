# 7.2 Positivity analysis — weekly births and exposure ----
#
# Tablas por semana de gestación: nacimientos, pretérminos, % PTB y estadísticos
# de exposición (Kriging) para PM2.5, NO2 y O3 (cohorte completa y solo PTB).
# Tabla adicional: exposición por semana gestacional (desde la primera a la última).
#
# Output: 02_Output/Descriptives/Table_positivity_by_week.xlsx

rm(list = ls())

## Settings ----
source("00_Code/0.1 Settings.R")
source("00_Code/0.2 Packages.R")

data_inp <- "01_Data/Output/"
data_out <- "02_Output/Descriptives/"
week_start <- 28L

pollutants <- list(
  pm25_krg = list(label = "PM2.5", unit = "ug/m3"),
  no2_krg  = list(label = "NO2",   unit = "ppbv"),
  o3_krg   = list(label = "O3",    unit = "ppbv")
)

dir.create(data_out, recursive = TRUE, showWarnings = FALSE)

## Remove previous outputs from this script ----
old_outputs <- c(
  "Table_positivity_PT_births.xlsx",
  "Table_positivity_by_week.xlsx",
  "Positivity_joint_expo_lag_weeks28_36.png",
  "Positivity_joint_expo_lag_by_pollutant.png",
  "Natural_course_coefs_weeks28_36.png"
)
invisible(lapply(file.path(data_out, old_outputs), function(f) {
  if (file.exists(f)) file.remove(f)
}))

## Load data ----
births_weeks <- rio::import(paste0(data_inp, "births_2010_2020_exposure_weeks.RData"))

births_weeks <- births_weeks |>
  dplyr::mutate(
    birth_preterm = if ("birth_preterm" %in% names(births_weeks)) {
      birth_preterm
    } else {
      as.integer(weeks < 37L)
    }
  )

week_all_start <- min(births_weeks$week_gest_num, na.rm = TRUE)
week_all_end <- max(births_weeks$week_gest_num, na.rm = TRUE)

weeks_dat <- births_weeks |>
  dplyr::filter(week_gest_num >= week_start)

week_end <- max(weeks_dat$week_gest_num, na.rm = TRUE)
message(
  "Cohorte semanal (PTB/positividad): semanas ", week_start, "–", week_end,
  " | filas: ", nrow(weeks_dat)
)
message(
  "Exposición completa: semanas ", week_all_start, "–", week_all_end,
  " | filas: ", nrow(births_weeks)
)

## Tabla por semana de gestación ----

summarise_exposure_by_week <- function(df, poll_var, prefix) {
  df |>
    dplyr::group_by(week_gest_num) |>
    dplyr::summarise(
      mean   = mean(.data[[poll_var]], na.rm = TRUE),
      median = stats::median(.data[[poll_var]], na.rm = TRUE),
      min    = min(.data[[poll_var]], na.rm = TRUE),
      max    = max(.data[[poll_var]], na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::rename(
      !!paste0(prefix, "_mean")   := mean,
      !!paste0(prefix, "_median") := median,
      !!paste0(prefix, "_min")    := min,
      !!paste0(prefix, "_max")    := max
    )
}

build_exposure_table <- function(df) {
  out <- df |>
    dplyr::distinct(week_gest_num) |>
    dplyr::arrange(week_gest_num)

  for (poll_var in names(pollutants)) {
    prefix <- sub("_krg$", "", poll_var)
    out <- out |>
      dplyr::left_join(
        summarise_exposure_by_week(df, poll_var, prefix),
        by = "week_gest_num"
      )
  }

  out |>
    dplyr::rename(semana_gestacion = week_gest_num) |>
    dplyr::mutate(
      dplyr::across(
        dplyr::matches("^(pm25|no2|o3)_"),
        ~ round(.x, 2)
      )
    )
}

build_week_table <- function(df, expo_df) {
  out <- df |>
    dplyr::group_by(week_gest_num) |>
    dplyr::summarise(
      n_nacimientos = dplyr::n(),
      n_preterminos = sum(birth_preterm == 1L, na.rm = TRUE),
      pct_preterminos = if (dplyr::n() > 0) {
        100 * sum(birth_preterm == 1L, na.rm = TRUE) / dplyr::n()
      } else {
        NA_real_
      },
      .groups = "drop"
    )

  for (poll_var in names(pollutants)) {
    prefix <- sub("_krg$", "", poll_var)
    out <- out |>
      dplyr::left_join(
        summarise_exposure_by_week(expo_df, poll_var, prefix),
        by = "week_gest_num"
      )
  }

  out |>
    dplyr::rename(semana_gestacion = week_gest_num) |>
    dplyr::arrange(semana_gestacion) |>
    dplyr::mutate(
      pct_preterminos = round(pct_preterminos, 2),
      dplyr::across(
        dplyr::matches("^(pm25|no2|o3)_"),
        ~ round(.x, 2)
      )
    )
}

table_by_week <- build_week_table(weeks_dat, weeks_dat)

weeks_ptb <- weeks_dat |>
  dplyr::filter(birth_preterm == 1L)

table_by_week_ptb <- build_week_table(weeks_dat, weeks_ptb)

table_exposure_all_weeks <- build_exposure_table(births_weeks)

readme_notes <- tibble::tribble(
  ~section, ~description,
  "Semanas (Por_semana)", paste0("Semanas de gestación ", week_start, " a ", week_end, "."),
  "Exposicion_por_semana", paste0(
    "Exposición Kriging por semana gestacional ", week_all_start, " a ", week_all_end,
    " (cohorte completa)."
  ),
  "n_nacimientos", "Nacimientos con observación de exposición en esa semana (una fila por persona-semana).",
  "n_preterminos", "Entre ellos, nacimientos con semanas de gestación al parto < 37.",
  "pct_preterminos", "Porcentaje de pretérminos en la semana (2 decimales).",
  "Por_semana", "Exposición (Kriging) calculada sobre toda la cohorte en cada semana (desde sem. 28).",
  "Por_semana_PTB", "Misma estructura; exposición calculada solo entre pretérminos (parto < 37 sem).",
  "Unidades", "PM2.5: ug/m3; NO2 y O3: ppbv."
)

out_xlsx <- paste0(data_out, "Table_positivity_by_week.xlsx")

writexl::write_xlsx(
  list(
    Readme = readme_notes,
    Exposicion_por_semana = table_exposure_all_weeks,
    Por_semana = table_by_week,
    Por_semana_PTB = table_by_week_ptb
  ),
  path = out_xlsx
)

message("Tabla guardada: ", out_xlsx)
print(table_exposure_all_weeks)
print(table_by_week)
print(table_by_week_ptb)

beepr::beep(8)
