# 11.3 Exposure–perinatal outcome models — tablas Cox ----
#
# Uso:
#   Rscript "00_Code/11.3 Exposure_PO_Tables_Models.R"
#
# Requiere: 02_Output/Exposure_PO/Models/Exposure_models_PO_cox.xlsx
# Salida:   02_Output/Exposure_PO/Tables/Tab_Exposure_PO.xlsx
#           02_Output/Exposure_PO/Tables/Tab_Exposure_PO_IQR.xlsx

source("00_Code/0.1 Settings.R")

install_load <- function(packages) {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg, repos = "https://cloud.r-project.org")
    }
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
}

install_load(c("readxl", "dplyr", "tidyr", "writexl"))
source("00_Code/11.0 Exposure_PO_functions.R")

data_out <- "02_Output/Exposure_PO/"
dir_models <- file.path(data_out, "Models")
dir_tables <- file.path(data_out, "Tables")
dir.create(dir_tables, recursive = TRUE, showWarnings = FALSE)

path_results <- file.path(dir_models, "Exposure_models_PO_cox.xlsx")
if (!file.exists(path_results)) {
  stop("Ejecute primero 11.1 Exposure_PO_Models.R")
}

models_cox <- readxl::read_excel(path_results, sheet = "cox_models")

prepare_po_table_data <- function(models_df) {
  models_df |>
    dplyr::filter(
      (model_type == "single" & tiempo == "tot") |
        model_type == "t1_t2_t3"
    ) |>
    dplyr::arrange(
      dependent_var, contaminante, tipo, adjustment,
      model_type, exposure_scale, term
    ) |>
    dplyr::group_by(
      dependent_var, contaminante, tipo, adjustment,
      model_type, exposure_scale
    ) |>
    dplyr::mutate(
      exposure = dplyr::case_when(
        model_type == "single" & tiempo == "tot" ~ "Overall",
        model_type == "t1_t2_t3" & dplyr::row_number() == 1L ~ "T1",
        model_type == "t1_t2_t3" & dplyr::row_number() == 2L ~ "T2",
        model_type == "t1_t2_t3" & dplyr::row_number() == 3L ~ "T3",
        TRUE ~ NA_character_
      )
    ) |>
    dplyr::ungroup() |>
    dplyr::filter(!is.na(.data$exposure))
}

outcomes_order <- c(
  "birth_preterm",
  "birth_very_preterm",
  "birth_moderately_preterm",
  "birth_late_preterm",
  "lbw",
  "tlbw",
  "sga"
)

outcomes_labels <- c(
  birth_preterm = "Preterm birth",
  birth_very_preterm = "Very preterm birth",
  birth_moderately_preterm = "Moderately preterm birth",
  birth_late_preterm = "Late preterm birth",
  lbw = "Low birth weight",
  tlbw = "Very low birth weight",
  sga = "Small for gestational age"
)

exposure_order <- c("Overall", "T1", "T2", "T3")

col_order <- c(
  "outcome", "exposure",
  "pm25_Unadjusted", "pm25_Adjusted",
  "no2_Unadjusted", "no2_Adjusted",
  "o3_Unadjusted", "o3_Adjusted"
)

col_labels <- c(
  outcome = "Outcome",
  exposure = "Exposure",
  pm25_Unadjusted = "PM2.5 unadjusted",
  pm25_Adjusted = "PM2.5 adjusted",
  no2_Unadjusted = "NO2 unadjusted",
  no2_Adjusted = "NO2 adjusted",
  o3_Unadjusted = "O3 unadjusted",
  o3_Adjusted = "O3 adjusted"
)

build_table_all <- function(data_all) {
  data_all <- data_all |>
    dplyr::mutate(
      exposure = factor(.data$exposure, levels = exposure_order),
      dependent_var = factor(.data$dependent_var, levels = outcomes_order),
      value_fmt = format_po_effect_ci(.data$estimate, .data$conf.low, .data$conf.high),
      col_name = paste0(.data$contaminante, "_", .data$adjustment)
    )

  data_wide <- data_all |>
    dplyr::select(
      outcome = .data$dependent_var,
      exposure = .data$exposure,
      col_name = .data$col_name,
      value_fmt = .data$value_fmt
    ) |>
    tidyr::pivot_wider(
      names_from = .data$col_name,
      values_from = .data$value_fmt
    ) |>
    dplyr::arrange(.data$outcome, .data$exposure) |>
    dplyr::mutate(
      outcome = outcomes_labels[as.character(.data$outcome)],
      exposure = as.character(.data$exposure)
    )

  col_order_avail <- col_order[col_order %in% names(data_wide)]
  data_wide <- data_wide |> dplyr::select(dplyr::all_of(col_order_avail))

  for (i in seq_along(col_labels)) {
    old_nm <- names(col_labels)[[i]]
    new_nm <- col_labels[[i]]
    if (old_nm %in% names(data_wide)) {
      names(data_wide)[names(data_wide) == old_nm] <- new_nm
    }
  }

  data_wide
}

write_scale_tables <- function(scale, suffix = "") {
  table_data_cox <- prepare_po_table_data(models_cox) |>
    dplyr::filter(.data$exposure_scale == scale)

  tab_hr <- build_table_all(table_data_cox)

  out_path <- file.path(
    dir_tables,
    if (nzchar(suffix)) paste0("Tab_Exposure_PO", suffix, ".xlsx") else "Tab_Exposure_PO.xlsx"
  )
  writexl::write_xlsx(list(HR = tab_hr), path = out_path)
  message("Guardado: ", out_path)
}

write_scale_tables("raw", suffix = "")
write_scale_tables("iqr", suffix = "_IQR")
