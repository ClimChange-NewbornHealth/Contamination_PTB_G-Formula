# 7.2 Positivity analysis — support for g-formula interventions ----
#
# Quantifies how many preterm births (PTB) have observed exposure support under
# each cap / proportional-reduction scenario, by full pregnancy, trimester,
# and gestational weeks 28–36 (g-formula risk window).
#
# Outputs (02_Output/Descriptives/):
#   Table_positivity_PT_births.xlsx
#   Positivity_joint_expo_lag_weeks28_36.png
#   Positivity_joint_expo_lag_by_pollutant.png
#   Natural_course_coefs_weeks28_36.png

rm(list = ls())

## Settings ----
source("00_Code/0.1 Settings.R")
source("00_Code/0.2 Packages.R")

data_inp <- "01_Data/Output/"
data_out <- "02_Output/Descriptives/"
dir_gform_models <- "02_Output/G-Form/Models/"

risk_weeks <- 28:36
trimester_breaks <- list(
  t1 = 1:13,
  t2 = 14:26,
  t3 = 27:44
)

cap_values <- c(5, 10, 15, 20)
pct_values <- c(0.10, 0.20, 0.30)

pollutant_specs <- list(
  pm25 = list(
    label = "PM2.5",
    unit = "ug/m3",
    vars = c(pm25_krg = "Kriging", pm25_idw = "IDW")
  ),
  no2 = list(
    label = "NO2",
    unit = "ppbv",
    vars = c(no2_krg = "Kriging", no2_idw = "IDW")
  ),
  o3 = list(
    label = "O3",
    unit = "ppbv",
    vars = c(o3_krg = "Kriging", o3_idw = "IDW")
  )
)

all_poll_vars <- unlist(lapply(pollutant_specs, function(x) names(x$vars)))

severe_scenarios <- list(
  pm25_krg = list(type = "cap", value = 5),
  no2_krg = list(type = "cap", value = 5),
  o3_krg = list(type = "pct", value = 0.30)
)

dir.create(data_out, recursive = TRUE, showWarnings = FALSE)

format_pct <- function(x) {
  formatC(x, format = "f", digits = 2, decimal.mark = ".")
}

format_n <- function(x) {
  format(x, big.mark = " ", scientific = FALSE, trim = TRUE)
}

## Load data ----
births <- rio::import(paste0(data_inp, "births_2010_2020_exposure.RData"))
births_weeks <- rio::import(paste0(data_inp, "births_2010_2020_exposure_weeks.RData"))

glimpse(births)
glimpse(births_weeks)

if (!"birth_preterm" %in% names(births)) {
  stop("birth_preterm no encontrado en births_2010_2020_exposure.RData")
}

ptb_ids <- births |>
  dplyr::filter(.data$birth_preterm == 1L) |>
  dplyr::pull(id)

n_births <- nrow(births)
n_ptb <- length(ptb_ids)

message("Cohorte: ", n_births, " nacimientos | PTB: ", n_ptb)

births_ptb <- births |>
  dplyr::filter(.data$id %in% ptb_ids) |>
  dplyr::select(id, weeks, birth_preterm)

lagged_path <- paste0(data_inp, "births_2010_2020_exposure_weeks_lagged.RData")

if (file.exists(lagged_path)) {
  message("Cargando exposición semanal con lag precomputado...")
  load(lagged_path) # data_long
  weeks_ptb <- data_long |>
    dplyr::filter(.data$id %in% ptb_ids)
  weeks_risk_all <- data_long |>
    dplyr::filter(.data$week_gest_num %in% risk_weeks)
  rm(data_long)
} else {
  message("Archivo lagged no encontrado; calculando lag (puede tardar)...")

  compute_lagged_exposure <- function(df, idbase, pollulant, week = week_gest_num) {
    poll_name <- if (is.character(pollulant) && length(pollulant) == 1L) {
      pollulant
    } else {
      rlang::as_name(rlang::ensym(pollulant))
    }
    lag_name <- paste0(poll_name, "_lagged")

    df |>
      dplyr::mutate(
        .id_tmp = {{ idbase }},
        .poll_tmp = .data[[poll_name]],
        .week_tmp = {{ week }}
      ) |>
      dplyr::arrange(.id_tmp, .week_tmp) |>
      dplyr::group_by(.id_tmp) |>
      dplyr::mutate(
        .lag_tmp = purrr::map_dbl(dplyr::row_number(), function(i) {
          if (is.na(.week_tmp[i]) || .week_tmp[i] == 0) return(NA_real_)
          past_rows <- which(.week_tmp < .week_tmp[i])
          if (length(past_rows) == 0) return(NA_real_)
          weights <- 1 / (.week_tmp[i] - .week_tmp[past_rows])
          exposures <- .poll_tmp[past_rows]
          sum(weights * exposures, na.rm = TRUE)
        })
      ) |>
      dplyr::ungroup() |>
      dplyr::select(-.id_tmp, -.poll_tmp, -.week_tmp) |>
      dplyr::rename(!!lag_name := .lag_tmp)
  }

  add_lags_for_pollutants <- function(df, poll_vars) {
    purrr::reduce(
      poll_vars,
      .init = df,
      .f = function(acc, v) {
        out <- compute_lagged_exposure(acc, idbase = id, pollulant = v, week = week_gest_num)
        lag_col <- paste0(v, "_lagged")
        dplyr::bind_cols(acc, dplyr::select(out, dplyr::all_of(lag_col)))
      }
    )
  }

  weeks_ptb <- births_weeks |>
    dplyr::filter(.data$id %in% ptb_ids) |>
    dplyr::left_join(births_ptb, by = "id")
  weeks_ptb <- add_lags_for_pollutants(weeks_ptb, all_poll_vars)

  weeks_risk_all <- births_weeks |>
    dplyr::filter(.data$week_gest_num %in% risk_weeks)
  weeks_risk_all <- add_lags_for_pollutants(weeks_risk_all, names(severe_scenarios))
}

rm(births_weeks)
gc()

weeks_ptb <- weeks_ptb |>
  dplyr::mutate(
    trimester = dplyr::case_when(
      .data$week_gest_num <= 13L ~ "t1",
      .data$week_gest_num <= 26L ~ "t2",
      TRUE ~ "t3"
    )
  )

## Helpers: positividad vectorizada ----

summarise_support <- function(in_support_vec, n_den = n_ptb) {
  n_in <- sum(in_support_vec, na.rm = TRUE)
  n_out <- n_den - n_in
  pct_in <- if (n_den > 0) round(100 * n_in / n_den, 2) else NA_real_
  list(n_in = n_in, n_out = n_out, pct_in = pct_in)
}

person_max_exposure <- function(df, poll_col, week_idx = NULL) {
  dat <- df
  if (!is.null(week_idx)) {
    dat <- dat |> dplyr::filter(.data$week_gest_num %in% week_idx)
  }
  dat |>
    dplyr::group_by(.data$id) |>
    dplyr::summarise(
      max_expo = max(.data[[poll_col]], na.rm = TRUE),
      .groups = "drop"
    )
}

build_cumulative_support <- function(df, poll_col, delivery_weeks) {
  hist <- df |>
    dplyr::arrange(.data$id, .data$week_gest_num) |>
    dplyr::group_by(.data$id) |>
    dplyr::mutate(
      .expo = .data[[poll_col]],
      .cummax = cummax(replace(.expo, is.na(.expo), -Inf)),
      .cumpos = cumsum(is.finite(.expo) & .expo > 0) > 0L
    ) |>
    dplyr::ungroup()

  cap_rows <- hist |>
    dplyr::filter(.data$week_gest_num %in% risk_weeks) |>
    dplyr::select(.data$id, .data$week_gest_num, .data$.cummax)

  pct_rows <- hist |>
    dplyr::filter(.data$week_gest_num %in% risk_weeks) |>
    dplyr::select(.data$id, .data$week_gest_num, .data$.cumpos)

  eligible <- delivery_weeks |>
    dplyr::select(.data$id, .data$weeks)

  list(
    cap = cap_rows |>
      dplyr::left_join(eligible, by = "id") |>
      dplyr::filter(.data$weeks >= .data$week_gest_num),
    pct = pct_rows |>
      dplyr::left_join(eligible, by = "id") |>
      dplyr::filter(.data$weeks >= .data$week_gest_num)
  )
}


## Tabla principal de positividad (PTB) ----

build_scenario_rows <- function() {
  rows <- list()

  for (poll_key in names(pollutant_specs)) {
    spec <- pollutant_specs[[poll_key]]
    for (poll_var in names(spec$vars)) {
      method <- spec$vars[[poll_var]]

      max_full <- person_max_exposure(weeks_ptb, poll_var)
      max_t1 <- person_max_exposure(weeks_ptb, poll_var, trimester_breaks$t1)
      max_t2 <- person_max_exposure(weeks_ptb, poll_var, trimester_breaks$t2)
      max_t3 <- person_max_exposure(weeks_ptb, poll_var, trimester_breaks$t3)
      cum_support <- build_cumulative_support(weeks_ptb, poll_var, births_ptb)

      scenario_defs <- c(
        stats::setNames(
          paste0("< ", cap_values, " ", spec$unit),
          paste0("cap_", cap_values)
        ),
        stats::setNames(
          paste0("-", pct_values * 100, "%"),
          paste0("pct_", pct_values * 100)
        )
      )

      for (sc_name in names(scenario_defs)) {
        is_cap <- grepl("^cap_", sc_name)
        sc_value <- if (is_cap) {
          as.numeric(sub("^cap_", "", sc_name))
        } else {
          as.numeric(sub("^pct_", "", sc_name)) / 100
        }

        if (is_cap) {
          s_full <- summarise_support(max_full$max_expo > sc_value)
          s_t1 <- summarise_support(max_t1$max_expo > sc_value)
          s_t2 <- summarise_support(max_t2$max_expo > sc_value)
          s_t3 <- summarise_support(max_t3$max_expo > sc_value)
        } else {
          pct_flags <- is.finite(max_full$max_expo) & max_full$max_expo > 0
          s_full <- summarise_support(pct_flags)
          pct_t1 <- is.finite(max_t1$max_expo) & max_t1$max_expo > 0
          s_t1 <- summarise_support(pct_t1)
          pct_t2 <- is.finite(max_t2$max_expo) & max_t2$max_expo > 0
          s_t2 <- summarise_support(pct_t2)
          pct_t3 <- is.finite(max_t3$max_expo) & max_t3$max_expo > 0
          s_t3 <- summarise_support(pct_t3)
        }

        row <- tibble::tibble(
          pollutant = spec$label,
          method = method,
          pollutant_var = poll_var,
          scenario = scenario_defs[[sc_name]],
          scenario_type = if (is_cap) "cap" else "pct",
          scenario_value = sc_value,
          n_ptb_total = n_ptb,
          n_ptb_in_full = s_full$n_in,
          n_ptb_out_full = s_full$n_out,
          pct_ptb_in_full = format_pct(s_full$pct_in),
          n_ptb_in_t1 = s_t1$n_in,
          n_ptb_out_t1 = s_t1$n_out,
          pct_ptb_in_t1 = format_pct(s_t1$pct_in),
          n_ptb_in_t2 = s_t2$n_in,
          n_ptb_out_t2 = s_t2$n_out,
          pct_ptb_in_t2 = format_pct(s_t2$pct_in),
          n_ptb_in_t3 = s_t3$n_in,
          n_ptb_out_t3 = s_t3$n_out,
          pct_ptb_in_t3 = format_pct(s_t3$pct_in)
        )

        for (w in risk_weeks) {
          if (is_cap) {
            flags <- cum_support$cap |>
              dplyr::filter(.data$week_gest_num == w)
            n_at_risk <- nrow(flags)
            n_in <- sum(flags$.cummax > sc_value, na.rm = TRUE)
          } else {
            flags <- cum_support$pct |>
              dplyr::filter(.data$week_gest_num == w)
            n_at_risk <- nrow(flags)
            n_in <- sum(flags$.cumpos, na.rm = TRUE)
          }
          n_out <- n_at_risk - n_in
          pct_in <- if (n_at_risk > 0) round(100 * n_in / n_at_risk, 2) else NA_real_

          row[[paste0("w", w, "_n_ptb_at_risk")]] <- n_at_risk
          row[[paste0("w", w, "_n_ptb_in")]] <- n_in
          row[[paste0("w", w, "_n_ptb_out")]] <- n_out
          row[[paste0("w", w, "_pct_ptb_in")]] <- format_pct(pct_in)
        }

        rows[[length(rows) + 1L]] <- row
      }
    }
  }

  dplyr::bind_rows(rows)
}

positivity_summary <- build_scenario_rows()

## Detalle largo por semana 28–36 ----

positivity_by_week <- positivity_summary |>
  tidyr::pivot_longer(
    cols = dplyr::matches("^w(2[89]|3[0-6])_"),
    names_to = c("week", ".value"),
    names_pattern = "^w([0-9]+)_(.*)$"
  ) |>
  dplyr::mutate(week = as.integer(.data$week)) |>
  dplyr::select(
    pollutant, method, pollutant_var, scenario_type, scenario, scenario_value,
    week, n_ptb_at_risk, n_ptb_in, n_ptb_out, pct_ptb_in
  )

## Cohorte (referencia) ----

cohort_summary <- tibble::tibble(
  metric = c("All births", "Preterm births (birth_preterm = 1)"),
  n = c(n_births, n_ptb),
  pct_of_births = c(
    "100.00",
    format_pct(if (n_births > 0) 100 * n_ptb / n_births else NA_real_)
  )
)

readme_notes <- tibble::tribble(
  ~section, ~description,
  "Denominator (summary)", "All preterm births in analytic cohort (birth_preterm = 1).",
  "Cap scenarios", "PTB in support if ANY gestational week has observed exposure > cap (min(X, c) binds).",
  "Pct scenarios", "PTB in support if ANY gestational week has observed exposure > 0 (reduction applies).",
  "Full pregnancy", "Weeks 1 through last observed gestational week per birth.",
  "Trimesters", "T1: weeks 1-13; T2: 14-26; T3: 27+.",
  "Weeks 28-36", "Among PTB still at risk at week w (delivery weeks >= w); support if exposure history up to w is affected.",
  "pct columns", "Percent of PTB in support; 2 decimals; decimal separator '.'",
  "Joint plots", "Person-weeks at weeks 28-36; severe scenarios: PM2.5/NO2 cap 5, O3 -30%.",
  "Coefficients", "From cached natural_course_{pollutant}.rds (weeks 28-36); not re-estimated."
)

out_xlsx <- paste0(data_out, "Table_positivity_PT_births.xlsx")

writexl::write_xlsx(
  list(
    Readme = readme_notes,
    Cohort_summary = cohort_summary,
    Positivity_summary = positivity_summary,
    Positivity_by_week = positivity_by_week
  ),
  path = out_xlsx
)

message("Tabla guardada: ", out_xlsx)

rm(weeks_ptb)
gc()

## Gráficos: distribución conjunta expo vs lag (semanas 28–36) ----

plot_joint_one <- function(poll_var, spec, scenario) {
  lag_col <- paste0(poll_var, "_lagged")
  dat <- weeks_risk_all |>
    dplyr::filter(!is.na(.data[[poll_var]]), !is.na(.data[[lag_col]]))

  if (!nrow(dat)) return(NULL)

  if (scenario$type == "cap") {
    cap <- scenario$value
    dat <- dat |>
      dplyr::mutate(
        support = .data[[poll_var]] > cap,
        support_lab = dplyr::if_else(.data$support, "Above cap (in support)", "At/below cap")
      )
    ref_expo <- cap
    ref_note <- paste0("Cap = ", cap)
  } else {
    pct <- scenario$value
    dat <- dat |>
      dplyr::mutate(
        support = is.finite(.data[[poll_var]]) & .data[[poll_var]] > 0,
        support_lab = dplyr::if_else(.data$support, "Reduction applies", "Zero exposure")
      )
    ref_expo <- NA_real_
    ref_note <- paste0("-", pct * 100, "% reduction")
  }

  ggplot2::ggplot(
    dat,
    ggplot2::aes(
      x = .data[[poll_var]],
      y = .data[[lag_col]],
      colour = .data$support_lab
    )
  ) +
    ggplot2::geom_point(alpha = 0.12, size = 0.6) +
    ggplot2::scale_colour_manual(values = c(
      "Above cap (in support)" = "#d62728",
      "At/below cap" = "#1f77b4",
      "Reduction applies" = "#d62728",
      "Zero exposure" = "#1f77b4"
    )) +
    ggplot2::facet_wrap(~ week_gest_num, ncol = 3, labeller = ggplot2::labeller(
      week_gest_num = function(x) paste0("Week ", x)
    )) +
    ggplot2::labs(
      title = paste0(spec$label, " (", spec$vars[[poll_var]], ") — ", ref_note),
      subtitle = paste0(
        "Observed person-weeks (weeks 28–36); n = ",
        format_n(nrow(dat))
      ),
      x = paste0("Weekly exposure (", spec$unit, ")"),
      y = "Weighted lagged exposure (L)",
      colour = NULL
    ) +
    ggplot2::theme_light(base_size = 11) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      legend.position = "bottom",
      strip.text = ggplot2::element_text(size = 9)
    ) +
    {
      if (scenario$type == "cap" && is.finite(ref_expo)) {
        ggplot2::geom_vline(xintercept = ref_expo, linetype = "dashed", colour = "grey30")
      }
    }
}

poll_plots <- list()
for (poll_var in names(severe_scenarios)) {
  poll_key <- sub("_(krg|idw)$", "", poll_var)
  spec <- pollutant_specs[[poll_key]]
  p <- plot_joint_one(poll_var, spec, severe_scenarios[[poll_var]])
  if (!is.null(p)) poll_plots[[poll_var]] <- p
}

if (length(poll_plots)) {
  p_joint <- patchwork::wrap_plots(poll_plots, ncol = 1)
  ggplot2::ggsave(
    filename = paste0(data_out, "Positivity_joint_expo_lag_weeks28_36.png"),
    plot = p_joint,
    width = 14,
    height = 5 * length(poll_plots),
    dpi = 300
  )
  message("Guardado: Positivity_joint_expo_lag_weeks28_36.png")
}

make_focus_scatter <- function(poll_var, week_focus = 32L) {
  poll_key <- sub("_(krg|idw)$", "", poll_var)
  spec <- pollutant_specs[[poll_key]]
  scenario <- severe_scenarios[[poll_var]]
  lag_col <- paste0(poll_var, "_lagged")

  dat <- weeks_risk_all |>
    dplyr::filter(
      .data$week_gest_num == week_focus,
      !is.na(.data[[poll_var]]),
      !is.na(.data[[lag_col]])
    )

  if (!nrow(dat)) return(NULL)

  if (scenario$type == "cap") {
    dat <- dat |> dplyr::mutate(in_support = .data[[poll_var]] > scenario$value)
  } else {
    dat <- dat |> dplyr::mutate(in_support = .data[[poll_var]] > 0)
  }

  ggplot2::ggplot(dat, ggplot2::aes(x = .data[[poll_var]], y = .data[[lag_col]])) +
    ggplot2::geom_point(
      ggplot2::aes(colour = .data$in_support),
      alpha = 0.25,
      size = 0.8
    ) +
    ggplot2::scale_colour_manual(
      values = c("FALSE" = "#4C78A8", "TRUE" = "#E45756"),
      labels = c("FALSE" = "Not in support", "TRUE" = "In support"),
      name = NULL
    ) +
    ggplot2::labs(
      title = paste0(spec$label, " — week ", week_focus, " (severe scenario)"),
      x = paste0("exposicion_", week_focus, " (", spec$unit, ")"),
      y = paste0("exposicion_lagged_", week_focus)
    ) +
    ggplot2::theme_light() +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), legend.position = "bottom")
}

focus_plots <- lapply(names(severe_scenarios), make_focus_scatter)
focus_plots <- focus_plots[!vapply(focus_plots, is.null, logical(1))]

if (length(focus_plots)) {
  p_focus <- patchwork::wrap_plots(focus_plots, ncol = 3)
  ggplot2::ggsave(
    filename = paste0(data_out, "Positivity_joint_expo_lag_by_pollutant.png"),
    plot = p_focus,
    width = 15,
    height = 5,
    dpi = 300
  )
  message("Guardado: Positivity_joint_expo_lag_by_pollutant.png")
}

rm(weeks_risk_all)
gc()

## Coeficientes natural course (semanas 28–36) — un contaminante a la vez ----

extract_natural_course_coefs <- function(pollutant, dir_models, weeks = risk_weeks) {
  path <- file.path(dir_models, paste0("natural_course_", pollutant, ".rds"))
  if (!file.exists(path)) {
    warning("No se encontró: ", path)
    return(NULL)
  }

  message("Cargando coeficientes desde ", path, " ...")
  cached <- readRDS(path)
  on.exit({
    rm(cached)
    gc()
  }, add = TRUE)

  if (is.null(cached$model_store)) {
    warning("model_store vacío en ", path)
    return(NULL)
  }

  rows <- list()
  for (w in weeks) {
    ms <- cached$model_store[[as.character(w)]]
    if (is.null(ms) || is.null(ms$model)) next

    coefs <- stats::coef(ms$model)
    vc <- tryCatch(stats::vcov(ms$model), error = function(e) NULL)
    se <- if (!is.null(vc)) sqrt(diag(vc)) else rep(NA_real_, length(coefs))

    for (nm in names(coefs)) {
      rows[[length(rows) + 1L]] <- tibble::tibble(
        pollutant_var = pollutant,
        risk_week = w,
        term = nm,
        beta = unname(coefs[[nm]]),
        se = unname(se[[nm]]),
        hr = exp(unname(coefs[[nm]])),
        n_model = ms$model$n
      )
    }
  }

  if (!length(rows)) return(NULL)
  dplyr::bind_rows(rows)
}

krg_pollutants <- c("pm25_krg", "no2_krg", "o3_krg")
coef_tables <- vector("list", length(krg_pollutants))
names(coef_tables) <- krg_pollutants

for (i in seq_along(krg_pollutants)) {
  coef_tables[[i]] <- extract_natural_course_coefs(
    krg_pollutants[i],
    dir_models = dir_gform_models
  )
  gc()
}

coef_all <- dplyr::bind_rows(coef_tables)

coef_expo_lag <- coef_all |>
  dplyr::filter(
    grepl("^exposicion_[0-9]+$", .data$term) |
      grepl("^exposicion_lagged_[0-9]+$", .data$term)
  ) |>
  dplyr::mutate(
    beta = round(.data$beta, 6),
    se = round(.data$se, 6),
    hr = round(.data$hr, 6)
  ) |>
  dplyr::arrange(.data$pollutant_var, .data$risk_week, .data$term)

if (nrow(coef_expo_lag)) {
  coef_plot_dat <- coef_expo_lag |>
    dplyr::mutate(
      term_type = dplyr::if_else(
        grepl("^exposicion_lagged_", .data$term),
        "Lagged exposure (beta)",
        "Weekly exposure (beta)"
      ),
      pollutant = dplyr::case_when(
        grepl("^pm25", .data$pollutant_var) ~ "PM2.5",
        grepl("^no2", .data$pollutant_var) ~ "NO2",
        grepl("^o3", .data$pollutant_var) ~ "O3",
        TRUE ~ .data$pollutant_var
      )
    )

  p_coef <- ggplot2::ggplot(
    coef_plot_dat,
    ggplot2::aes(x = .data$risk_week, y = .data$beta, colour = .data$term_type)
  ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 1.8) +
    ggplot2::facet_wrap(~ pollutant, ncol = 3, scales = "free_y") +
    ggplot2::labs(
      title = "Natural-course Cox coefficients (weeks 28–36)",
      subtitle = "Cached models (natural_course_*.rds); beta on log-hazard scale",
      x = "Risk week",
      y = expression(hat(beta)),
      colour = NULL
    ) +
    ggplot2::theme_light() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      legend.position = "bottom"
    )

  ggplot2::ggsave(
    filename = paste0(data_out, "Natural_course_coefs_weeks28_36.png"),
    plot = p_coef,
    width = 12,
    height = 5,
    dpi = 300
  )
  message("Guardado: Natural_course_coefs_weeks28_36.png")
} else {
  message(
    "Coeficientes no exportados: no hay natural_course_*.rds legibles en ",
    dir_gform_models
  )
}

## Anexar coeficientes al Excel ----

writexl::write_xlsx(
  list(
    Readme = readme_notes,
    Cohort_summary = cohort_summary,
    Positivity_summary = positivity_summary,
    Positivity_by_week = positivity_by_week,
    Natural_course_coefs_28_36 = coef_expo_lag,
    Natural_course_coefs_all_terms = coef_all
  ),
  path = out_xlsx
)

message("Tabla actualizada (con coeficientes): ", out_xlsx)
print(head(positivity_summary, 10))

beepr::beep(8)
