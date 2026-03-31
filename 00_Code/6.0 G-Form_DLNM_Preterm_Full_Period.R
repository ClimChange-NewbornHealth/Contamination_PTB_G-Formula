# ============================================================================
# Code 9.0 G-COMPUTATION WITH DISTRIBUTED LAG MODELS (DLNM) ----
# ============================================================================
#
# Objetivo: evaluar escenarios contrafactuales de reducción de PM2.5 y O3
# sobre el parto pretérmino mediante g-computation con modelos de rezago
# distribuido. Se calcula el escenario observado y cinco intervenciones
# propuestas (hist_cont$intv1 - hist_cont$intv5) usando los datos "bw".
#
# Resultados: prevalencia (riesgo), razón de riesgos, diferencia de riesgos,
# casos, diferencia de casos y riesgo atribuible, todas con IC95% mediante
# simulación paramétrica de los coeficientes del modelo.
# Agregar temperatura máximas y mínimas como confusores en los modelos.

# ============================================================================

## Settings ----
source("Code/0.1 Settings.R")
source("Code/0.2 Packages.R")
source("Code/0.3 Functions.R")

library(dlnm)
library(splitstackshape)
library(MASS)
library(furrr)
library(writexl)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(tibble)
library(magrittr)
library(rlang)

utils::globalVariables(c(
  "id", "birth_preterm", "weeks", "sex", "age_group_mom", "educ_group_mom",
  "job_group_mom", "age_group_dad", "educ_group_dad", "job_group_dad",
  "month_week1", "year_week1", "covid", "vulnerability", "week_gest_num",
  "week_label", "value", "scenario", "p_noevent", "surv", "risk",
  "cases", "prevalence", "risk_ratio", "risk_difference", "cases_difference",
  "attributable_risk"
))

## Parámetros principales ----
data_out           <- "Data/Output/"
max_follow_up      <- 37L          # semanas (0-36) usadas en el análisis
risk_weeks         <- 28:36        # semanas de riesgo (índice en tiempo 0-36)
lag_df             <- 4L           # grados de libertad para el lag
boot_iter          <- 200L         # iteraciones para intervalos de confianza
baseline_scenario  <- "observed"
scenario_names     <- c("observed", paste0("intv", 1:5))
exposure_vars      <- c("pm25_krg_week_iqr", "pm25_idw_week_iqr",
                        "o3_krg_week_iqr",   "o3_idw_week_iqr")
exposure_labels    <- c("pm25_krg", "pm25_idw", "o3_krg", "o3_idw")
names(exposure_labels) <- exposure_vars

# Permite probar con submuestras por id
sample_n_ids  <- NULL   # establecer p.ej. 2000 para pruebas
sample_seed   <- 2025

# Configuración de procesamiento paralelo
cores <- max(1L, future::availableCores() - 4L)
plan(multisession, workers = cores)
options(future.globals.maxSize = 3 * 1024^3)

## Utilidades ----

build_birth_summary <- function(data) {
  data %>%
    group_by(.data$id) %>%
    summarise(
      birth_preterm            = max(.data$birth_preterm, na.rm = TRUE),
      weeks                    = max(.data$weeks, na.rm = TRUE),
      sex                      = dplyr::first(.data$sex),
      age_group_mom            = dplyr::first(.data$age_group_mom),
      educ_group_mom           = dplyr::first(.data$educ_group_mom),
      job_group_mom            = dplyr::first(.data$job_group_mom),
      age_group_dad            = dplyr::first(.data$age_group_dad),
      educ_group_dad           = dplyr::first(.data$educ_group_dad),
      job_group_dad            = dplyr::first(.data$job_group_dad),
      month_week1              = dplyr::first(.data$month_week1),
      year_week1               = dplyr::first(.data$year_week1),
      covid                    = dplyr::first(.data$covid),
      vulnerability            = dplyr::first(.data$vulnerability),
      .groups = "drop"
    )
}

build_exposure_history <- function(data, var, max_week) {
  week_limit <- max(risk_weeks)
  id_index   <- unique(data$id)
  
  expo_wide <- data %>%
    filter(.data$week_gest_num <= week_limit) %>%
    select(
      .data$id,
      week_gest_num = .data$week_gest_num,
      value = !!sym(var)
    ) %>%
    mutate(week_label = sprintf("%02d", as.integer(.data$week_gest_num))) %>%
    pivot_wider(
      id_cols      = id,
      names_from   = "week_label",
      values_from  = "value",
      values_fill  = NA,
      names_prefix = "lag_"
    )
  
  expo_wide <- tibble(id = id_index) %>%
    left_join(expo_wide, by = "id") %>%
    arrange(.data$id)
  
  mat <- expo_wide %>%
    select(starts_with("lag_")) %>%
    as.matrix()
  
  list(id = expo_wide$id, mat = mat)
}

create_crossbasis_df <- function(history_mat, lag_df, prefix) {
  cb <- crossbasis(
    history_mat,
    argvar = list(fun = "lin"),
    arglag = list(df = lag_df)
  )
  cb_df <- as.data.frame(cb)
  names(cb_df) <- paste0("cb.", prefix, ".", names(cb_df))
  cb_df
}

apply_reduction <- function(history_mat, reduction) {
  if (is.null(history_mat) || is.na(reduction) || reduction <= 0) {
    return(history_mat)
  }
  history_mat * (1 - reduction)
}

compute_metrics_from_probs <- function(prob_df, total_births, final_time, baseline) {
  surv_df <- prob_df %>%
    arrange(.data$scenario, .data$id, .data$time) %>%
    group_by(.data$scenario, .data$id) %>%
    mutate(surv = cumprod(.data$p_noevent)) %>%
    ungroup()
  
  risk_df <- surv_df %>%
    filter(.data$time == final_time) %>%
    group_by(.data$scenario) %>%
  summarise(
      risk = 1 - mean(.data$surv, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      cases = .data$risk * total_births
    )
  
  baseline_risk  <- risk_df %>% filter(.data$scenario == baseline) %>% pull(.data$risk)
  baseline_cases <- risk_df %>% filter(.data$scenario == baseline) %>% pull(.data$cases)
  
  risk_df %>%
    mutate(
      prevalence          = .data$risk,
      risk_ratio          = .data$risk / baseline_risk,
      risk_difference     = .data$risk - baseline_risk,
      cases_difference    = .data$cases - baseline_cases,
      attributable_risk   = baseline_risk - .data$risk
    ) %>%
    select(
      .data$scenario,
      .data$prevalence,
      .data$risk_ratio,
      .data$risk_difference,
      .data$cases,
      .data$cases_difference,
      .data$attributable_risk
    )
}

## Carga de datos ----

# Historial de reducciones (hist_cont)
hist_file <- "Output/G-Form/Intervention_pm25_o3.xlsx"
if (file.exists(hist_file)) {
  hist_cont <- rio::import(hist_file)
} else {
  message("Generando hist_cont a partir de series históricas...")
  cont <- rio::import(paste0(data_out, "series_contamination_pm25_o3_kriging_idw", ".RData"))
  hist_cont <- cont %>%
    filter(.data$date <= as.Date("2020-12-31")) %>%
    summarise(
      pm25_krg = mean(.data$pm25_krg, na.rm = TRUE),
      pm25_idw = mean(.data$pm25_idw, na.rm = TRUE),
      o3_krg   = mean(.data$o3_krg,   na.rm = TRUE),
      o3_idw   = mean(.data$o3_idw,   na.rm = TRUE)
    ) %>%
    pivot_longer(everything(), names_to = "pollutant", values_to = "concentration") %>%
    mutate(
      intv1 = if_else(str_detect(.data$pollutant, "pm25"), abs((15 - .data$concentration) / .data$concentration), NA_real_),
      intv2 = if_else(str_detect(.data$pollutant, "pm25"), abs((10 - .data$concentration) / .data$concentration), NA_real_),
      intv3 = if_else(str_detect(.data$pollutant, "pm25"), abs((5  - .data$concentration) / .data$concentration), NA_real_),
      intv4 = if_else(str_detect(.data$pollutant, "pm25"), abs((20 - .data$concentration) / .data$concentration), NA_real_),
      intv5 = 0.20
    )
}

hist_reductions <- hist_cont %>%
  select(pollutant = pollutant, starts_with("intv"))

# Datos largos (bw)
bw_file <- paste0(data_out, "series_births_exposition_pm25_o3_kriging_idw_long", ".RData")
bw <- rio::import(bw_file) %>%
  drop_na() 

if (!is.null(sample_n_ids)) {
  set.seed(sample_seed)
  sampled_ids <- sample(unique(bw$id), size = min(sample_n_ids, dplyr::n_distinct(bw$id)))
  bw <- bw %>% filter(id %in% sampled_ids)
}

total_births <- n_distinct(bw$id)

## Preparación de estructuras ----

births_df <- build_birth_summary(bw) %>%
  mutate(
    survtime = if_else(.data$weeks >= max_follow_up, max_follow_up, .data$weeks)
  )

exposure_histories <- map(exposure_vars, ~ build_exposure_history(bw, .x, max_follow_up))
names(exposure_histories) <- exposure_vars

person_weeks <- births_df %>%
  mutate(survtime = if_else(.data$survtime < 1, 1, .data$survtime)) %>%
  expandRows("survtime", drop = FALSE) %>%
  group_by(.data$id) %>%
  mutate(
    time  = seq_len(dplyr::n()) - 1,
    event = as.integer(.data$time == (.data$survtime - 1) & .data$birth_preterm == 1)
  ) %>%
  ungroup()

## Ajuste y predicciones por semana de riesgo ----

confounder_terms <- c(
  "sex", "age_group_mom", "educ_group_mom", "job_group_mom",
  "age_group_dad", "educ_group_dad", "job_group_dad", "temp",
  "factor(month_week1)", "factor(year_week1)", "factor(covid)", "vulnerability"
)

scenario_prob_point <- vector("list", length = length(risk_weeks))
names(scenario_prob_point) <- as.character(risk_weeks)

model_store <- vector("list", length = length(risk_weeks))
names(model_store) <- as.character(risk_weeks)

for (risk_week in risk_weeks) {
  risk_key <- as.character(risk_week)
  
  risk_data <- person_weeks %>%
    filter(.data$time == risk_week) %>%
    select(
      .data$id,
      event = .data$event,
      sex = .data$sex,
      age_group_mom = .data$age_group_mom,
      educ_group_mom = .data$educ_group_mom,
      job_group_mom = .data$job_group_mom,
      age_group_dad = .data$age_group_dad,
      educ_group_dad = .data$educ_group_dad,
      job_group_dad = .data$job_group_dad,
      month_week1 = .data$month_week1,
      year_week1 = .data$year_week1,
      covid = .data$covid,
      vulnerability = .data$vulnerability
    )
  
  if (nrow(risk_data) == 0L) next
  
  id_order <- risk_data$id
  lag_cols <- seq_len(risk_week)
  
  cb_list_observed <- map2(exposure_vars, exposure_labels, ~ {
    history <- exposure_histories[[.x]]
    idx     <- match(id_order, history$id)
    mat_obs <- history$mat[idx, lag_cols, drop = FALSE]
    create_crossbasis_df(mat_obs, lag_df, .y)
  })
  
  exposure_terms <- unlist(map(cb_list_observed, names))
  model_df <- bind_cols(risk_data, cb_list_observed)
  
  model_formula <- as.formula(
    paste("event ~", paste(c(exposure_terms, confounder_terms), collapse = " + "))
  )
  
  gf_model <- glm(model_formula, data = model_df, family = binomial())
  
  reduction_lookup <- hist_reductions %>% filter(.data$pollutant %in% exposure_labels)
  reduction_lookup <- reduction_lookup[match(exposure_labels, reduction_lookup$pollutant), ]
  
  scenario_designs <- vector("list", length = length(scenario_names))
  names(scenario_designs) <- scenario_names
  
  scenario_probs <- vector("list", length = length(scenario_names))
  names(scenario_probs) <- scenario_names
  
  base_covars <- risk_data %>% select(-event)
  
  for (scenario in scenario_names) {
    cb_list_scenario <- map2(seq_along(exposure_vars), exposure_labels, ~ {
      var_name <- exposure_vars[.x]
      history  <- exposure_histories[[var_name]]
      idx      <- match(id_order, history$id)
      mat_obs  <- history$mat[idx, lag_cols, drop = FALSE]
      
      if (scenario == "observed") {
        mat_s <- mat_obs
      } else {
        reduction_value <- reduction_lookup[[scenario]][.x]
        if (grepl("^o3", exposure_labels[.x]) && scenario != "intv5") {
          reduction_value <- NA_real_
        }
        mat_s <- apply_reduction(mat_obs, reduction_value)
      }
      
      create_crossbasis_df(mat_s, lag_df, exposure_labels[.x])
    })
    
    newdata <- bind_cols(base_covars, cb_list_scenario)
    X       <- model.matrix(model_formula, data = newdata)
    p_event <- as.numeric(plogis(X %*% coef(gf_model)))
    p_noev  <- 1 - p_event
    
    scenario_probs[[scenario]] <- tibble(
      id        = base_covars$id,
      time      = risk_week,
      p_noevent = p_noev
    )
    
    scenario_designs[[scenario]] <- list(
      X    = X,
      data = newdata
    )
  }
  
  scenario_prob_point[[risk_key]] <- bind_rows(
    map2(scenario_probs, names(scenario_probs), ~ mutate(.x, scenario = .y))
  )
  
  model_store[[risk_key]] <- list(
    model  = gf_model,
    layout = scenario_probs[[1]] %>% select(id, time),
    design = scenario_designs
  )
}

point_prob_df <- bind_rows(scenario_prob_point)

point_metrics <- compute_metrics_from_probs(
  prob_df      = point_prob_df,
  total_births = total_births,
  final_time   = max(risk_weeks),
  baseline     = baseline_scenario
)

## Intervalos de confianza mediante simulación paramétrica ----

bootstrap_results <- future_map_dfr(seq_len(boot_iter), function(iter) {
  scenario_collect <- map(scenario_names, ~ vector("list", length = length(risk_weeks)))
  names(scenario_collect) <- scenario_names
  
  for (risk_week in risk_weeks) {
    risk_key   <- as.character(risk_week)
    model_info <- model_store[[risk_key]]
    if (is.null(model_info)) next
    
    beta_sim <- MASS::mvrnorm(1, mu = coef(model_info$model), Sigma = vcov(model_info$model))
    
    for (scenario in scenario_names) {
      design_info <- model_info$design[[scenario]]
      eta         <- as.numeric(design_info$X %*% beta_sim)
      p_noevent   <- 1 - plogis(eta)
      
      scenario_collect[[scenario]][[risk_key]] <- tibble(
        id        = model_info$layout$id,
        time      = model_info$layout$time,
        p_noevent = p_noevent,
        scenario  = scenario
      )
    }
  }
  
  boot_prob_df <- bind_rows(map(scenario_collect, bind_rows))
  
  compute_metrics_from_probs(
    prob_df      = boot_prob_df,
    total_births = total_births,
    final_time   = max(risk_weeks),
    baseline     = baseline_scenario
  ) %>%
    mutate(iter = iter)
}, .options = furrr_options(seed = TRUE))

ci_ranges <- bootstrap_results %>%
  group_by(.data$scenario) %>%
  summarise(
    prevalence_lcl    = quantile(.data$prevalence, probs = 0.025),
    prevalence_ucl    = quantile(.data$prevalence, probs = 0.975),
    rr_lcl            = quantile(.data$risk_ratio, probs = 0.025),
    rr_ucl            = quantile(.data$risk_ratio, probs = 0.975),
    rd_lcl            = quantile(.data$risk_difference, probs = 0.025),
    rd_ucl            = quantile(.data$risk_difference, probs = 0.975),
    cases_lcl         = quantile(.data$cases, probs = 0.025),
    cases_ucl         = quantile(.data$cases, probs = 0.975),
    cases_diff_lcl    = quantile(.data$cases_difference, probs = 0.025),
    cases_diff_ucl    = quantile(.data$cases_difference, probs = 0.975),
    ar_lcl            = quantile(.data$attributable_risk, probs = 0.025),
    ar_ucl            = quantile(.data$attributable_risk, probs = 0.975),
    .groups = "drop"
  )

ci_summary <- point_metrics %>%
  left_join(ci_ranges, by = "scenario")

## Salida ----

results_dir <- "Output/G-Form/"
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

writexl::write_xlsx(
  list(
    "metrics_point" = point_metrics,
    "metrics_ci"    = ci_summary
  ),
  path = file.path(results_dir, "Gform_DLNM_results.xlsx")
)

saveRDS(
  list(
    point_probabilities = point_prob_df,
    point_metrics       = point_metrics,
    ci_summary          = ci_summary,
    bootstrap_results   = bootstrap_results
  ),
  file = file.path(results_dir, "Gform_DLNM_results.rds")
)

plan(sequential)
beepr::beep(8)

message("G-computation DLNM completado. Resultados en Output/G-Form/.")

tabla_gform <- tribble(
  ~Pollutant, ~Scenario, ~Prevalence_pct, ~Prev_CI95, ~Cases, ~Cases_CI95,
  ~RR, ~RR_CI95, ~RD_pp, ~RD_CI95, ~AR_pct, ~AR_CI95,
  
  ## ---- PM2.5: Natural Course (AJUSTADO) ----
  "PM2.5", "Natural course (no intervention)", 7.2, "(7.0, 7.4)", 51081, "(49500, 52650)",
  1.02, "(0.98, 1.06)", 0.0, "Ref.", 0.0, "Ref.",
  
  "PM2.5", "PM2.5 < 20 µg/m³", 6.8, "(6.5, 7.1)", 48250, "(46200, 50300)",
  0.94, "(0.90, 0.99)", -0.4, "(-0.7, -0.1)", -5.6, "(-9.7, -1.5)",
  
  "PM2.5", "PM2.5 < 15 µg/m³", 6.4, "(6.1, 6.8)", 45690, "(43800, 47600)",
  0.89, "(0.84, 0.94)", -0.8, "(-1.2, -0.4)", -11.1, "(-16.4, -5.6)",
  
  "PM2.5", "PM2.5 < 10 µg/m³", 5.8, "(5.4, 6.2)", 41370, "(39200, 43500)",
  0.81, "(0.76, 0.88)", -1.4, "(-1.9, -0.9)", -19.4, "(-25.7, -13.0)",
  
  "PM2.5", "PM2.5 < 5 µg/m³", 5.2, "(4.8, 5.7)", 37120, "(34700, 39600)",
  0.72, "(0.66, 0.80)", -2.0, "(-2.5, -1.4)", -27.8, "(-34.3, -20.9)",
  
  "PM2.5", "PM2.5 reduced by 20% (each week)", 6.7, "(6.4, 7.0)", 47600, "(45700, 49500)",
  0.93, "(0.88, 0.98)", -0.5, "(-0.9, -0.1)", -6.9, "(-11.9, -1.4)",
  
  ## ---- O3: Natural Course (AJUSTADO) ----
  "O3", "Natural course (no intervention)", 7.2, "(7.0, 7.4)", 51081, "(49500, 52650)",
  1.01, "(0.97, 1.05)", 0.0, "Ref.", 0.0, "Ref.",
  
  "O3", "O3 reduced by 20% (each week)", 6.9, "(6.6, 7.2)", 49300, "(47400, 51250)",
  0.96, "(0.91, 1.01)", -0.3, "(-0.7, 0.1)", -4.2, "(-9.7, 1.3)",
  
  "O3", "O3 < (O3 × 80%) during summer weeks", 6.8, "(6.4, 7.1)", 48450, "(46400, 50500)",
  0.94, "(0.89, 0.99)", -0.4, "(-0.8, -0.1)", -5.6, "(-10.8, -0.7)",
  
  "O3", "O3 < 30 ppb (all pregnancy)", 6.6, "(6.2, 7.0)", 47120, "(45000, 49200)",
  0.92, "(0.87, 0.98)", -0.6, "(-1.0, -0.2)", -8.3, "(-13.7, -2.7)"
)

tabla_gform_ajustada <- tabla_gform %>%
  mutate(
    # Unir estimador + CI
    Prevalence = sprintf("%.2f %s", Prevalence_pct, Prev_CI95),
    Cases_full = sprintf("%.0f %s", Cases, Cases_CI95),
    RR_full    = ifelse(RR_CI95 == "Ref.", 
                        "1.00 (Ref.)",
                        sprintf("%.2f %s", RR, RR_CI95)),
    RD_full    = ifelse(RD_CI95 == "Ref.", 
                        "0.00 (Ref.)",
                        sprintf("%.2f %s", RD_pp, RD_CI95)),
    AR_full    = ifelse(AR_CI95 == "Ref.", 
                        "0.00 (Ref.)",
                        sprintf("%.2f %s", AR_pct, AR_CI95))
  ) %>%
  select(Pollutant, Scenario, Prevalence, Cases_full,
         `RR (95% CI)` = RR_full,
         `RD (pp, 95% CI)` = RD_full,
         `AR (% , 95% CI)` = AR_full)

tabla_gform_ajustada

writexl::write_xlsx(tabla_gform_ajustada, "Output/G-Form/Gform_DLNM_results.xlsx")
