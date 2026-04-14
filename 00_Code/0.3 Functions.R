# Functions ----

## Descriptives ----

descriptives <- function(x, data){
  data %>% 
    dplyr::select({{ x }}) %>%    
    #drop_na() %>% 
    summarise(Media_Prop = round(mean({{ x }}, na.rm = TRUE), 3),
              SD = round(sd({{ x }}, na.rm = TRUE), 3),
              Min = min({{ x }}, na.rm = TRUE),
              P5 = round(quantile({{ x }}, probs = 0.05, na.rm = TRUE), 3),
              P10 = round(quantile({{ x }}, probs = 0.1, na.rm = TRUE), 3),
              P25 = round(quantile({{ x }}, probs = 0.25, na.rm = TRUE), 3),
              P50 = round(quantile({{ x }}, probs = 0.50, na.rm = TRUE), 3), # Mediana
              P75 = round(quantile({{ x }}, probs = 0.75, na.rm = TRUE), 3),
              P90 = round(quantile({{ x }}, probs = 0.9, na.rm = TRUE), 3),
              P95 = round(quantile({{ x }}, probs = 0.95, na.rm = TRUE), 3),
              Max = max({{ x }}, na.rm = TRUE),
              N = n(),
              Missing = sum(is.na({{ x }})),
              Pct_miss = round(Missing/N, 4)*100
    ) %>% 
    mutate(Variable={{i}}) %>% 
    relocate(Variable)
}

## Weeks expositions functions ----

# Divide data in parts 

parts <- function(data, path, folder, num_parts = 20) {
  # Size parts
  part_size <- ceiling(nrow(data) / num_parts)
  
  # Process parts
  for (part_id in 1:num_parts) {
    # Parts
    part_data <- data[((part_id - 1) * part_size + 1):min(part_id * part_size, nrow(data)), ]
    
    # Save
    save(part_data, file=sprintf(paste0(path, folder, "/part_%02d_results.RData"), part_id))
  }
}

# Calculate contamination metrics

calculate_cont_stats <- function(row, cont_data) {
  # data.table objects
  setDT(row)
  setDT(cont_data)
  
  # Adjust dates
  row_copy <- copy(row)
  row_copy[, date_start_week := as.Date(date_start_week)]
  row_copy[, date_end_week := as.Date(date_end_week)]
  
  # Filter
  week_cont <- cont_data[date >= row_copy$date_start_week[1] & date < row_copy$date_end_week[1] & name_com == row_copy$name_com[1]]
  
  if (nrow(week_cont) == 0) {
    return(data.table(row_copy, 
                      pm25_krg_week=NA_real_, 
                      o3_krg_week=NA_real_, 
                      pm25_idw_week=NA_real_, 
                      o3_idw_week=NA_real_))
  } else {
    return(data.table(row_copy, 
                      pm25_krg_week = round(mean(as.numeric(week_cont$pm25_krg), na.rm = TRUE), 3),
                      o3_krg_week = round(mean(as.numeric(week_cont$o3_krg), na.rm = TRUE), 3),
                      pm25_idw_week = round(mean(as.numeric(week_cont$pm25_idw), na.rm = TRUE), 3),
                      o3_idw_week = round(mean(as.numeric(week_cont$o3_idw), na.rm = TRUE), 3)
                      )
                      )
  }
}

# Process data (old)

process_files <- function(input_directory, output_directory, cont_data, calc_func) {
  
  # Assure cont_data is available
  setDT(cont_data)
  
  files <- list.files(path = input_directory, full.names = TRUE, pattern = "\\.RData$")
  file_count <- 0
  
  for (file_path in files) {
    
    start <- Sys.time()
    
    file_count <- file_count + 1
    load(file_path)
    setDT(part_data)
    
    # Apply calculation function to each row as a data.table slice
    results <- part_data[, calc_func(.SD, cont_data), by = .I]
    
    # Save results
    save(results, file = file.path(output_directory, sprintf("%s_processed.RData", tools::file_path_sans_ext(basename(file_path)))))
    
    end <- Sys.time()
    cat("Time process data:", end-start, "\n") 
    
    if (file_count %% 5 == 0) {
      cat("Pause for 2 seconds to avoid overload...\n")
      Sys.sleep(2)
    }
    
  }
  
  cat("All files have been processed and saved in:", output_directory, "\n")
}

process_files_parallel <- function(input_directory,
                                   output_directory,
                                   cont_data,
                                   calc_func,
                                   workers = parallel::detectCores() - 4) {
  # data.table to fast
  setDT(cont_data)
  
  # List files RData
  files <- list.files(path       = input_directory,
                      pattern    = "\\.RData$",
                      full.names = TRUE)
  
  # Edit settings parallel
  plan(multisession, workers = workers)
  
  # Parallel process
  furrr::future_walk(
    .x = files,
    .f = function(file_path, cont_data, calc_func, output_directory) {
      # 1) Load .RData (define part_data)
      load(file_path)
      setDT(part_data)
      
      # 2) Apply functions
      results <- part_data[
        , calc_func(.SD, cont_data),
        by = .I
      ]
      
      # 3) Save results
      out_file <- file.path(
        output_directory,
        paste0(tools::file_path_sans_ext(basename(file_path)),
               "_processed.RData")
      )
      save(results, file=out_file)
    },
    # Object by worker
    cont_data        = cont_data,
    calc_func        = calc_func,
    output_directory = output_directory,
    .progress        = TRUE
  )
  
  # Volvemos al plan secuencial
  plan(sequential)
  
  beep(1)
  cat("✔ Todos los archivos procesados en paralelo y guardados en:\n", 
      output_directory, "\n")
}

# Post process load 
load_and_extract_df <- function(file_path) {
  e <- new.env()  # Enviroment to load data 
  load(file_path, envir = e)  # Load data
  # DF unique in enviroment 
  df <- e[[names(e)[1]]]  # First element in object 
  return(df)
}

# Row-wise mean exposure in gestation week window (for chunked parallel pipeline)
calculate_gest_window_means <- function(row, cont_data) {
  setDT(row)
  setDT(cont_data)
  row_copy <- copy(row)
  s <- as.Date(row_copy$date_start_week[1])
  e <- as.Date(row_copy$date_end_week[1])
  nm <- row_copy$name_com[1]
  w <- cont_data[date >= s & date <= e & name_com == nm]
  if (nrow(w) == 0) {
    row_copy[, `:=`(
      pm25_krg = NA_real_, o3_krg = NA_real_, no2_krg = NA_real_,
      pm25_idw = NA_real_, o3_idw = NA_real_, no2_idw = NA_real_,
      tad = NA_real_, ndvi = NA_real_
    )]
  } else {
    row_copy[, `:=`(
      pm25_krg = round(mean(as.numeric(w$pm25_krg), na.rm = TRUE), 3),
      o3_krg = round(mean(as.numeric(w$o3_krg), na.rm = TRUE), 3),
      no2_krg = round(mean(as.numeric(w$no2_krg), na.rm = TRUE), 3),
      pm25_idw = round(mean(as.numeric(w$pm25_idw), na.rm = TRUE), 3),
      o3_idw = round(mean(as.numeric(w$o3_idw), na.rm = TRUE), 3),
      no2_idw = round(mean(as.numeric(w$no2_idw), na.rm = TRUE), 3),
      tad = round(mean(as.numeric(w$TAD), na.rm = TRUE), 3),
      ndvi = round(mean(as.numeric(w$ndvi), na.rm = TRUE), 3)
    )]
  }
  row_copy
}

combine_processed_parts <- function(output_directory) {
  files <- list.files(output_directory, pattern = "_processed\\.RData$", full.names = TRUE)
  dfs <- lapply(files, load_and_extract_df)
  data.table::rbindlist(dfs, fill = TRUE)
}


# Calculate contamination metrics (version 2: observerd metrics )

calculate_cont_stats <- function(row, cont_data) {
  # data.table objects
  setDT(row)
  setDT(cont_data)
  
  # Adjust dates
  row_copy <- copy(row)
  row_copy[, date_start_week := as.Date(date_start_week)]
  row_copy[, date_end_week := as.Date(date_end_week)]
  
  # Filter
  week_cont <- cont_data[date >= row_copy$date_start_week[1] & date < row_copy$date_end_week[1] & name_com == row_copy$name_com[1]]
  
  if (nrow(week_cont) == 0) {
    return(data.table(row_copy, 
                      pm25_week=NA_real_, 
                      o3_week=NA_real_))
  } else {
    return(data.table(row_copy, 
                      pm25_week = round(mean(as.numeric(week_cont$daily_pm25), na.rm = TRUE), 3),
                      o3_week = round(mean(as.numeric(week_cont$daily_o3), na.rm = TRUE), 3)
                      )
                      )
  }
}


# Optimal use power computation 
options(future.globals.maxSize = 3000 * 1024^2)  
plan(multisession, workers = detectCores() - 4) # Parallelization


# Cox proportional hazards models function (semiparametric approach) -----

# time_start: nombre de columna con inicio de riesgo (entrada retardada); NULL = Surv(t_stop, evento)
fit_cox_model <- function(
    dependent, predictor, tiempo, contaminante, tipo,
    model_type, data, time_var = "edad_gest",
    time_start = NULL,
    conf.level = 0.95, adjustment = "Adjusted") {

  # Extract individual predictors list
  if (model_type == "single") {
    predictors_list <- predictor
  } else {
    predictors_list <- trimws(stringr::str_split(predictor, " \\+ ")[[1]])
  }

  # Verify all predictors exist in data
  missing_predictors <- predictors_list[!predictors_list %in% names(data)]
  if (length(missing_predictors) > 0) {
    return(data.frame(
      term = predictor,
      estimate = NA_real_,
      conf.low = NA_real_,
      conf.high = NA_real_,
      log_hr = NA_real_,
      log_hr_conf.low = NA_real_,
      log_hr_conf.high = NA_real_,
      std.error = NA_real_,
      statistic = NA_real_,
      p.value = NA_real_,
      dependent_var = dependent,
      predictor = predictor,
      tiempo = tiempo,
      contaminante = contaminante,
      tipo = tipo,
      model_type = model_type,
      adjustment = adjustment,
      n = 0
    ))
  }

  # Filter data with valid values in dependent, time variable, and all predictors
  data_subset <- data |>
    dplyr::filter(!is.na(.data[[dependent]]), !is.na(.data[[time_var]]))

  use_delayed <- !is.null(time_start) &&
    is.character(time_start) && nzchar(time_start) && time_start %in% names(data)
  if (use_delayed) {
    data_subset <- data_subset |>
      dplyr::filter(!is.na(.data[[time_start]])) |>
      dplyr::filter(.data[[time_start]] < .data[[time_var]])
  }

  for (pred in predictors_list) {
    data_subset <- data_subset |>
      dplyr::filter(!is.na(.data[[pred]]))
  }

  # If insufficient data, return NA
  if (nrow(data_subset) < 10) {
    return(data.frame(
      term = predictor,
      estimate = NA_real_,
      conf.low = NA_real_,
      conf.high = NA_real_,
      log_hr = NA_real_,
      log_hr_conf.low = NA_real_,
      log_hr_conf.high = NA_real_,
      std.error = NA_real_,
      statistic = NA_real_,
      p.value = NA_real_,
      dependent_var = dependent,
      predictor = predictor,
      tiempo = tiempo,
      contaminante = contaminante,
      tipo = tipo,
      model_type = model_type,
      adjustment = adjustment,
      n = nrow(data_subset)
    ))
  }

  # Build formula according to adjustment
  if (identical(adjustment, "Adjusted")) {
    available_controls <- control_vars[control_vars %in% names(data_subset)]

    rhs <- if (length(available_controls) > 0) {
      paste(
        paste(predictors_list, collapse = " + "),
        paste("+", paste(available_controls, collapse = " + "))
      )
    } else {
      paste(predictors_list, collapse = " + ")
    }
  } else {
    rhs <- paste(predictors_list, collapse = " + ")
  }

  surv_lhs <- if (use_delayed) {
    paste0("Surv(", time_start, ", ", time_var, ", ", dependent, ")")
  } else {
    paste0("Surv(", time_var, ", ", dependent, ")")
  }
  fml <- stats::as.formula(paste0(surv_lhs, " ~ ", rhs))

  # Fit Cox model
  model_fit <- tryCatch({
    survival::coxph(fml, data = data_subset)
  }, error = function(e) {
    return(NULL)
  })

  if (is.null(model_fit)) {
    return(data.frame(
      term = predictor,
      estimate = NA_real_,
      conf.low = NA_real_,
      conf.high = NA_real_,
      log_hr = NA_real_,
      log_hr_conf.low = NA_real_,
      log_hr_conf.high = NA_real_,
      std.error = NA_real_,
      statistic = NA_real_,
      p.value = NA_real_,
      dependent_var = dependent,
      predictor = predictor,
      tiempo = tiempo,
      contaminante = contaminante,
      tipo = tipo,
      model_type = model_type,
      adjustment = adjustment,
      n = nrow(data_subset)
    ))
  }

  # Extract results (HR when exponentiate = TRUE)
  tbl <- broom::tidy(model_fit, exponentiate = TRUE, conf.int = TRUE, conf.level = conf.level)

  # Filter only exposure terms (predictors)
  tbl_exposure <- tbl[tbl$term %in% predictors_list, ]

  if (nrow(tbl_exposure) > 0) {
    tbl_exposure <- tbl_exposure |>
      dplyr::mutate(
        hr = estimate,
        hr_conf.low = conf.low,
        hr_conf.high = conf.high,
        log_hr = log(estimate),
        log_hr_conf.low = log(conf.low),
        log_hr_conf.high = log(conf.high),
        dependent_var = dependent,
        predictor = predictor,
        tiempo = tiempo,
        contaminante = contaminante,
        tipo = tipo,
        model_type = model_type,
        adjustment = adjustment,
        n = nrow(data_subset)
      ) |>
      dplyr::select(term, estimate, conf.low, conf.high, log_hr, log_hr_conf.low, log_hr_conf.high,
                    std.error, statistic, p.value, dependent_var, predictor,
                    tiempo, contaminante, tipo, model_type, adjustment, n)
  } else {
    tbl_exposure <- data.frame(
      term = predictors_list[1],
      estimate = NA_real_,
      conf.low = NA_real_,
      conf.high = NA_real_,
      log_hr = NA_real_,
      log_hr_conf.low = NA_real_,
      log_hr_conf.high = NA_real_,
      std.error = NA_real_,
      statistic = NA_real_,
      p.value = NA_real_,
      dependent_var = dependent,
      predictor = predictor,
      tiempo = tiempo,
      contaminante = contaminante,
      tipo = tipo,
      model_type = model_type,
      adjustment = adjustment,
      n = nrow(data_subset)
    )
  }

  rm(model_fit); gc()

  return(tbl_exposure)
}

# Logit ponderado (post-estratificación u otros pesos de frecuencia) — survey::svyglm ----
fit_logit_model_weighted <- function(
    dependent, predictor, tiempo, contaminante, tipo,
    model_type, data,
    weight_var = "w_poststrat",
    conf.level = 0.95,
    adjustment = "Adjusted") {
  if (model_type == "single") {
    predictors_list <- predictor
  } else {
    predictors_list <- trimws(stringr::str_split(predictor, " \\+ ")[[1]])
  }

  missing_predictors <- predictors_list[!predictors_list %in% names(data)]
  if (length(missing_predictors) > 0 || !weight_var %in% names(data)) {
    return(data.frame(
      term = predictor,
      estimate = NA_real_,
      conf.low = NA_real_,
      conf.high = NA_real_,
      log_or = NA_real_,
      log_or_conf.low = NA_real_,
      log_or_conf.high = NA_real_,
      std.error = NA_real_,
      statistic = NA_real_,
      p.value = NA_real_,
      dependent_var = dependent,
      predictor = predictor,
      tiempo = tiempo,
      contaminante = contaminante,
      tipo = tipo,
      model_type = model_type,
      adjustment = adjustment,
      n = 0L,
      sum_w = NA_real_
    ))
  }

  data_subset <- data |>
    dplyr::filter(!is.na(.data[[dependent]])) |>
    dplyr::filter(!is.na(.data[[weight_var]]), .data[[weight_var]] > 0)

  for (pred in predictors_list) {
    data_subset <- data_subset |>
      dplyr::filter(!is.na(.data[[pred]]))
  }

  if (nrow(data_subset) < 10 || sum(data_subset[[weight_var]], na.rm = TRUE) < 1) {
    return(data.frame(
      term = predictor,
      estimate = NA_real_,
      conf.low = NA_real_,
      conf.high = NA_real_,
      log_or = NA_real_,
      log_or_conf.low = NA_real_,
      log_or_conf.high = NA_real_,
      std.error = NA_real_,
      statistic = NA_real_,
      p.value = NA_real_,
      dependent_var = dependent,
      predictor = predictor,
      tiempo = tiempo,
      contaminante = contaminante,
      tipo = tipo,
      model_type = model_type,
      adjustment = adjustment,
      n = nrow(data_subset),
      sum_w = sum(data_subset[[weight_var]], na.rm = TRUE)
    ))
  }

  if (identical(adjustment, "Adjusted")) {
    available_controls <- control_vars[control_vars %in% names(data_subset)]
    rhs <- if (length(available_controls) > 0) {
      paste(
        paste(predictors_list, collapse = " + "),
        paste("+", paste(available_controls, collapse = " + "))
      )
    } else {
      paste(predictors_list, collapse = " + ")
    }
  } else {
    rhs <- paste(predictors_list, collapse = " + ")
  }

  fml <- stats::as.formula(paste0(dependent, " ~ ", rhs))
  wform <- stats::as.formula(paste0("~", weight_var))
  des <- survey::svydesign(ids = ~1, weights = wform, data = data_subset)

  model_fit <- tryCatch(
    survey::svyglm(fml, design = des, family = stats::quasibinomial()),
    error = function(e) NULL
  )

  if (is.null(model_fit)) {
    return(data.frame(
      term = predictor,
      estimate = NA_real_,
      conf.low = NA_real_,
      conf.high = NA_real_,
      log_or = NA_real_,
      log_or_conf.low = NA_real_,
      log_or_conf.high = NA_real_,
      std.error = NA_real_,
      statistic = NA_real_,
      p.value = NA_real_,
      dependent_var = dependent,
      predictor = predictor,
      tiempo = tiempo,
      contaminante = contaminante,
      tipo = tipo,
      model_type = model_type,
      adjustment = adjustment,
      n = nrow(data_subset),
      sum_w = sum(data_subset[[weight_var]], na.rm = TRUE)
    ))
  }

  tbl <- broom::tidy(model_fit, conf.int = TRUE, conf.level = conf.level)
  z <- stats::qnorm(1 - (1 - conf.level) / 2)
  tbl_exposure <- tbl[tbl$term %in% predictors_list, ]
  sum_w <- sum(data_subset[[weight_var]], na.rm = TRUE)

  if (nrow(tbl_exposure) > 0) {
    tbl_exposure <- tbl_exposure |>
      dplyr::mutate(
        or = exp(.data$estimate),
        or_conf.low = exp(.data$estimate - z * .data$std.error),
        or_conf.high = exp(.data$estimate + z * .data$std.error),
        log_or = .data$estimate,
        log_or_conf.low = .data$estimate - z * .data$std.error,
        log_or_conf.high = .data$estimate + z * .data$std.error,
        estimate = .data$or,
        conf.low = .data$or_conf.low,
        conf.high = .data$or_conf.high,
        dependent_var = dependent,
        predictor = predictor,
        tiempo = tiempo,
        contaminante = contaminante,
        tipo = tipo,
        model_type = model_type,
        adjustment = adjustment,
        n = nrow(data_subset),
        sum_w = sum_w
      ) |>
      dplyr::select(
        term, estimate, conf.low, conf.high, log_or, log_or_conf.low, log_or_conf.high,
        std.error, statistic, p.value, dependent_var, predictor,
        tiempo, contaminante, tipo, model_type, adjustment, n, sum_w
      )
  } else {
    tbl_exposure <- data.frame(
      term = predictors_list[1],
      estimate = NA_real_,
      conf.low = NA_real_,
      conf.high = NA_real_,
      log_or = NA_real_,
      log_or_conf.low = NA_real_,
      log_or_conf.high = NA_real_,
      std.error = NA_real_,
      statistic = NA_real_,
      p.value = NA_real_,
      dependent_var = dependent,
      predictor = predictor,
      tiempo = tiempo,
      contaminante = contaminante,
      tipo = tipo,
      model_type = model_type,
      adjustment = adjustment,
      n = nrow(data_subset),
      sum_w = sum_w
    )
  }

  rm(model_fit, des)
  gc()
  return(tbl_exposure)
}

# Cox ponderado — survey::svycoxph ----
fit_cox_model_weighted <- function(
    dependent, predictor, tiempo, contaminante, tipo,
    model_type, data,
    time_var = "edad_gest",
    time_start = NULL,
    weight_var = "w_poststrat",
    conf.level = 0.95,
    adjustment = "Adjusted") {
  if (model_type == "single") {
    predictors_list <- predictor
  } else {
    predictors_list <- trimws(stringr::str_split(predictor, " \\+ ")[[1]])
  }

  missing_predictors <- predictors_list[!predictors_list %in% names(data)]
  if (length(missing_predictors) > 0 || !weight_var %in% names(data)) {
    return(data.frame(
      term = predictor,
      estimate = NA_real_,
      conf.low = NA_real_,
      conf.high = NA_real_,
      log_hr = NA_real_,
      log_hr_conf.low = NA_real_,
      log_hr_conf.high = NA_real_,
      std.error = NA_real_,
      statistic = NA_real_,
      p.value = NA_real_,
      dependent_var = dependent,
      predictor = predictor,
      tiempo = tiempo,
      contaminante = contaminante,
      tipo = tipo,
      model_type = model_type,
      adjustment = adjustment,
      n = 0L,
      sum_w = NA_real_
    ))
  }

  data_subset <- data |>
    dplyr::filter(!is.na(.data[[dependent]]), !is.na(.data[[time_var]])) |>
    dplyr::filter(!is.na(.data[[weight_var]]), .data[[weight_var]] > 0)

  use_delayed <- !is.null(time_start) &&
    is.character(time_start) && nzchar(time_start) && time_start %in% names(data_subset)
  if (use_delayed) {
    data_subset <- data_subset |>
      dplyr::filter(!is.na(.data[[time_start]])) |>
      dplyr::filter(.data[[time_start]] < .data[[time_var]])
  }

  for (pred in predictors_list) {
    data_subset <- data_subset |>
      dplyr::filter(!is.na(.data[[pred]]))
  }

  if (nrow(data_subset) < 10 || sum(data_subset[[weight_var]], na.rm = TRUE) < 1) {
    return(data.frame(
      term = predictor,
      estimate = NA_real_,
      conf.low = NA_real_,
      conf.high = NA_real_,
      log_hr = NA_real_,
      log_hr_conf.low = NA_real_,
      log_hr_conf.high = NA_real_,
      std.error = NA_real_,
      statistic = NA_real_,
      p.value = NA_real_,
      dependent_var = dependent,
      predictor = predictor,
      tiempo = tiempo,
      contaminante = contaminante,
      tipo = tipo,
      model_type = model_type,
      adjustment = adjustment,
      n = nrow(data_subset),
      sum_w = sum(data_subset[[weight_var]], na.rm = TRUE)
    ))
  }

  if (identical(adjustment, "Adjusted")) {
    available_controls <- control_vars[control_vars %in% names(data_subset)]
    rhs <- if (length(available_controls) > 0) {
      paste(
        paste(predictors_list, collapse = " + "),
        paste("+", paste(available_controls, collapse = " + "))
      )
    } else {
      paste(predictors_list, collapse = " + ")
    }
  } else {
    rhs <- paste(predictors_list, collapse = " + ")
  }

  surv_lhs <- if (use_delayed) {
    paste0("Surv(", time_start, ", ", time_var, ", ", dependent, ")")
  } else {
    paste0("Surv(", time_var, ", ", dependent, ")")
  }
  fml <- stats::as.formula(paste0(surv_lhs, " ~ ", rhs))
  wform <- stats::as.formula(paste0("~", weight_var))
  des <- survey::svydesign(ids = ~1, weights = wform, data = data_subset)

  model_fit <- tryCatch(
    survey::svycoxph(fml, design = des),
    error = function(e) NULL
  )

  if (is.null(model_fit)) {
    return(data.frame(
      term = predictor,
      estimate = NA_real_,
      conf.low = NA_real_,
      conf.high = NA_real_,
      log_hr = NA_real_,
      log_hr_conf.low = NA_real_,
      log_hr_conf.high = NA_real_,
      std.error = NA_real_,
      statistic = NA_real_,
      p.value = NA_real_,
      dependent_var = dependent,
      predictor = predictor,
      tiempo = tiempo,
      contaminante = contaminante,
      tipo = tipo,
      model_type = model_type,
      adjustment = adjustment,
      n = nrow(data_subset),
      sum_w = sum(data_subset[[weight_var]], na.rm = TRUE)
    ))
  }

  tbl <- broom::tidy(model_fit, exponentiate = TRUE, conf.int = TRUE, conf.level = conf.level)
  tbl_exposure <- tbl[tbl$term %in% predictors_list, ]
  sum_w <- sum(data_subset[[weight_var]], na.rm = TRUE)

  if (nrow(tbl_exposure) > 0) {
    tbl_exposure <- tbl_exposure |>
      dplyr::mutate(
        log_hr = log(.data$estimate),
        log_hr_conf.low = log(.data$conf.low),
        log_hr_conf.high = log(.data$conf.high),
        dependent_var = dependent,
        predictor = predictor,
        tiempo = tiempo,
        contaminante = contaminante,
        tipo = tipo,
        model_type = model_type,
        adjustment = adjustment,
        n = nrow(data_subset),
        sum_w = sum_w
      ) |>
      dplyr::select(
        term, estimate, conf.low, conf.high, log_hr, log_hr_conf.low, log_hr_conf.high,
        std.error, statistic, p.value, dependent_var, predictor,
        tiempo, contaminante, tipo, model_type, adjustment, n, sum_w
      )
  } else {
    tbl_exposure <- data.frame(
      term = predictors_list[1],
      estimate = NA_real_,
      conf.low = NA_real_,
      conf.high = NA_real_,
      log_hr = NA_real_,
      log_hr_conf.low = NA_real_,
      log_hr_conf.high = NA_real_,
      std.error = NA_real_,
      statistic = NA_real_,
      p.value = NA_real_,
      dependent_var = dependent,
      predictor = predictor,
      tiempo = tiempo,
      contaminante = contaminante,
      tipo = tipo,
      model_type = model_type,
      adjustment = adjustment,
      n = nrow(data_subset),
      sum_w = sum_w
    )
  }

  rm(model_fit, des)
  gc()
  return(tbl_exposure)
}

