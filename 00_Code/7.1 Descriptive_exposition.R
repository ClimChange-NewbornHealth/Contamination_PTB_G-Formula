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

# Season based on astronomical change dates (Southern Hemisphere):
# 21/03 Summer -> Fall, 21/06 Fall -> Winter, 21/09 Winter -> Spring, 21/12 Spring -> Summer
get_season_from_date <- function(x_date) {
  md <- lubridate::month(x_date) * 100 + lubridate::day(x_date)
  dplyr::case_when(
    md >= 1221 | md < 321 ~ "Summer",
    md >= 321 & md < 621 ~ "Fall",
    md >= 621 & md < 921 ~ "Winter",
    TRUE ~ "Spring"
  )
}

exposure <- exposure |>
  mutate(
    date = as.Date(date),
    season = get_season_from_date(date)
  )

glimpse(exposure)

## Global descriptive pollulants ----

summarise_exposure <- function(df, pm_col, o3_col, no2_col, method_label) {
  df |>
    group_by(com, name_com, lat, long, sup) |>
    summarise(
      `PM2.5_Mean` = mean(.data[[pm_col]], na.rm = TRUE),
      `PM2.5_Min` = min(.data[[pm_col]], na.rm = TRUE),
      `PM2.5_Max` = max(.data[[pm_col]], na.rm = TRUE),
      `O3_Mean` = mean(.data[[o3_col]], na.rm = TRUE),
      `O3_Min` = min(.data[[o3_col]], na.rm = TRUE),
      `O3_Max` = max(.data[[o3_col]], na.rm = TRUE),
      `NO2_Mean` = mean(.data[[no2_col]], na.rm = TRUE),
      `NO2_Min` = min(.data[[no2_col]], na.rm = TRUE),
      `NO2_Max` = max(.data[[no2_col]], na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(
      Method = method_label,
      `Zip code` = com,
      Zip = name_com,
      Lat = lat,
      Lon = long,
      `Km²` = sup,
      .before = 1
    ) |>
    dplyr::select(
      `Zip code`,
      Zip,
      Lat,
      Lon,
      `Km²`,
      Method,
      `PM2.5_Mean`,
      `PM2.5_Min`,
      `PM2.5_Max`,
      `O3_Mean`,
      `O3_Min`,
      `O3_Max`,
      `NO2_Mean`,
      `NO2_Min`,
      `NO2_Max`
    ) |>
    mutate(
      across(
        -c(`Zip code`, Zip, Method),
        ~ formatC(round(.x, 2), format = "f", digits = 2, decimal.mark = ".")
      )
    ) |>
    arrange(`Zip code`)
}

tab_krg <- summarise_exposure(
  exposure,
  "pm25_ok_pred",
  "o3_ok_pred",
  "no2_ok_pred",
  "Kriging (ordinary)"
)

tab_idw <- summarise_exposure(
  exposure,
  "pm25_idw_pred",
  "o3_idw_pred",
  "no2_idw_pred",
  "IDW"
)

## Paper tab 
tab_publication <- tab_krg |>
  bind_rows(tab_idw) |> 
  arrange(`Zip code`)

## Save results
out_xlsx <- paste0(data_out, "Table_exposure_commune_PM25_O3_summary.xlsx")

writexl::write_xlsx(
  list(
    `Kriging_OK` = tab_krg,
    IDW = tab_idw,
    `Kriging_IDW_label` = tab_publication
  ),
  path = out_xlsx
)

## Histograms with intervention ----

build_hist_panel <- function(df, method_tag = c("KRG", "IDW")) {
  method_tag <- match.arg(method_tag)

  if (method_tag == "KRG") {
    pm_col <- "pm25_ok_pred"
    o3_col <- "o3_ok_pred"
    no2_col <- "no2_ok_pred"
  } else {
    pm_col <- "pm25_idw_pred"
    o3_col <- "o3_idw_pred"
    no2_col <- "no2_idw_pred"
  }

  cont_long <- df |>
    transmute(
      date = as.Date(date),
      season = season,
      PM2.5 = .data[[pm_col]],
      O3 = .data[[o3_col]],
      NO2 = .data[[no2_col]]
    ) |>
    pivot_longer(cols = c("PM2.5", "O3", "NO2"), names_to = "pollutant", values_to = "value")

  thresholds <- tibble::tribble(
    ~pollutant, ~x, ~label,
    "PM2.5", 15, "WHO: 15",
    "PM2.5", 50, "Chile: 50 µg/m³",
    "O3", 51, "WHO: 51",
    "O3", 61, "Chile: 61 ppbv",
    "NO2", 13, "WHO: 13",
    "NO2", 53, "Chile: 53 ppbv"
  )

  pollutant_cols <- list(
    "PM2.5" = c(fill = "#F4A261", color = "#C46D1A"),
    "O3"    = c(fill = "#2A9D8F", color = "#1C6E64"),
    "NO2"   = c(fill = "#8E7CC3", color = "#5D4A99")
  )

  make_single_hist <- function(pollutant_name, period_name, panel_letter) {
    dat <- cont_long |>
      filter(
        pollutant == pollutant_name,
        if (period_name == "Overall") TRUE else season == period_name
      )

    thr <- thresholds |>
      filter(pollutant == pollutant_name)

    x_range <- range(dat$value, na.rm = TRUE)
    x_span <- diff(x_range)
    if (!is.finite(x_span) || x_span <= 0) x_span <- 1
    x_pad <- if (pollutant_name == "O3") 0.28 * x_span else 0.18 * x_span
    x_text <- thr$x + 0.02 * x_span
    max_x <- max(c(dat$value, thr$x), na.rm = TRUE)

    x_label <- if (pollutant_name == "PM2.5") {
      expression("Concentration (" * mu * "g/" * m^3 * ")")
    } else {
      "Concentration (ppbv)"
    }

    pollutant_md <- switch(
      pollutant_name,
      "PM2.5" = "PM<sub>2.5</sub>",
      "O3" = "O<sub>3</sub>",
      "NO2" = "NO<sub>2</sub>"
    )
    title_md <- paste0(panel_letter, ". ", pollutant_md, " ", period_name)

    ggplot(dat, aes(x = value)) +
      geom_histogram(
        fill = pollutant_cols[[pollutant_name]]["fill"],
        color = pollutant_cols[[pollutant_name]]["fill"],
        alpha = 0.55,
        binwidth = 0.5
      ) +
      geom_vline(data = thr, aes(xintercept = x), linewidth = 0.5, linetype = "longdash", color = "black") +
      geom_text(
        data = thr,
        aes(x = x_text, y = Inf, label = label),
        hjust = 0,
        vjust = 2.1,
        size = 3
      ) +
      labs(
        title = title_md,
        x = x_label,
        y = "Frequency"
      ) +
      scale_x_continuous(
        limits = c(min(x_range, na.rm = TRUE), max_x + x_pad),
        labels = scales::label_number(decimal.mark = ".", big.mark = "")
      ) +
      theme_light() +
      coord_cartesian(clip = "off") +
      theme(
        panel.grid = element_blank(),
        legend.position = "none",
        plot.title = ggtext::element_markdown(size = 10),
        plot.margin = margin(t = 6, r = 8, b = 6, l = 6)
      )
  }

  pA <- make_single_hist("PM2.5", "Overall", "A")
  pB <- make_single_hist("PM2.5", "Winter", "B")
  pC <- make_single_hist("NO2", "Overall", "C")
  pD <- make_single_hist("NO2", "Winter", "D")
  pE <- make_single_hist("O3", "Overall", "E")
  pF <- make_single_hist("O3", "Summer", "F")

  compiled <- ggpubr::ggarrange(
    pA, pB,
    pC, pD,
    pE, pF,
    ncol = 2,
    nrow = 3,
    align = "hv"
  )

  ggsave(
    filename = paste0(data_out, "Histogram_", method_tag, "_panel_compiled.png"),
    plot = compiled,
    res = 300,
    width = 32,
    height = 25,
    units = "cm",
    bg = "white",
    scale = 0.9,
    device = ragg::agg_png
  )

  ggsave(
    filename = paste0(data_out, "Histogram_", method_tag, "_A_PM25_Overall.png"),
    plot = pA, res = 300, width = 16, height = 12, units = "cm", bg = "white", device = ragg::agg_png
  )
  ggsave(
    filename = paste0(data_out, "Histogram_", method_tag, "_B_PM25_Winter.png"),
    plot = pB, res = 300, width = 16, height = 12, units = "cm", bg = "white", device = ragg::agg_png
  )
  ggsave(
    filename = paste0(data_out, "Histogram_", method_tag, "_C_O3_Overall.png"),
    plot = pE, res = 300, width = 16, height = 12, units = "cm", bg = "white", device = ragg::agg_png
  )
  ggsave(
    filename = paste0(data_out, "Histogram_", method_tag, "_D_O3_Summer.png"),
    plot = pF, res = 300, width = 16, height = 12, units = "cm", bg = "white", device = ragg::agg_png
  )
  ggsave(
    filename = paste0(data_out, "Histogram_", method_tag, "_E_NO2_Overall.png"),
    plot = pC, res = 300, width = 16, height = 12, units = "cm", bg = "white", device = ragg::agg_png
  )
  ggsave(
    filename = paste0(data_out, "Histogram_", method_tag, "_F_NO2_Winter.png"),
    plot = pD, res = 300, width = 16, height = 12, units = "cm", bg = "white", device = ragg::agg_png
  )

  invisible(list(compiled = compiled, A = pA, B = pB, C = pC, D = pD, E = pE, F = pF))
}

plots_krg <- build_hist_panel(exposure, "KRG")
plots_idw <- build_hist_panel(exposure, "IDW")


## Histograms with intervention by municipality ----

build_hist_facet_municipality <- function(df, pollutant = c("PM2.5", "NO2", "O3"), method = c("KRG", "IDW")) {
  pollutant <- match.arg(pollutant)
  method <- match.arg(method)

  col_map <- list(
    KRG = c("PM2.5" = "pm25_ok_pred", "NO2" = "no2_ok_pred", "O3" = "o3_ok_pred"),
    IDW = c("PM2.5" = "pm25_idw_pred", "NO2" = "no2_idw_pred", "O3" = "o3_idw_pred")
  )

  season_target <- ifelse(pollutant == "O3", "Summer", "Winter")
  pollutant_col <- col_map[[method]][[pollutant]]

  comuna_levels <- df |>
    distinct(com, name_com) |>
    arrange(com) |>
    pull(name_com)

  cont_data_aux <- df |>
    transmute(
      com = com,
      name_com = factor(name_com, levels = comuna_levels),
      season = season,
      value = .data[[pollutant_col]]
    ) |>
    filter(!is.na(value))

  dat_overall <- cont_data_aux |>
    mutate(season_plot = "Overall")

  dat_season <- cont_data_aux |>
    filter(season == season_target) |>
    mutate(season_plot = season_target)

  plot_data <- bind_rows(dat_overall, dat_season) |>
    mutate(season_plot = factor(season_plot, levels = c("Overall", season_target)))

  season_col <- switch(
    pollutant,
    "PM2.5" = "#F4A261",
    "NO2" = "#8E7CC3",
    "O3" = "#2A9D8F"
  )
  cols <- c("Overall" = "#BDBDBD", stats::setNames(season_col, season_target))

  thr <- switch(
    pollutant,
    "PM2.5" = tibble::tribble(
      ~x, ~type,
      50, "Chile guideline (50 µg/m³)",
      15, "WHO guideline (15 µg/m³)"
    ),
    "NO2" = tibble::tribble(
      ~x, ~type,
      53, "Chile guideline (53 ppbv)",
      13, "WHO guideline (13 ppbv)"
    ),
    "O3" = tibble::tribble(
      ~x, ~type,
      61, "Chile guideline (61 ppbv)",
      51, "WHO guideline (51 ppbv)"
    )
  )

  x_lab <- if (pollutant == "PM2.5") {
    expression("Concentration (" * mu * "g/" * m^3 * ")")
  } else {
    "Concentration (ppbv)"
  }

  p <- ggplot(plot_data, aes(x = value, fill = season_plot, color = season_plot)) +
    geom_histogram(position = "identity", alpha = 0.35, binwidth = 0.5) +
    scale_fill_manual(values = cols, name = NULL, breaks = c("Overall", season_target)) +
    scale_color_manual(values = cols, name = NULL, breaks = c("Overall", season_target)) +
    labs(
      x = x_lab,
      y = "Frequency"
    ) +
    geom_vline(data = thr, aes(xintercept = x, linetype = type), linewidth = 0.5, color = "black") +
    scale_linetype_manual(
      values = setNames(
        c("longdash", "dotdash"),
        c(thr$type[[1]], thr$type[[2]])
      ),
      breaks = c(thr$type[[2]], thr$type[[1]]),
      name = NULL
    ) +
    guides(
      fill = guide_legend(order = 1),
      color = guide_legend(order = 1),
      linetype = guide_legend(order = 2)
    ) +
    facet_wrap(~name_com, scales = "free", ncol = 5) +
    scale_x_continuous(labels = scales::label_number(decimal.mark = ".", big.mark = "")) +
    theme_light() +
    theme(
      strip.background = element_rect(fill = "white", color = "black"),
      strip.text = element_text(color = "black"),
      panel.grid = element_blank(),
      legend.position = "top",
      legend.title = element_text(),
      legend.box = "horizontal",
      legend.spacing.x = unit(0.3, "cm"),
      legend.spacing.y = unit(0, "cm"),
      legend.margin = margin(t = 0, r = 0, b = 0, l = 0)
    )

  out_file <- paste0(
    data_out,
    "Histogram_FACET_",
    gsub("\\.", "", pollutant),
    "_",
    method,
    ".png"
  )

  ggsave(
    filename = out_file,
    plot = p,
    res = 300,
    width = 20,
    height = 25,
    units = "cm",
    scaling = 0.7,
    bg = "white",
    device = ragg::agg_png
  )

  invisible(p)
}

facet_pm_krg <- build_hist_facet_municipality(exposure, pollutant = "PM2.5", method = "KRG")
facet_pm_idw <- build_hist_facet_municipality(exposure, pollutant = "PM2.5", method = "IDW")
facet_no2_krg <- build_hist_facet_municipality(exposure, pollutant = "NO2", method = "KRG")
facet_no2_idw <- build_hist_facet_municipality(exposure, pollutant = "NO2", method = "IDW")
facet_o3_krg <- build_hist_facet_municipality(exposure, pollutant = "O3", method = "KRG")
facet_o3_idw <- build_hist_facet_municipality(exposure, pollutant = "O3", method = "IDW")

## Annual-season summary table (Summer/Winter): mean, min and max by pollutant and estimator ----

annual_summary_long <- exposure |>
  mutate(
    year = lubridate::year(date),
    season = factor(season, levels = c("Summer", "Winter"))
  ) |>
  filter(season %in% c("Summer", "Winter")) |>
  group_by(year, season) |>
  summarise(
    pm25_krg_mean = mean(pm25_ok_pred, na.rm = TRUE),
    pm25_krg_min = min(pm25_ok_pred, na.rm = TRUE),
    pm25_krg_max = max(pm25_ok_pred, na.rm = TRUE),
    pm25_idw_mean = mean(pm25_idw_pred, na.rm = TRUE),
    pm25_idw_min = min(pm25_idw_pred, na.rm = TRUE),
    pm25_idw_max = max(pm25_idw_pred, na.rm = TRUE),
    no2_krg_mean = mean(no2_ok_pred, na.rm = TRUE),
    no2_krg_min = min(no2_ok_pred, na.rm = TRUE),
    no2_krg_max = max(no2_ok_pred, na.rm = TRUE),
    no2_idw_mean = mean(no2_idw_pred, na.rm = TRUE),
    no2_idw_min = min(no2_idw_pred, na.rm = TRUE),
    no2_idw_max = max(no2_idw_pred, na.rm = TRUE),
    o3_krg_mean = mean(o3_ok_pred, na.rm = TRUE),
    o3_krg_min = min(o3_ok_pred, na.rm = TRUE),
    o3_krg_max = max(o3_ok_pred, na.rm = TRUE),
    o3_idw_mean = mean(o3_idw_pred, na.rm = TRUE),
    o3_idw_min = min(o3_idw_pred, na.rm = TRUE),
    o3_idw_max = max(o3_idw_pred, na.rm = TRUE),
    .groups = "drop"
  ) |>
  pivot_longer(
    cols = -c(year, season),
    names_to = c("pollutant", "estimator", "stat"),
    names_pattern = "(pm25|no2|o3)_(krg|idw)_(mean|min|max)",
    values_to = "value"
  ) |>
  mutate(
    pollutant = recode(pollutant, pm25 = "PM2.5", no2 = "NO2", o3 = "O3"),
    estimator = recode(estimator, krg = "Kriging", idw = "IDW"),
    stat = recode(stat, mean = "Mean", min = "Min", max = "Max")
  ) |>
  pivot_wider(names_from = stat, values_from = value) |>
  arrange(year, season, pollutant, estimator) |>
  mutate(
    across(
      c(Mean, Min, Max),
      ~ formatC(round(.x, 2), format = "f", digits = 2, decimal.mark = ".")
    )
  )

writexl::write_xlsx(
  list(annual_summary_long = annual_summary_long),
  path = paste0(data_out, "Table_Annual_Summary_Contaminants_Estimators.xlsx")
)

## Annual-season summary table (Summer/Winter) by municipality ----

mun_summary_long <- exposure |>
  mutate(
    season = factor(season, levels = c("Summer", "Winter"))
  ) |>
  filter(season %in% c("Summer", "Winter")) |>
  group_by(com, name_com, season) |>
  summarise(
    pm25_krg_mean = mean(pm25_ok_pred, na.rm = TRUE),
    pm25_krg_min = min(pm25_ok_pred, na.rm = TRUE),
    pm25_krg_max = max(pm25_ok_pred, na.rm = TRUE),
    pm25_idw_mean = mean(pm25_idw_pred, na.rm = TRUE),
    pm25_idw_min = min(pm25_idw_pred, na.rm = TRUE),
    pm25_idw_max = max(pm25_idw_pred, na.rm = TRUE),
    no2_krg_mean = mean(no2_ok_pred, na.rm = TRUE),
    no2_krg_min = min(no2_ok_pred, na.rm = TRUE),
    no2_krg_max = max(no2_ok_pred, na.rm = TRUE),
    no2_idw_mean = mean(no2_idw_pred, na.rm = TRUE),
    no2_idw_min = min(no2_idw_pred, na.rm = TRUE),
    no2_idw_max = max(no2_idw_pred, na.rm = TRUE),
    o3_krg_mean = mean(o3_ok_pred, na.rm = TRUE),
    o3_krg_min = min(o3_ok_pred, na.rm = TRUE),
    o3_krg_max = max(o3_ok_pred, na.rm = TRUE),
    o3_idw_mean = mean(o3_idw_pred, na.rm = TRUE),
    o3_idw_min = min(o3_idw_pred, na.rm = TRUE),
    o3_idw_max = max(o3_idw_pred, na.rm = TRUE),
    .groups = "drop"
  ) |>
  pivot_longer(
    cols = -c(com, name_com, season),
    names_to = c("pollutant", "estimator", "stat"),
    names_pattern = "(pm25|no2|o3)_(krg|idw)_(mean|min|max)",
    values_to = "value"
  ) |>
  mutate(
    pollutant = recode(pollutant, pm25 = "PM2.5", no2 = "NO2", o3 = "O3"),
    estimator = recode(estimator, krg = "Kriging", idw = "IDW"),
    stat = recode(stat, mean = "Mean", min = "Min", max = "Max")
  ) |>
  pivot_wider(names_from = stat, values_from = value) |>
  arrange(com, name_com, season, pollutant, desc(estimator)) |>
  mutate(
    across(
      c(Mean, Min, Max),
      ~ formatC(round(.x, 2), format = "f", digits = 2, decimal.mark = ".")
    )
  )

writexl::write_xlsx(
  list(mun_summary_long = mun_summary_long),
  path = paste0(data_out, "Table_Municipality_Summary_Contaminants_Estimators.xlsx")
)


## Time distribution plots (daily mean across municipalities) ----

cont_data_mean <- exposure |>
  group_by(date) |>
  summarise(
    pm25_krg = mean(pm25_ok_pred, na.rm = TRUE),
    pm25_idw = mean(pm25_idw_pred, na.rm = TRUE),
    no2_krg = mean(no2_ok_pred, na.rm = TRUE),
    no2_idw = mean(no2_idw_pred, na.rm = TRUE),
    o3_krg = mean(o3_ok_pred, na.rm = TRUE),
    o3_idw = mean(o3_idw_pred, na.rm = TRUE),
    .groups = "drop"
  )

lab_ugm3 <- expression("Concentration (" * mu * "g/" * m^3 * ")")
lab_ppb <- "Concentration (ppbv)"

gvars <- list(
  list(name = "pm25_krg", var = "pm25_krg", title_expr = expression("A. PM"[2.5] * " - Kriging")),
  list(name = "pm25_idw", var = "pm25_idw", title_expr = expression("B. PM"[2.5] * " - IDW")),
  list(name = "no2_krg", var = "no2_krg", title_expr = expression("C. NO"[2] * " - Kriging")),
  list(name = "no2_idw", var = "no2_idw", title_expr = expression("D. NO"[2] * " - IDW")),
  list(name = "o3_krg", var = "o3_krg", title_expr = expression("E. O"[3] * " - Kriging")),
  list(name = "o3_idw", var = "o3_idw", title_expr = expression("F. O"[3] * " - IDW"))
)

plots_time <- list()

for (i in seq_along(gvars)) {
  v <- gvars[[i]]
  ylab <- if (startsWith(v$var, "o3") || startsWith(v$var, "no2")) lab_ppb else lab_ugm3

  p <- ggplot(cont_data_mean, aes(x = date, y = .data[[v$var]])) +
    geom_point(size = 0.5, alpha = 0.1) +
    geom_smooth(method = "loess", span = 0.05, se = TRUE, linewidth = 0.6, color = "#2F6DF6") +
    labs(title = v$title_expr, x = NULL, y = ylab) +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
    theme_light() +
    theme(
      plot.title = element_text(),
      panel.grid = element_blank(),
      strip.background = element_rect(fill = "white"),
      strip.text = element_text(size = 11, color = "black", hjust = 0),
      axis.text.y = element_text(size = 9),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8)
    )

  plots_time[[i]] <- p
  names(plots_time)[i] <- v$name
}

fig_time_final <- ggpubr::ggarrange(
  plots_time[[1]],
  plots_time[[2]],
  plots_time[[3]],
  plots_time[[4]],
  plots_time[[5]],
  plots_time[[6]],
  ncol = 2,
  nrow = 3,
  align = "hv"
)

ggsave(
  filename = paste0(data_out, "Time_distribution_pm25_no2_o3.png"),
  plot = fig_time_final,
  res = 300,
  width = 20,
  height = 22,
  units = "cm",
  scaling = 0.9,
  bg = "white",
  device = ragg::agg_png
)


