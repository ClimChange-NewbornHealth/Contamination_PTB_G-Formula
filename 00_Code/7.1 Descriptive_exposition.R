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

ptb_mun <- rio::import(paste0(data_inp, "Data_births_ptb_mun.xlsx")) |>
  transmute(
    `Zip code` = com,
    `% PTB` = formatC(round(ptb * 100, 2), format = "f", digits = 2, decimal.mark = ".")
  )

add_ptb_column <- function(tab) {
  tab |>
    left_join(ptb_mun, by = "Zip code") |>
    relocate(`% PTB`, .after = Zip)
}

tab_krg <- add_ptb_column(tab_krg)
tab_idw <- add_ptb_column(tab_idw)

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


## Map with exposures ----

com_codes_rm <- chilemapas::codigos_territoriales |>
  dplyr::filter(codigo_region == 13) |>
  dplyr::mutate(codigo_comuna = as.numeric(codigo_comuna))

com_suburb <- c(
  unique(com_codes_rm$codigo_comuna[com_codes_rm$nombre_provincia == "Santiago"]),
  13201L
)

exposure_urb <- exposure |>
  dplyr::filter(com %in% com_suburb)

map_pollutant_palettes <- list(
  "PM2.5" = c("#00E400", "#FFFF00", "#FF7E00", "#FF0000", "#8F3F97", "#7E0023"),
  "NO2" = c("#7D8C96", "#9E8B7E", "#8D6E63", "#6D4C41", "#4E342E", "#1B120F"),
  "O3" = c("#FFFFE5", "#FFF9C4", "#FFEB3B", "#FFD54F", "#FFB300", "#F57C00")
)

map_pollutant_specs <- list(
  list(
    label = "PM2.5",
    krg_col = "pm25_ok_pred",
    idw_col = "pm25_idw_pred",
    legend_title = expression("PM"[2.5] * " mean (" * mu * "g/m"^3 * ")"),
    xlab = expression("PM"[2.5] * " daily concentration (" * mu * "g/" * m^3 * ")")
  ),
  list(
    label = "NO2",
    krg_col = "no2_ok_pred",
    idw_col = "no2_idw_pred",
    legend_title = expression("NO"[2] * " mean (ppbv)"),
    xlab = expression("NO"[2] * " daily concentration (ppbv)")
  ),
  list(
    label = "O3",
    krg_col = "o3_ok_pred",
    idw_col = "o3_idw_pred",
    legend_title = expression("O"[3] * " mean (ppbv)"),
    xlab = expression("O"[3] * " daily concentration (ppbv)")
  )
)

exposure_reduction_factor <- 0.80

map_fill_alpha <- 0.58
exposure_target_quantile <- 0.80
warning_icon_path <- file.path(data_out, "assets", "warning_icon.png")

warning_icon_is_valid <- function(path) {
  if (!file.exists(path)) {
    return(FALSE)
  }
  tryCatch({
    png::readPNG(path)
    TRUE
  }, error = function(e) FALSE)
}

ensure_warning_icon_png <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (warning_icon_is_valid(path)) {
    return(invisible(path))
  }
  grDevices::png(path, width = 128, height = 128, bg = "transparent")
  grid::grid.newpage()
  grid::grid.polygon(
    x = grid::unit(c(0.12, 0.88, 0.5), "npc"),
    y = grid::unit(c(0.08, 0.08, 0.92), "npc"),
    gp = grid::gpar(fill = "#FFCC00", col = "#1A1A1A", lwd = 3)
  )
  grid::grid.text(
    "!",
    x = grid::unit(0.5, "npc"),
    y = grid::unit(0.46, "npc"),
    gp = grid::gpar(fontsize = 42, fontface = "bold", col = "#1A1A1A")
  )
  grDevices::dev.off()
  invisible(path)
}

add_warning_icon_layer <- function(
    p,
    warn_pts,
    icon_path,
    width = 0.042,
    height = 0.026) {
  if (is.null(warn_pts) || !nrow(warn_pts)) {
    return(p)
  }
  coords <- sf::st_coordinates(warn_pts)
  img <- png::readPNG(icon_path)
  half_lon <- width / 2
  half_lat <- height / 2
  for (i in seq_len(nrow(coords))) {
    p <- p + ggplot2::annotation_raster(
      img,
      xmin = coords[i, 1] - half_lon,
      xmax = coords[i, 1] + half_lon,
      ymin = coords[i, 2] - half_lat,
      ymax = coords[i, 2] + half_lat
    )
  }
  p
}

density_below_cut_label_coords <- function(values, p_cut, x_min) {
  dens <- stats::density(values)
  idx <- which(dens$x <= p_cut)
  if (!length(idx)) {
    return(list(x = (x_min + p_cut) / 2, y = max(dens$y, na.rm = TRUE) * 0.35))
  }
  w <- dens$y[idx]
  x_center <- stats::weighted.mean(dens$x[idx], w)
  y_center <- stats::weighted.mean(dens$y[idx], w) / 2
  list(x = x_center, y = y_center)
}

build_map_sf_urb <- function(mun_means) {
  geo <- chilemapas::mapa_comunas |>
    dplyr::mutate(codigo_comuna = as.numeric(codigo_comuna)) |>
    dplyr::filter(codigo_comuna %in% com_suburb) |>
    dplyr::left_join(mun_means, by = c("codigo_comuna" = "com"))

  if (!inherits(geo, "sf")) {
    geo <- sf::st_as_sf(geo)
  }
  sf::st_transform(geo, 4326)
}

fetch_rm_basemap <- function(lims, zoom = 11) {
  if (!requireNamespace("maptiles", quietly = TRUE)) {
    warning("maptiles not available; maps will render without basemap.")
    return(NULL)
  }
  bounds <- sf::st_bbox(
    c(lims$xmin, lims$ymin, lims$xmax, lims$ymax),
    crs = sf::st_crs(4326)
  )
  providers <- c(
    "Esri.WorldGrayCanvas",
    "OpenStreetMap",
    "CartoDB.PositronNoLabels"
  )
  for (prov in providers) {
    out <- try(
      maptiles::get_tiles(
        bounds,
        provider = prov,
        zoom = zoom,
        crop = TRUE,
        retina = FALSE
      ),
      silent = TRUE
    )
    if (!inherits(out, "try-error")) {
      message("Basemap OK (", prov, ", zoom ", zoom, ").")
      return(out)
    }
  }
  warning("Basemap download failed for all providers.")
  NULL
}

build_map_context_urb <- function(map_sf) {
  bb <- sf::st_bbox(map_sf)
  pad <- 0.02
  lims <- list(
    xmin = bb$xmin - pad,
    xmax = bb$xmax + pad,
    ymin = bb$ymin - pad,
    ymax = bb$ymax + pad
  )

  comunas_lim <- chilemapas::mapa_comunas |>
    dplyr::mutate(codigo_comuna = as.numeric(codigo_comuna)) |>
    dplyr::filter(codigo_comuna %in% com_suburb) |>
    sf::st_as_sf() |>
    sf::st_transform(4326)

  urb_boundary <- sf::st_union(comunas_lim) |>
    sf::st_boundary()

  list(
    lims = lims,
    map_base = fetch_rm_basemap(lims),
    comunas_lim = comunas_lim,
    urb_boundary = urb_boundary
  )
}

pollutant_value_limits <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) {
    return(c(0, 1))
  }
  r_min <- min(x)
  r_max <- max(x)
  if (!is.finite(r_max) || r_max <= r_min) {
    r_max <- r_min + 1
  }
  c(r_min, r_max)
}

legend_breaks_five <- function(x) {
  rng <- range(x, na.rm = TRUE)
  seq(rng[[1L]], rng[[2L]], length.out = 5)
}

legend_labels_1dec <- function(breaks) {
  vapply(breaks, function(b) {
    format(round(b, 1), nsmall = 1, decimal.mark = ".", trim = TRUE)
  }, character(1))
}

guide_fill_map_exposure <- function() {
  ggplot2::guide_colorbar(
    barwidth = 12,
    barheight = 0.35,
    nbin = 5,
    label.position = "bottom",
    title.position = "top",
    direction = "horizontal"
  )
}

panel_tag_theme <- function(show_legend = TRUE) {
  theme_light(base_size = 9) +
    theme(
      legend.position = if (show_legend) "top" else "none",
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 7),
      legend.margin = margin(t = 0, b = 4),
      plot.margin = margin(t = 8, r = 6, b = 0, l = 6),
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      plot.title = element_blank(),
      plot.tag = element_text(size = 11, face = "bold", hjust = 0, vjust = 1),
      plot.tag.position = c(0.01, 0.99)
    )
}

theme_density_panel <- function() {
  panel_tag_theme(show_legend = FALSE) +
    theme(
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 9),
      plot.margin = margin(t = 28, r = 6, b = 4, l = 6)
    )
}

commune_warning_points <- function(map_sf, value_col, daily_values) {
  vals <- map_sf[[value_col]]
  thr <- stats::quantile(
    daily_values,
    probs = exposure_target_quantile,
    na.rm = TRUE,
    type = 7
  )
  flagged <- map_sf[!is.na(vals) & vals > thr, , drop = FALSE]
  if (!nrow(flagged)) {
    return(list(points = NULL, threshold = thr))
  }
  pts <- sf::st_point_on_surface(flagged)
  list(points = pts, threshold = thr)
}

plot_commune_mean_map <- function(
    map_sf,
    value_col,
    palette_colors,
    legend_title,
    map_ctx,
    daily_values,
    panel_tag = NULL) {
  dat <- sf::st_make_valid(map_sf)
  dat$map_value <- dat[[value_col]]
  lims <- map_ctx$lims
  fill_limits <- pollutant_value_limits(dat$map_value)
  brks <- legend_breaks_five(dat$map_value)
  warn <- commune_warning_points(dat, "map_value", daily_values)

  layer_basemap <- function() {
    if (!is.null(map_ctx$map_base) && requireNamespace("tidyterra", quietly = TRUE)) {
      tidyterra::geom_spatraster_rgb(data = map_ctx$map_base, maxcell = 5e5)
    } else {
      ggplot2::annotate(
        "rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "grey92"
      )
    }
  }

  p <- ggplot2::ggplot() +
    layer_basemap() +
    ggplot2::geom_sf(
      data = map_ctx$comunas_lim,
      fill = NA,
      color = "gray40",
      linewidth = 0.35,
      inherit.aes = FALSE
    ) +
    ggplot2::geom_sf(
      data = map_ctx$urb_boundary,
      color = "gray15",
      linewidth = 0.9,
      inherit.aes = FALSE
    ) +
    ggplot2::geom_sf(
      data = dat,
      ggplot2::aes(fill = map_value),
      color = "white",
      linewidth = 0.5,
      alpha = map_fill_alpha
    )

  p <- add_warning_icon_layer(
    p,
    warn_pts = warn$points,
    icon_path = warning_icon_path,
    width = 0.042,
    height = 0.026
  )

  p <- p +
    ggplot2::scale_fill_gradientn(
      colors = palette_colors,
      limits = fill_limits,
      oob = scales::squish,
      name = legend_title,
      na.value = "gray90",
      breaks = brks,
      labels = legend_labels_1dec,
      guide = guide_fill_map_exposure()
    ) +
    ggplot2::coord_sf(
      crs = sf::st_crs(4326),
      expand = FALSE,
      xlim = c(lims$xmin, lims$xmax),
      ylim = c(lims$ymin, lims$ymax)
    ) +
    ggplot2::labs(tag = panel_tag) +
    panel_tag_theme(show_legend = TRUE)

  if (requireNamespace("ggspatial", quietly = TRUE)) {
    p <- p +
      ggspatial::annotation_scale(location = "bl", width_hint = 0.2) +
      ggspatial::annotation_north_arrow(
        location = "tr",
        height = grid::unit(0.8, "cm"),
        width = grid::unit(0.6, "cm"),
        style = ggspatial::north_arrow_fancy_orienteering()
      )
  }

  p
}

plot_daily_density_panel <- function(
    daily_values,
    palette_colors,
    fill_limits,
    xlab,
    panel_tag = NULL) {
  dat_natural <- tibble::tibble(value = daily_values)
  dat_natural <- dat_natural[is.finite(dat_natural$value), , drop = FALSE]
  dat_reduced <- tibble::tibble(value = dat_natural$value * exposure_reduction_factor)

  r_min <- fill_limits[1]
  r_max <- fill_limits[2]
  p80 <- stats::quantile(
    dat_natural$value,
    probs = exposure_target_quantile,
    na.rm = TRUE,
    type = 7
  )
  line_col <- palette_colors[[length(palette_colors)]]
  label_pos <- density_below_cut_label_coords(dat_natural$value, p80, r_min)

  ggplot2::ggplot() +
    ggplot2::geom_density(
      data = dat_natural,
      ggplot2::aes(x = value),
      fill = "#BDBDBD",
      color = "#616161",
      alpha = 0.55,
      linewidth = 0.65
    ) +
    ggplot2::geom_density(
      data = dat_reduced,
      ggplot2::aes(x = value, fill = ggplot2::after_stat(x)),
      color = line_col,
      alpha = 0.62,
      linewidth = 0.75
    ) +
    ggplot2::scale_fill_gradientn(
      colors = palette_colors,
      limits = c(r_min * exposure_reduction_factor, r_max * exposure_reduction_factor),
      guide = "none"
    ) +
    ggplot2::scale_x_continuous(
      limits = c(r_min, r_max),
      labels = scales::label_number(decimal.mark = ".", big.mark = "")
    ) +
    ggplot2::annotate(
      "text",
      x = label_pos$x,
      y = label_pos$y,
      label = "\u2190 20%",
      size = 4.2,
      fontface = "bold",
      hjust = 0.5,
      vjust = 0.5
    ) +
    ggplot2::labs(
      tag = panel_tag,
      x = xlab,
      y = "Density"
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    theme_density_panel()
}

build_map_exposure_panel <- function(
    df,
    method_tag = c("KRG", "IDW"),
    map_ctx = NULL) {
  method_tag <- match.arg(method_tag)

  mun_means <- df |>
    dplyr::group_by(com, name_com) |>
    dplyr::summarise(
      pm25_ok_pred = mean(pm25_ok_pred, na.rm = TRUE),
      no2_ok_pred = mean(no2_ok_pred, na.rm = TRUE),
      o3_ok_pred = mean(o3_ok_pred, na.rm = TRUE),
      pm25_idw_pred = mean(pm25_idw_pred, na.rm = TRUE),
      no2_idw_pred = mean(no2_idw_pred, na.rm = TRUE),
      o3_idw_pred = mean(o3_idw_pred, na.rm = TRUE),
      .groups = "drop"
    )

  map_sf <- build_map_sf_urb(mun_means)
  if (is.null(map_ctx)) {
    map_ctx <- build_map_context_urb(map_sf)
  }
  panel_tags <- LETTERS[1:6]

  map_plots <- list()
  density_plots <- list()

  for (i in seq_along(map_pollutant_specs)) {
    spec <- map_pollutant_specs[[i]]
    val_col <- if (method_tag == "KRG") spec$krg_col else spec$idw_col
    pal <- map_pollutant_palettes[[spec$label]]
    daily_limits_poll <- pollutant_value_limits(df[[val_col]])

    map_plots[[spec$label]] <- plot_commune_mean_map(
      map_sf = map_sf,
      value_col = val_col,
      palette_colors = pal,
      legend_title = spec$legend_title,
      map_ctx = map_ctx,
      daily_values = df[[val_col]],
      panel_tag = panel_tags[i]
    )

    density_plots[[spec$label]] <- plot_daily_density_panel(
      daily_values = df[[val_col]],
      palette_colors = pal,
      fill_limits = daily_limits_poll,
      xlab = spec$xlab,
      panel_tag = panel_tags[i + 3L]
    )
  }

  row_maps <- patchwork::wrap_plots(
    map_plots[["PM2.5"]],
    map_plots[["NO2"]],
    map_plots[["O3"]],
    ncol = 3,
    nrow = 1
  )
  row_density <- patchwork::wrap_plots(
    density_plots[["PM2.5"]],
    density_plots[["NO2"]],
    density_plots[["O3"]],
    ncol = 3,
    nrow = 1
  )
  panel <- row_maps / row_density +
    patchwork::plot_layout(heights = c(1.12, 0.88), guides = "keep")

  out_file <- paste0(
    data_out,
    "Map_Exposure_daily_mean_and_density_",
    method_tag,
    ".png"
  )

  ggplot2::ggsave(
    filename = out_file,
    plot = panel,
    width = 32,
    height = 22,
    units = "cm",
    res = 300,
    bg = "white",
    device = ragg::agg_png
  )

  message("Map exposure panel (", method_tag, ") saved: ", out_file)
  invisible(list(panel = panel, map_plots = map_plots, density_plots = density_plots))
}

ensure_warning_icon_png(warning_icon_path)

mun_means_map <- exposure_urb |>
  dplyr::group_by(com, name_com) |>
  dplyr::summarise(
    pm25_ok_pred = mean(pm25_ok_pred, na.rm = TRUE),
    .groups = "drop"
  )
map_sf_urb_panel <- build_map_sf_urb(mun_means_map)
map_ctx_urb <- build_map_context_urb(map_sf_urb_panel)

map_panel_krg <- build_map_exposure_panel(exposure_urb, "KRG", map_ctx = map_ctx_urb)
map_panel_idw <- build_map_exposure_panel(exposure_urb, "IDW", map_ctx = map_ctx_urb)

