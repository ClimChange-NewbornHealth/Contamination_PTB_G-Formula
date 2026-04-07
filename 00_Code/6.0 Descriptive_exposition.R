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

