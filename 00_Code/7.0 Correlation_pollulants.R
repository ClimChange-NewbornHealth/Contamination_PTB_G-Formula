# Code 7: IDW vs Kriging (density + correlation panels) ----

rm(list = ls())

## Settings ----
source("00_Code/0.1 Settings.R")
source("00_Code/0.2 Packages.R")

data_inp <- "01_Data/Output/"
data_out <- "02_Output/idw_vs_kriging/"
dir.create(data_out, recursive = TRUE, showWarnings = FALSE)

## Load exposure and prepare contaminants ----
exposure <- rio::import(paste0(data_inp, "Contamination_Climate_Data_2010_2020.RData"))

cont_data <- exposure |>
  dplyr::transmute(
    pm25_krg = pm25_ok_pred,
    pm25_idw = pm25_idw_pred,
    no2_krg = no2_ok_pred,
    no2_idw = no2_idw_pred,
    o3_krg = o3_ok_pred,
    o3_idw = o3_idw_pred
  )

## Correlation between contaminants ----

data_pm <- cont_data |>
  dplyr::select(pm25_krg, pm25_idw) |>
  tidyr::pivot_longer(cols = dplyr::everything(), names_to = "method", values_to = "value") |>
  dplyr::mutate(
    method = dplyr::recode(method, "pm25_krg" = "Kriging", "pm25_idw" = "IDW"),
    contaminant = "PM2.5"
  )

data_no2 <- cont_data |>
  dplyr::select(no2_krg, no2_idw) |>
  tidyr::pivot_longer(cols = dplyr::everything(), names_to = "method", values_to = "value") |>
  dplyr::mutate(
    method = dplyr::recode(method, "no2_krg" = "Kriging", "no2_idw" = "IDW"),
    contaminant = "NO2"
  )

data_o3 <- cont_data |>
  dplyr::select(o3_krg, o3_idw) |>
  tidyr::pivot_longer(cols = dplyr::everything(), names_to = "method", values_to = "value") |>
  dplyr::mutate(
    method = dplyr::recode(method, "o3_krg" = "Kriging", "o3_idw" = "IDW"),
    contaminant = "O3"
  )

data_dens <- dplyr::bind_rows(data_pm, data_no2, data_o3)

plot_density_panel <- function(data, contaminant_name, xlim_range = NULL, x_label = "") {
  data_sub <- dplyr::filter(data, .data$contaminant == contaminant_name)

  ggplot(data_sub, aes(x = value, fill = method, color = method)) +
    geom_density(alpha = 0.3, linewidth = 0.8) +
    scale_fill_manual(values = c("Kriging" = "#1f77b4", "IDW" = "#ff7f0e")) +
    scale_color_manual(values = c("Kriging" = "#1f77b4", "IDW" = "#ff7f0e")) +
    scale_x_continuous(labels = scales::label_number(decimal.mark = ".")) +
    scale_y_continuous(labels = scales::label_number(decimal.mark = ".")) +
    labs(x = x_label, y = "Density") +
    theme_light() +
    theme(
      panel.grid = element_blank(),
      legend.position = "top",
      legend.title = element_blank(),
      plot.title = element_text(size = 11, hjust = 0)
    ) +
    coord_cartesian(xlim = xlim_range)
}

g1 <- plot_density_panel(
  data_dens,
  contaminant_name = "PM2.5",
  xlim_range = c(0, 100),
  x_label = expression("Concentration (" * mu * "g/" * m^3 * ")")
) +
  labs(title = expression("A. PM"[2.5] * " (" * mu * "g/" * m^3 * ")"))

g2 <- plot_density_panel(
  data_dens,
  contaminant_name = "NO2",
  xlim_range = c(0, 80),
  x_label = "Concentration (ppbv)"
) +
  labs(title = expression("B. NO"[2] * " (ppbv)"))

g3 <- plot_density_panel(
  data_dens,
  contaminant_name = "O3",
  xlim_range = c(0, 40),
  x_label = "Concentration (ppbv)"
) +
  labs(title = expression("C. O"[3] * " (ppbv)"))

plt1 <- ggpubr::ggarrange(g1, g2, g3, ncol = 3, common.legend = TRUE)

ggsave(
  filename = paste0(data_out, "Distribution_krg_idw_pm25_no2_o3.png"),
  plot = plt1,
  res = 300,
  width = 30,
  height = 10,
  units = "cm",
  scaling = 0.9,
  bg = "white",
  device = ragg::agg_png
)

### Scatter-plots ----

make_cor_label <- function(x, y, data) {
  test <- cor.test(data[[x]], data[[y]], use = "complete.obs")
  r <- formatC(test$estimate, format = "f", digits = 2, decimal.mark = ".")
  p <- sub(",", ".", format.pval(test$p.value, digits = 3, eps = .001), fixed = TRUE)
  paste0("r = ", r, ", p = ", p)
}

plot_corr_panel <- function(data, x, y, xlab, ylab) {
  label_txt <- make_cor_label(x, y, data)

  ggplot(data, aes_string(x = x, y = y)) +
    stat_bin2d(bins = 100) +
    scale_fill_gradient(low = "#fee8c8", high = "#e34a33", name = "Count") +
    geom_smooth(method = "lm", formula = y ~ x, color = "#08519c", alpha = 0.3, linewidth = 1) +
    annotate("text", x = Inf, y = Inf, label = label_txt, hjust = 1.1, vjust = 1.5, size = 4) +
    labs(x = xlab, y = ylab) +
    scale_x_continuous(labels = scales::label_number(decimal.mark = ".")) +
    scale_y_continuous(labels = scales::label_number(decimal.mark = ".")) +
    theme_light() +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(size = 11, hjust = 0)
    )
}

g4 <- plot_corr_panel(
  data = cont_data,
  x = "pm25_krg",
  y = "pm25_idw",
  xlab = "Kriging",
  ylab = "IDW"
) +
  scale_x_continuous(limits = c(0, 100), labels = scales::label_number(decimal.mark = ".")) +
  scale_y_continuous(limits = c(0, 100), labels = scales::label_number(decimal.mark = ".")) +
  labs(title = expression("D. PM"[2.5] * " (" * mu * "g/" * m^3 * ")")) +
  theme(legend.position = "right")

g5 <- plot_corr_panel(
  data = cont_data,
  x = "no2_krg",
  y = "no2_idw",
  xlab = "Kriging",
  ylab = "IDW"
) +
  scale_x_continuous(limits = c(0, 80), labels = scales::label_number(decimal.mark = ".")) +
  scale_y_continuous(limits = c(0, 80), labels = scales::label_number(decimal.mark = ".")) +
  labs(title = expression("E. NO"[2] * " (ppbv)")) +
  theme(legend.position = "right")

g6 <- plot_corr_panel(
  data = cont_data,
  x = "o3_krg",
  y = "o3_idw",
  xlab = "Kriging",
  ylab = "IDW"
) +
  scale_x_continuous(limits = c(0, 35), labels = scales::label_number(decimal.mark = ".")) +
  scale_y_continuous(limits = c(0, 35), labels = scales::label_number(decimal.mark = ".")) +
  labs(title = expression("F. O"[3] * " (ppbv)")) +
  theme(legend.position = "right")

plt2 <- ggpubr::ggarrange(g4, g5, g6, ncol = 3, common.legend = TRUE, legend = "none")

ggsave(
  filename = paste0(data_out, "Correlation_krg_idw_pm25_no2_o3.png"),
  plot = plt2,
  res = 300,
  width = 30,
  height = 10,
  units = "cm",
  scaling = 0.9,
  bg = "white",
  device = ragg::agg_png
)

### Complete plot -----

fig_all <- ggpubr::ggarrange(plt1, plt2, ncol = 1)

ggsave(
  filename = paste0(data_out, "Descriptive_krg_idw_pm25_no2_o3.png"),
  plot = fig_all,
  res = 300,
  width = 30,
  height = 20,
  units = "cm",
  scaling = 0.9,
  bg = "white",
  device = ragg::agg_png
)

### Full correlation panel (7x7) using metan::corr_plot ----

corr_vars_krg <- exposure |>
  dplyr::transmute(
    PM25 = pm25_ok_pred,
    NO2 = no2_ok_pred,
    O3 = o3_ok_pred,
    TMAX = tmax,
    TMIN = tmin,
    TAD = TAD,
    NDVI = ndvi
  )

corr_vars_idw <- exposure |>
  dplyr::transmute(
    PM25 = pm25_idw_pred,
    NO2 = no2_idw_pred,
    O3 = o3_idw_pred,
    TMAX = tmax,
    TMIN = tmin,
    TAD = TAD,
    NDVI = ndvi
  )

panel_krg <- metan::corr_plot(
  corr_vars_krg,
  prob = 0.01,
  shape.point = 21,
  col.point = "gray35",
  fill.point = "gray70",
  size.point = 0.25,
  alpha.point = 0.12,
  maxsize = 6,
  minsize = 6,
  smooth = TRUE,
  size.line = 1,
  col.smooth = "#2F5FB3",
  confint = TRUE,
  col.sign = "white",
  alpha.sign = 0,
  col.up.panel = "black",
  col.lw.panel = "black",
  col.dia.panel = "black",
  bins = 35,
  pan.spacing = 0.08,
  lab.position = "tl",
  decimal.mark = "."
) +
  theme(plot.margin = margin(12, 12, 12, 12))

panel_idw <- metan::corr_plot(
  corr_vars_idw,
  prob = 0.01,
  shape.point = 21,
  col.point = "gray35",
  fill.point = "gray70",
  size.point = 0.25,
  alpha.point = 0.12,
  maxsize = 6,
  minsize = 6,
  smooth = TRUE,
  size.line = 1,
  col.smooth = "#2F5FB3",
  confint = TRUE,
  col.sign = "white",
  alpha.sign = 0,
  col.up.panel = "black",
  col.lw.panel = "black",
  col.dia.panel = "black",
  bins = 35,
  pan.spacing = 0.08,
  lab.position = "tl",
  decimal.mark = "."
) +
  theme(plot.margin = margin(12, 12, 12, 12))

ggsave(
  filename = paste0(data_out, "Correlation_panel_7x7_Kriging.png"),
  plot = panel_krg,
  res = 300,
  width = 24,
  height = 24,
  units = "cm",
  bg = "white",
  device = ragg::agg_png
)

ggsave(
  filename = paste0(data_out, "Correlation_panel_7x7_IDW.png"),
  plot = panel_idw,
  res = 300,
  width = 24,
  height = 24,
  units = "cm",
  bg = "white",
  device = ragg::agg_png
)


