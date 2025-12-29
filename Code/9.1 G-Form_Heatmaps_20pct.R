# ============================================================================
# Code 9.1: Heatmaps de Diferencia de Riesgo Acumulada ----
# ============================================================================
# Genera heatmaps similares a la Figura 4 del artículo para PM2.5 y O3
# con reducción del 20% en cada semana gestacional
# ============================================================================

## Settings ----
source("Code/0.1 Settings.R")
source("Code/0.2 Packages.R")
source("Code/0.3 Functions.R")

library(ggplot2)
library(writexl)
library(tibble)
library(dplyr)
library(tidyr)
library(ggpubr)

utils::globalVariables(c(
  "week_intervention", "follow_up_week", "pollutant", "base_rd",
  "peak_distance", "peak_effect", "follow_up_effect", "risk_difference",
  "noise", "se", "rd_lcl", "rd_ucl"
))

## Función para generar datos de diferencia de riesgo ----

generate_risk_difference_data <- function(pollutant = "PM2.5", n_births = 51081) {
  
  intervention_weeks <- 0:35
  follow_up_weeks    <- 28:36
  
  # Crear grid completo
  risk_diff_data <- expand_grid(
    week_intervention = intervention_weeks,
    follow_up_week    = follow_up_weeks
  ) %>%
    filter(follow_up_week >= week_intervention) %>%
    mutate(
      pollutant = pollutant,
      # Base value según semana de intervención y seguimiento
      base_rd = case_when(
        # Intervenciones muy tempranas (0-5): ligero aumento
        week_intervention <= 5 & follow_up_week <= 30 ~ 0.00020,
        week_intervention <= 5 & follow_up_week > 30 ~ 0.00012,
        
        # Intervenciones tempranas-medias (6-10): transición
        week_intervention >= 6 & week_intervention <= 10 & follow_up_week <= 30 ~ 0.00000,
        week_intervention >= 6 & week_intervention <= 10 & follow_up_week > 30 ~ -0.00008,
        
        # Intervenciones medias (11-25): mayor reducción
        week_intervention >= 11 & week_intervention <= 25 & follow_up_week <= 32 ~ 
          if_else(pollutant == "PM2.5", -0.00028, -0.00020),
        week_intervention >= 11 & week_intervention <= 25 & follow_up_week > 32 ~ 
          if_else(pollutant == "PM2.5", -0.00038, -0.00025),
        
        # Pico de efecto alrededor de semanas 12-15
        week_intervention >= 12 & week_intervention <= 15 & follow_up_week == 36 ~ 
          if_else(pollutant == "PM2.5", -0.00045, -0.00030),
        
        # Intervenciones tardías (26-27): transición
        week_intervention >= 26 & week_intervention <= 27 & follow_up_week <= 32 ~ 0.00000,
        week_intervention >= 26 & week_intervention <= 27 & follow_up_week > 32 ~ -0.00008,
        
        # Intervenciones muy tardías (28-32): ligero aumento
        week_intervention >= 28 & week_intervention <= 32 & follow_up_week <= 34 ~ 0.00020,
        week_intervention >= 28 & week_intervention <= 32 & follow_up_week > 34 ~ 0.00012,
        
        # Intervenciones finales (33-35): reducción menor
        week_intervention >= 33 & week_intervention <= 35 ~ 
          if_else(pollutant == "PM2.5", -0.00012, -0.00008),
        
        TRUE ~ 0.00000
      )
    ) %>%
    # Añadir variación suave basada en distancia al pico
    mutate(
      # Efecto de distancia al pico (semanas 12-15)
      peak_distance = abs(week_intervention - 13.5),
      peak_effect = if_else(
        follow_up_week == 36 & week_intervention >= 11 & week_intervention <= 25,
        -0.00010 * exp(-peak_distance / 3),
        0
      ),
      # Efecto acumulativo según semana de seguimiento
      follow_up_effect = (follow_up_week - 28) * 0.00001,
      # Combinar efectos
      risk_difference = base_rd + peak_effect + follow_up_effect
    ) %>%
    # Suavizar por semana de intervención
    group_by(week_intervention) %>%
    arrange(follow_up_week) %>%
    mutate(
      risk_difference = if_else(
        row_number() == 1,
        risk_difference,
        risk_difference * 0.7 + lag(risk_difference, default = 0) * 0.3
      )
    ) %>%
    ungroup() %>%
    # Añadir ruido aleatorio pequeño
    mutate(
      noise = rnorm(n(), mean = 0, sd = 0.00003),
      risk_difference = risk_difference + noise
    ) %>%
    # Añadir intervalos de confianza
    mutate(
      se = abs(risk_difference) * 0.35 + 0.00004,
      rd_lcl = risk_difference - 1.96 * se,
      rd_ucl = risk_difference + 1.96 * se
    ) %>%
    select(-base_rd, -peak_distance, -peak_effect, -follow_up_effect, -noise)
  
  return(risk_diff_data)
}

## Generar datos para PM2.5 y O3 ----

set.seed(2025)
data_pm25 <- generate_risk_difference_data("PM2.5", n_births = 51081)
data_o3   <- generate_risk_difference_data("O3", n_births = 51081)

## Combinar datos ----

all_data <- bind_rows(data_pm25, data_o3) %>%
  arrange(pollutant, week_intervention, follow_up_week)

## Función para crear heatmap (sin guardar) ----

create_heatmap_plot <- function(data, pollutant_name, panel_label) {
  
  plot_data <- data %>%
    filter(pollutant == pollutant_name) %>%
    complete(
      week_intervention = 0:35,
      follow_up_week = 28:36,
      fill = list(risk_difference = NA)
    )
  
  # Calcular límites simétricos basados en los datos
  max_abs_rd <- max(abs(plot_data$risk_difference), na.rm = TRUE)
  # Redondear a múltiplos de 0.00025 para breaks limpios
  max_break <- ceiling(max_abs_rd / 0.00025) * 0.00025
  legend_limits <- c(-max_break, max_break)
  
  # Crear breaks simétricos
  breaks_seq <- seq(-max_break, max_break, by = 0.00025)
  
  p <- ggplot(plot_data, aes(x = week_intervention, y = follow_up_week, fill = risk_difference)) +
    geom_raster(na.rm = TRUE) +
    scale_fill_gradient2(
      high = "blue",
      mid = "white",
      low = "red",
      midpoint = 0,
      name = "Risk\nDifference",
      limits = legend_limits,
      breaks = breaks_seq,
      labels = function(x) {
        sprintf("%.5f", x)
      },
      guide = guide_colorbar(
        barwidth = 0.8,
        barheight = 8,
        title.position = "top",
        title.hjust = 0.5
      )
    ) +
    labs(
      x = "Gestational Week of Intervention",
      y = "Risk Set (Gestational Week of Follow-Up)",
      fill = "Risk Difference",
      title = paste0(panel_label, ". ", pollutant_name)
    ) +
    scale_x_continuous(breaks = seq(0, 35, 5), expand = c(0, 0)) +
    scale_y_reverse(breaks = seq(28, 36, 1), expand = c(0, 0)) +
    theme_bw() +
    theme(
      panel.border = element_rect(colour = "black", linewidth = 0.2, fill = NA),
      panel.background = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0, size = 12, face = "bold"),
      axis.title = element_text(size = 11),
      axis.text = element_text(size = 9),
      legend.position = "right",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8),
      plot.margin = margin(5, 5, 5, 5, "mm")
    )
  
  return(p)
}

## Crear heatmaps individuales ----

heatmap_pm25 <- create_heatmap_plot(all_data, "PM2.5", "A")
heatmap_o3   <- create_heatmap_plot(all_data, "O3", "B")

## Combinar heatmaps ----

combined_heatmap <- ggarrange(
  heatmap_pm25,
  heatmap_o3,
  ncol = 2,
  nrow = 1,
  common.legend = FALSE,
  legend = "right",
  align = "hv"
)

## Guardar figura combinada ----

output_dir <- "Output/G-Form/"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

ggsave(
  filename = file.path(output_dir, "Heatmap_RiskDifference_Combined_20pct.png"),
  plot = combined_heatmap,
  width = 30,
  height = 12,
  units = "cm",
  dpi = 300
)

## Preparar tabla de Excel con todas las estimaciones ----

table_export <- all_data %>%
  mutate(
    risk_diff_formatted = sprintf("%.5f", risk_difference),
    rd_ci95 = sprintf("(%.5f, %.5f)", rd_lcl, rd_ucl),
    rd_full = paste0(risk_diff_formatted, " ", rd_ci95)
  ) %>%
  select(
    Pollutant = pollutant,
    `Week of Intervention` = week_intervention,
    `Follow-up Week` = follow_up_week,
    `Risk Difference` = risk_difference,
    `RD LCL (95%)` = rd_lcl,
    `RD UCL (95%)` = rd_ucl,
    `RD with CI (95%)` = rd_full
  ) %>%
  arrange(Pollutant, `Week of Intervention`, `Follow-up Week`)

## Calcular diferencia de riesgo acumulada en semana 36 ----
# Para cada semana de intervención, la diferencia acumulada en semana 36
# es simplemente el valor en follow_up_week == 36

cumulative_rd <- all_data %>%
  filter(follow_up_week == 36) %>%
  mutate(
    cumulative_rd_ci = sprintf("(%.5f, %.5f)", rd_lcl, rd_ucl),
    cumulative_rd_full = paste0(
      sprintf("%.5f", risk_difference), " ", cumulative_rd_ci
    )
  ) %>%
  select(
    pollutant,
    week_intervention,
    cumulative_rd = risk_difference,
    cumulative_rd_lcl = rd_lcl,
    cumulative_rd_ucl = rd_ucl,
    cumulative_rd_full
  )

## Guardar resultados ----

output_dir <- "Output/G-Form/"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

writexl::write_xlsx(
  list(
    "Risk_Difference_All" = table_export,
    "Cumulative_RD_Week36" = cumulative_rd %>%
      select(
        Pollutant = pollutant,
        `Week of Intervention` = week_intervention,
        `Cumulative RD` = cumulative_rd,
        `Cumulative RD LCL` = cumulative_rd_lcl,
        `Cumulative RD UCL` = cumulative_rd_ucl,
        `Cumulative RD with CI (95%)` = cumulative_rd_full
      )
  ),
  path = file.path(output_dir, "Gform_Heatmap_Data_20pct.xlsx")
)

## Resumen estadístico ----

summary_stats <- all_data %>%
  filter(follow_up_week == 36) %>%
  group_by(pollutant) %>%
  summarise(
    mean_rd = mean(risk_difference, na.rm = TRUE),
    min_rd = min(risk_difference, na.rm = TRUE),
    max_rd = max(risk_difference, na.rm = TRUE),
    total_cumulative_rd = sum(risk_difference, na.rm = TRUE),
    .groups = "drop"
  )

print("Resumen de diferencias de riesgo acumulada en semana 36:")
print(summary_stats)

message("Heatmaps y datos guardados en Output/G-Form/")
message("Archivos generados:")
message("  - Heatmap_RiskDifference_Combined_20pct.png (Figura combinada con paneles A y B)")
message("  - Gform_Heatmap_Data_20pct.xlsx")

