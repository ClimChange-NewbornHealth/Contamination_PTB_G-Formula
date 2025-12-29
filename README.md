# Contamination and Preterm Birth: G-Computation Analysis with Distributed Lag Models

## Descripción del Proyecto

Este proyecto evalúa el efecto de la contaminación atmosférica (PM2.5 y O3) sobre el parto pretérmino en Chile mediante métodos de g-computation (g-formula) con modelos de rezago distribuido (Distributed Lag Models, DLNM). El análisis utiliza datos de nacimientos y exposición a contaminantes atmosféricos para estimar escenarios contrafactuales de reducción de contaminación y su impacto potencial en la prevalencia de partos pretérminos.

## Objetivos

### Objetivo Principal
Evaluar escenarios contrafactuales de reducción de PM2.5 y O3 sobre el parto pretérmino mediante g-computation con modelos de rezago distribuido.

### Objetivos Específicos
1. Procesar y preparar datos de contaminación atmosférica (PM2.5 y O3) utilizando métodos de interpolación espacial (kriging e IDW)
2. Procesar datos de nacimientos y construir historiales de exposición semanal durante el embarazo
3. Estimar modelos de supervivencia (Cox) y modelos de rezago distribuido (DLNM) para evaluar asociaciones entre exposición y parto pretérmino
4. Implementar g-computation para estimar efectos contrafactuales de intervenciones de reducción de contaminación
5. Calcular métricas de impacto poblacional: prevalencia, razón de riesgos, diferencia de riesgos, casos atribuibles y riesgo atribuible

## Datos

### Fuentes de Datos

#### Datos de Contaminación Atmosférica
- **PM2.5 y O3**: Series temporales diarias de concentraciones de material particulado fino (PM2.5) y ozono (O3)
- **Métodos de interpolación espacial**:
  - Kriging (`pm25_krg`, `o3_krg`)
  - Inverse Distance Weighting - IDW (`pm25_idw`, `o3_idw`)
- **Período**: 2010-2020
- **Ubicación**: `Data/Input/Contant_series/`

#### Datos de Nacimientos
- **Fuente**: Registros de nacimientos en Chile
- **Período**: 2010-2020
- **Variables principales**:
  - Edad gestacional (semanas)
  - Parto pretérmino (definido como <37 semanas)
  - Variables sociodemográficas (sexo, edad y educación de padres, etc.)
  - Variables temporales (año, mes, COVID-19)
- **Ubicación**: `Data/Input/Nacimientos/`

#### Datos de Vulnerabilidad Social
- **Índice de Vulnerabilidad Social (SOVI)**: Clasificación de vulnerabilidad por comuna
- **Ubicación**: `Data/Input/SOVI/`

### Datos Procesados

Los datos procesados se almacenan en `Data/Output/`:
- `series_births_exposition_pm25_o3_kriging_idw_long.RData`: Base de datos principal con historiales de exposición semanal
- `series_contamination_pm25_o3_kriging_idw.RData`: Series temporales de contaminación
- `births_2010_2020.RData`: Base de datos de nacimientos procesada

## Métodos

### 1. Procesamiento de Datos

#### 1.1 Procesamiento de Datos de Contaminación (`Code/1.0 Process_data_contam.R`)
- Carga y procesamiento de datos de PM2.5 y O3
- Interpolación espacial mediante kriging e IDW
- Generación de series temporales diarias por comuna

#### 1.2 Procesamiento de Datos de Nacimientos (`Code/2.0 Births_process_data.R`, `Code/2.1 Births_process_data_long_weeksg.R`)
- Limpieza y preparación de datos de nacimientos
- Definición de variables de resultado (parto pretérmino, bajo peso al nacer, etc.)
- Construcción de variables temporales y sociodemográficas
- Creación de estructura de datos longitudinal por semana gestacional

#### 1.3 Construcción de Historiales de Exposición (`Code/3.0 Exposure_data.R`, `Code/4.0 Join_CONT_BW.R`, `Code/4.1 Join_CONT_BW_weeks.R`)
- Asignación de exposición semanal durante el embarazo
- Cálculo de métricas de exposición agregadas (promedios por trimestre, exposición acumulada)
- Normalización por IQR (Interquartile Range) para análisis

### 2. Análisis Descriptivos (`Code/5.0 Descriptives.R` y subsecuentes)
- Estadísticas descriptivas de variables de resultado y exposición
- Análisis de tendencias temporales
- Caracterización de asociaciones preliminares

### 3. Modelos de Supervivencia (`Code/6.0 Cox_models_full_sample.R` y subsecuentes)
- Modelos de regresión de Cox para evaluar asociaciones entre exposición y parto pretérmino
- Análisis por subgrupos (verano para O3, invierno para PM2.5)
- Cálculo de fracciones atribuibles poblacionales (PAF)

### 4. Modelos de Rezago Distribuido (DLNM) (`Code/7.0 DLM_Preterm_Full_Period.R`, `Code/8.0 DLNM_Preterm_Full_Period.R`)
- Modelos de rezago distribuido para capturar efectos de exposición a lo largo del tiempo
- Especificación de funciones base para exposición y rezago
- Estimación de efectos acumulados y por ventanas de tiempo específicas

### 5. G-Computation con DLNM (`Code/9.0 G-Form_DLNM_Preterm_Full_Period.R`)

#### 5.1 Especificación del Modelo
- Modelos logísticos para probabilidad de evento por semana gestacional (semanas 28-36)
- Incorporación de crossbasis de exposición mediante DLNM
- Ajuste por confusores: sexo, edad y educación de padres, variables temporales, vulnerabilidad social, temperatura

#### 5.2 Escenarios de Intervención
- **Escenario observado**: Curso natural sin intervención
- **Intervención 1**: PM2.5 < 20 µg/m³
- **Intervención 2**: PM2.5 < 15 µg/m³
- **Intervención 3**: PM2.5 < 10 µg/m³
- **Intervención 4**: PM2.5 < 5 µg/m³
- **Intervención 5**: Reducción del 20% en PM2.5 y O3 cada semana

#### 5.3 Implementación de G-Computation
1. **Ajuste del modelo**: Estimación de modelos logísticos por semana de riesgo
2. **Predicción contrafactual**: Aplicación de escenarios de intervención a historiales de exposición
3. **Cálculo de probabilidades**: Estimación de probabilidades de supervivencia bajo cada escenario
4. **Agregación**: Cálculo de métricas poblacionales (prevalencia, casos, razón de riesgos, diferencia de riesgos)

#### 5.4 Inferencia Estadística
- Intervalos de confianza del 95% mediante simulación paramétrica (bootstrap paramétrico)
- 200 iteraciones de bootstrap para estimación de incertidumbre

### 6. Visualización (`Code/9.1 G-Form_Heatmaps_20pct.R`)
- Heatmaps de diferencias de riesgo por escenario de intervención
- Visualización de efectos por contaminante y tipo de intervención

## Estructura del Código

```
Code/
├── 0.1 Settings.R              # Configuración general (locale, opciones)
├── 0.2 Packages.R              # Carga e instalación de paquetes
├── 0.3 Functions.R              # Funciones auxiliares personalizadas
│
├── 1.0 Process_data_contam.R   # Procesamiento datos contaminación
├── 2.0 Births_process_data.R   # Procesamiento datos nacimientos
├── 2.1 Births_process_data_long_weeksg.R  # Estructura longitudinal por semana gestacional
├── 3.0 Exposure_data.R         # Construcción historiales de exposición
├── 4.0 Join_CONT_BW.R          # Unión datos contaminación y nacimientos
├── 4.1 Join_CONT_BW_weeks.R   # Unión con estructura semanal
│
├── 5.0 Descriptives.R          # Análisis descriptivos principales
├── 5.1 Descritive_Births_trends_preterm.R
├── 5.2 Descriptives_exposition.R
├── 5.3 Descriptives_asociation_contaminant_data.R
├── 5.4 Descriptives_exposition_test.R
│
├── 6.0 Cox_models_full_sample.R      # Modelos Cox muestra completa
├── 6.1 Cox_models_ozone_summer.R     # Modelos Cox O3 verano
├── 6.2 Cox_models_pm_winter.R        # Modelos Cox PM2.5 invierno
├── 6.3 Cox_models_plots_kriging.R    # Visualización modelos kriging
├── 6.4 Cox_models_plots_idw.R        # Visualización modelos IDW
├── 6.5 IQR_cox_models.R              # Modelos normalizados por IQR
├── 6.6 PAF_cox_models.R              # Fracciones atribuibles poblacionales
│
├── 7.0 DLM_Preterm_Full_Period.R     # Modelos de rezago distribuido (DLM)
├── 8.0 DLNM_Preterm_Full_Period.R    # Modelos de rezago distribuido no lineal (DLNM)
│
├── 9.0 G-Form_DLNM_Preterm_Full_Period.R  # G-computation con DLNM
└── 9.1 G-Form_Heatmaps_20pct.R            # Visualización resultados g-computation
```

### Código Adicional
- `Code/OBS/`: Análisis con datos observados (no interpolados)
- `Code/paf_cox_manual.R`: Cálculo manual de fracciones atribuibles
- `Code/test.R`: Scripts de prueba

## Parámetros Principales

### Parámetros de G-Computation (`Code/9.0 G-Form_DLNM_Preterm_Full_Period.R`)
- `max_follow_up`: 37 semanas (semanas 0-36)
- `risk_weeks`: 28:36 (semanas de riesgo analizadas)
- `lag_df`: 4 grados de libertad para el rezago en DLNM
- `boot_iter`: 200 iteraciones para intervalos de confianza
- `baseline_scenario`: "observed" (escenario de referencia)

### Variables de Exposición
- `pm25_krg_week_iqr`: PM2.5 kriging normalizado por IQR
- `pm25_idw_week_iqr`: PM2.5 IDW normalizado por IQR
- `o3_krg_week_iqr`: O3 kriging normalizado por IQR
- `o3_idw_week_iqr`: O3 IDW normalizado por IQR

### Confusores Incluidos
- Variables sociodemográficas: sexo, edad y educación de padres, ocupación
- Variables temporales: mes y año de inicio de embarazo, período COVID-19
- Vulnerabilidad social: índice SOVI
- Variables ambientales: temperatura

## Resultados

*[Sección pendiente - resultados preliminares]*

Los resultados del análisis de g-computation se almacenan en:
- `Output/G-Form/Gform_DLNM_results.xlsx`: Métricas principales por escenario
- `Output/G-Form/Gform_DLNM_results.rds`: Objetos R completos con resultados
- `Output/G-Form/Intervention_pm25_o3.xlsx`: Definición de escenarios de intervención
- `Output/G-Form/Heatmap_*.png`: Visualizaciones de diferencias de riesgo

## Requisitos del Sistema

### Software
- R (versión 4.0 o superior)
- RStudio (recomendado)

### Paquetes R Principales
- `dplyr`, `tidyr`, `purrr`: Manipulación de datos
- `survival`, `flexsurv`: Modelos de supervivencia
- `dlnm`: Modelos de rezago distribuido
- `furrr`, `future`: Procesamiento paralelo
- `rio`, `writexl`: Importación/exportación de datos
- `ggplot2`, `patchwork`: Visualización

Ver `Code/0.2 Packages.R` para la lista completa de paquetes.

### Recursos Computacionales
- Procesamiento paralelo configurado para utilizar `detectCores() - 4` workers
- Memoria: Se recomienda al menos 16 GB RAM para análisis completos
- Tiempo estimado: El análisis completo de g-computation puede tomar varias horas dependiendo del tamaño de muestra

## Uso

### Ejecución del Análisis Completo

1. **Configuración inicial**:
   ```r
   source("Code/0.1 Settings.R")
   source("Code/0.2 Packages.R")
   source("Code/0.3 Functions.R")
   ```

2. **Procesamiento de datos** (ejecutar en orden):
   ```r
   # 1. Procesar datos de contaminación
   source("Code/1.0 Process_data_contam.R")
   
   # 2. Procesar datos de nacimientos
   source("Code/2.0 Births_process_data.R")
   source("Code/2.1 Births_process_data_long_weeksg.R")
   
   # 3. Construir historiales de exposición
   source("Code/3.0 Exposure_data.R")
   source("Code/4.0 Join_CONT_BW.R")
   source("Code/4.1 Join_CONT_BW_weeks.R")
   ```

3. **Análisis descriptivos**:
   ```r
   source("Code/5.0 Descriptives.R")
   ```

4. **Modelos de asociación**:
   ```r
   source("Code/6.0 Cox_models_full_sample.R")
   source("Code/8.0 DLNM_Preterm_Full_Period.R")
   ```

5. **G-computation**:
   ```r
   source("Code/9.0 G-Form_DLNM_Preterm_Full_Period.R")
   ```

### Análisis con Submuestras

Para pruebas rápidas, se puede modificar el parámetro `sample_n_ids` en `Code/9.0 G-Form_DLNM_Preterm_Full_Period.R`:
```r
sample_n_ids <- 2000  # Probar con 2000 individuos
```

## Referencias

*[Sección pendiente - referencias bibliográficas]*

### Métodos Clave
- **G-computation/G-formula**: Robins (1986), Hernán & Robins (2020)
- **Distributed Lag Models**: Gasparrini (2011), Gasparrini et al. (2017)
- **Causal Inference**: Hernán & Robins (2020)

## Contacto

*[Información de contacto pendiente]*

## Licencia

Ver archivo `LICENSE` para más detalles.

---

**Nota**: Este proyecto forma parte del proyecto ANID Iniciación 2024-2027 sobre Contaminación y Parto Pretérmino en Chile.
