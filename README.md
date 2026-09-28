# Estimating the Perinatal Health Benefits of Hypothetical Pollution Interventions in Santiago, Chile Using Parametric G-Computation :factory: :baby:

![GitHub Repo stars](https://img.shields.io/github/stars/ClimChange-NewbornHealth/Contamination_PTB_G-Formula)
![GitHub watchers](https://img.shields.io/github/watchers/ClimChange-NewbornHealth/Contamination_PTB_G-Formula)
![GitHub forks](https://img.shields.io/github/forks/ClimChange-NewbornHealth/Contamination_PTB_G-Formula)
![GitHub commit activity](https://img.shields.io/github/commit-activity/t/ClimChange-NewbornHealth/Contamination_PTB_G-Formula)
![GitHub contributors](https://img.shields.io/github/contributors/ClimChange-NewbornHealth/Contamination_PTB_G-Formula)
![GitHub last commit](https://img.shields.io/github/last-commit/ClimChange-NewbornHealth/Contamination_PTB_G-Formula)
![GitHub language count](https://img.shields.io/github/languages/count/ClimChange-NewbornHealth/Contamination_PTB_G-Formula)
![GitHub top language](https://img.shields.io/github/languages/top/ClimChange-NewbornHealth/Contamination_PTB_G-Formula)
![GitHub License](https://img.shields.io/github/license/ClimChange-NewbornHealth/Contamination_PTB_G-Formula)
![GitHub repo file or directory count](https://img.shields.io/github/directory-file-count/ClimChange-NewbornHealth/Contamination_PTB_G-Formula)
![GitHub code size in bytes](https://img.shields.io/github/languages/code-size/ClimChange-NewbornHealth/Contamination_PTB_G-Formula)

## :moneybag: Funding

**FONDECYT Nº 11240322**: Climate change and urban health: how air pollution, temperature, and city structure relate to preterm birth

Additional support: **(CR)²**, Chile, FONDAP/ANID 1523A0002

## :busts_in_silhouette: Research Team

:mailbox_with_mail: **Estela Blanco** (<estela.blanco@uc.cl>), **Principal Investigator / Corresponding Author**

:mailbox_with_mail: **José Daniel Conejeros** (<jdconejeros@uc.cl>), **Research Assistant / Repository Manager**

**Research Collaborators**: Ismael Bravo, Felipe Cornejo, Axel Osses, & Tarik Benmarhnia

## :pushpin: Publication

*Work in progress.*

---

## :dart: Project Overview

### Background

Ambient air pollution is a major environmental risk factor for adverse perinatal outcomes. Preterm birth (delivery before 37 completed weeks of gestation) remains a leading cause of neonatal mortality and long-term morbidity. While epidemiological evidence links PM₂.₅, NO₂, and O₃ to preterm birth, fewer studies have translated observed exposure–response associations into population-level estimates of the health benefits of realistic pollution-reduction scenarios, particularly in Latin American urban settings with heterogeneous exposure patterns.

### Objective

To estimate the population-level impact of hypothetical reductions in PM₂.₅, NO₂, and O₃ on the cumulative risk of preterm birth among singleton births in urban Santiago, Chile (2010–2020), using distributed-lag Cox models combined with parametric g-computation.

### Methods (as reported in the manuscript)

- **Design**: Retrospective population-based cohort (DEIS birth records, 2010–2020).
- **Study area**: Urban conurbation of Santiago (32 municipalities in the Province of Santiago plus Puente Alto; 33 municipalities total).
- **Analytic sample**: **713,918** singleton live births with gestational age ≥28 weeks at delivery, complete covariates, plausible birthweight-for-gestational age (Alexander et al., 1996), and no fixed-cohort bias (gestation start ≥ January 1, 2010; delivery on or before March 20, 2020). **51,081** preterm births (**7.2%**). Exclusions in Figure 1.
- **Exposure**: Daily PM₂.₅, NO₂, and O₃ from the national monitoring network; municipality-day concentrations via **ordinary kriging** to municipal administrative centers (primary). **IDW** in parallel for sensitivity. Values below detection limits set to the limit. Weekly series \(X_{iw}\) (weeks 1–44) and weekly mean temperature; full-pregnancy NDVI (MODIS MOD13Q1, Kalman-imputed daily values); municipal SOVI (low, medium-low, medium-high).
- **Outcome**: Preterm birth (<37 completed weeks).
- **Distributed-lag models (main text)**: For each gestational week \(w = 1,\ldots,44\), separate Cox models on completed gestational age with delayed entry at week 28, Efron ties, time-weighted lag through week \(w-1\) (equation 1 in the manuscript), and **multipollutant** adjustment (concurrent weekly co-pollutants plus weekly and lagged terms for the target pollutant). Covariates: newborn sex; parental age, education, occupation; conception month and year; COVID-19 indicator (deliveries after March 1, 2020); SOVI; weekly temperature; full-pregnancy NDVI. Single-pollutant DLMs and trimester/overall averages are in the supplement (Tables S5–S10, Figure S6).
- **G-computation (main text)**: Parametric g-formula using the **multipollutant** natural-course Cox models (risk weeks 28–36). Counterfactual histories apply a **20% proportional reduction** (\(\pi = 0.20\)) to all gestational weeks of the target pollutant, with lag terms recomputed; co-pollutants, temperature, and covariates stay at observed values. Discrete-time risks from Breslow baseline hazards; global effects at week 36 (prevalence, expected cases, RR, RD, AR, PAF). **Single-week** 20% interventions for critical-window heatmaps (Figure 5). **95% CIs**: 2.5th–97.5th percentiles of a **parametric bootstrap** (250 replicates) resampling Cox coefficients (supplementary methods S1, Table S1).
- **Supplementary g-computation**: 10%, 20%, and 30% proportional reductions under **single-pollutant** weekly Cox models (not reported in the main tables).

### Key findings (main text)

- **Exposure (2010–2020)**: Mean daily kriging concentrations about **25.9 µg/m³** (PM₂.₅), **21.3 ppbv** (NO₂), and **14.0 ppbv** (O₃); strong winter–summer contrasts (Figure 2; Table S2). Kriging and IDW agreed closely (e.g., \(r = 0.96\) for PM₂.₅; Figure S3).
- **DLM (Figure 3)**: Multipollutant week-specific HRs per 1-unit increase were near the null but differed by timing (e.g., O₃ elevated in early pregnancy and week 34; NO₂ inverse in early gestation and positive around weeks 25–32; PM₂.₅ small positive associations in selected early weeks).
- **G-formula, 20% reduction (Table 2)**: Natural-course cumulative risk at week 36 **6.84%** (95% CI 5.43 to 8.90), **48,863** expected preterm births.
  - **O₃**: **6.58%**; RR **0.961** (0.945 to 0.977); RD **−0.27** pp (−0.44 to −0.14); PAF **3.88%** (2.32 to 5.60); about **1,900** fewer preterm births.
  - **NO₂**: **6.96%**; RR **1.017** (1.001 to 1.030); RD **+0.11** pp (0.01 to 0.21); PAF **−1.67%** (−3.00 to −0.07).
  - **PM₂.₅**: **6.85%**; RR **1.001** (0.988 to 1.010); RD **+0.01** pp (−0.09 to 0.07); PAF **−0.08%** (−0.98 to 1.21).
- **Timing (Figure 5)**: Single-week 20% O₃ reductions lowered cumulative risk for all intervention weeks, largest late in pregnancy (e.g., week 34: **−0.046** pp at week 36). NO₂ single-week reductions mostly increased risk slightly; PM₂.₅ risk differences were near zero.

---

## :information_source: Additional analyses in this repository (not in the main manuscript)

The repository implements the full data pipeline and many **secondary and sensitivity analyses** that support the supplement or internal robustness checks. They are **not** summarized in the main-text tables or figures above.

| Topic | Scripts (examples) | Outputs (examples) |
|--------|-------------------|---------------------|
| Single-pollutant weekly DLM | `9.0`, `9.1` | `02_Output/Models/DLM_models_krg.png`, `DLM_cox_*` |
| Multipollutant weekly DLM (manuscript Figure 3) | `13.0`, `13.1` | `DLM_multi_models_krg.png`, `DLM_multi_cox_*` |
| Period-average Cox (trimester / full pregnancy) | `11.x`, `15.x` | `02_Output/Exposure_PO/`, `Exposure_PO_Multi/` |
| G-computation, 17 scenarios (caps + 10/20/30%), single-pollutant weekly models | `10.0`–`10.4` | `02_Output/G-Form/` |
| G-computation, **20% multipollutant** (manuscript Table 2, Figures 4–5) | `14.0`–`14.4` | `02_Output/G-Form-Multi/` |
| G-computation on trimester- or full-pregnancy exposure aggregates | `12.x`, `16.x` | `02_Output/G-Form-Period/`, `G-Form-Period-Multi/` |
| Exposure descriptives beyond Figure 2 (maps, time series, annual tables) | `7.0`, `7.1`, `8.0` | `02_Output/Descriptives/` |
| Positivity diagnostics by gestational week | `7.2` | `Table_positivity_by_week.xlsx` |
| IDW vs kriging comparisons | All `*_idw*` outputs | Supplementary Figures S1, S3, S6; Tables S8, S10 |

**Practical note:** Default code settings may use **500** bootstrap replicates (`GFORM_DEFAULTS$boot_iter` in `10.1`); the manuscript reports **250**. Reproduce published CIs with `GFORM_BOOT_ITER=250`. Cap-threshold interventions (`min(X, c)`) and 10%/30% proportional reductions are implemented under `10.x` but are **not** part of the main-text g-computation results (only 20% multipollutant reductions are).

---

## ![R](https://skillicons.dev/icons?i=r) Code Structure

### Setup Scripts

- `00_Code/0.1 Settings.R`: global settings and locale
- `00_Code/0.2 Packages.R`: package installation and loading (main pipeline)
- `00_Code/0.2 Packages_gform.R`: packages for g-computation pipeline
- `00_Code/0.3 Functions.R`: custom helper functions

### Data Processing Scripts

- `00_Code/1.0 Pollution_process_data.R`: load and clean interpolated PM₂.₅, NO₂, O₃ series
- `00_Code/2.0 Births_process_data.R`: birth data cleaning, cohort definition, exclusions
- `00_Code/3.0 NDVI_EarthEngine_commune_extraction.py`: NDVI extraction (Google Earth Engine)
- `00_Code/3.1 Temp_NDVI_data.R`: temperature and NDVI processing
- `00_Code/4.0 Climate_data_generate.R`: climate data generation
- `00_Code/5.0 Exposure_data_births.R`: weekly gestational exposure histories
- `00_Code/6.0 Join_full_data.R`: merge pollution, climate, and birth data
- `00_Code/8.0 Correlation_pollulants.R`: pollutant correlation analysis

### Descriptive Analysis

- `00_Code/7.0 Descriptive_births.R`: birth and preterm trends
- `00_Code/7.1 Descriptive_exposition.R`: exposure descriptives (Table S2 source, Figure 2, repository maps and time series)
- `00_Code/7.2 Positivity_analysis.R`: positivity by gestational week

### Statistical Models (manuscript and supplement)

- `00_Code/9.0 DLM_pollution.R`, `9.1 DLM_plots.R`: single-pollutant weekly DLMs (supplement)
- `00_Code/13.0 DLM_multi_pollution.R`, `13.1 DLM_multi_plots.R`: **multipollutant weekly DLMs (Figure 3)**
- `00_Code/11.x`, `15.x`: period-average exposure Cox models (Tables S5–S6)

### G-Computation: manuscript (multipollutant, 20%)

- `00_Code/14.0 G-Form_multi_functions.R`: multipollutant weekly Cox + g-formula
- `00_Code/14.1 G-Form_multi_build_interventions.R`: counterfactual histories (pct20 × 3 pollutants)
- `00_Code/14.2 G-Form_multi_models.R`: models, bootstrap, single-week heatmaps
- `00_Code/14.3 G-Form_multi_plots.R`: **Figures 4–5** (`02_Output/G-Form-Multi/Figures/`)
- `00_Code/14.4 G-Form_multi_table.R`: **Table 2** summary exports

### G-Computation: extended scenarios (repository only)

- `00_Code/10.0`–`10.4`: 17 global scenarios (caps and 10/20/30% reductions) under single-pollutant weekly models (`02_Output/G-Form/`)
- `00_Code/12.x`, `16.x`: g-formula on trimester/full-pregnancy exposure aggregates
- `00_Code/intervention.txt`: scenario list reference

### G-Computation Intervention Registry (`10.x`, repository extension)

Scenarios are defined in `GFORM_INTERVENTION_REGISTRY` (`00_Code/10.1 G-Form_functions.R`). Stage 1 (`10.0`) writes one RDS per scenario to `02_Output/G-Form/Interventions/`; Stage 2 (`10.2`) runs models and bootstrap by `intervention_number` (1–17) or via `GFORM_INTERVENTIONS`. **These 17 scenarios are not reported in the main manuscript** (see table for stubs used in sensitivity work and composite figures under `G-Form/Figures/`).

| # | Pollutant | Scenario | Registry ID | Output stub |
|---|-----------|----------|-------------|-------------|
| 1 | PM₂.₅ | −20% all weeks | `pm25_krg_pct20` | `pm25_pct20` |
| 2 | NO₂ | −20% all weeks | `no2_krg_pct20` | `no2_pct20` |
| 3 | O₃ | −20% all weeks | `o3_krg_pct20` | `o3_pct20` |
| 4 | PM₂.₅ | < 20 µg/m³ | `pm25_krg_lt20` | `pm25_lt20` |
| 5 | PM₂.₅ | < 5 µg/m³ | `pm25_krg_lt5` | `pm25_lt5` |
| 6 | NO₂ | < 20 ppbv | `no2_krg_lt20` | `no2_lt20` |
| 7 | NO₂ | < 5 ppbv | `no2_krg_lt5` | `no2_lt5` |
| 8 | PM₂.₅ | < 15 µg/m³ | `pm25_krg_lt15` | `pm25_lt15` |
| 9 | PM₂.₅ | < 10 µg/m³ | `pm25_krg_lt10` | `pm25_lt10` |
| 10 | NO₂ | < 15 ppbv | `no2_krg_lt15` | `no2_lt15` |
| 11 | NO₂ | < 10 ppbv | `no2_krg_lt10` | `no2_lt10` |
| 12 | PM₂.₅ | −10% all weeks | `pm25_krg_pct10` | `pm25_pct10` |
| 13 | PM₂.₅ | −30% all weeks | `pm25_krg_pct30` | `pm25_pct30` |
| 14 | NO₂ | −10% all weeks | `no2_krg_pct10` | `no2_pct10` |
| 15 | NO₂ | −30% all weeks | `no2_krg_pct30` | `no2_pct30` |
| 16 | O₃ | −10% all weeks | `o3_krg_pct10` | `o3_pct10` |
| 17 | O₃ | −30% all weeks | `o3_krg_pct30` | `o3_pct30` |

*Cap semantics*: `< 5` means exposure is fixed at 5 units when the observed weekly value exceeds 5 (`min(X, c)`), not a subtraction of 5 from observed concentrations.

---

## :chart_with_upwards_trend: Manuscript Figures and Tables (September 2026)

### Figure 1. Flowchart of Analytical Sample Construction (2010–2020)

Starting from **2,557,140** singleton births in Chile (2010–2020), sequential exclusions yielded **713,918** births in urban Santiago (**51,081** preterm, **7.2%**). See `03_Paper/Gformula_Santiago_PTB_29092026.docx`.

### Table 1. Cohort Characteristics

Descriptive statistics for the full analytic cohort and preterm births (sex, parental characteristics, SOVI). Generated from the cleaned cohort (`2.0`); values in the manuscript Table 1.

### Figure 2. Distribution of Daily PM₂.₅, NO₂, and O₃ (Kriging)

![](/02_Output/Descriptives/Map_Exposure_daily_mean_and_density_KRG.png)

*Note*: Municipality-day concentrations from ordinary kriging. IDW analogue: Supplementary Figure S1 (`Histogram_IDW_panel_compiled.png`).

### Figure 3. Hazard Ratios for Preterm Birth by Gestational Week (Multipollutant DLM)

![](/02_Output/Models/DLM_multi_models_krg.png)

*Note*: Multipollutant distributed-lag Cox models (kriging); HRs per 1-unit increase, week *t* adjusted for concurrent co-pollutants and time-weighted lag through *t*−1. N = 713,918. Single-pollutant plots: `DLM_models_krg.png` (Figure S6 / Tables S8, S10).

### Table 2. G-Formula Estimates Under 20% Weekly Exposure Reduction

| Scenario | Prevalence at week 36 (95% CI) | RR (95% CI) | RD, pp (95% CI) | PAF, % (95% CI) |
|----------|--------------------------------|-------------|-----------------|-----------------|
| Natural course | 6.84 (5.43; 8.90) | 1.00 | 0.00 | 0.00 |
| PM₂.₅ −20% | 6.85 (5.43; 8.93) | 1.001 (0.988; 1.010) | 0.01 (−0.09; 0.07) | −0.08 (−0.98; 1.21) |
| NO₂ −20% | 6.96 (5.63; 9.36) | 1.017 (1.001; 1.030) | 0.11 (0.01; 0.21) | −1.67 (−3.00; −0.07) |
| O₃ −20% | 6.58 (5.22; 8.86) | 0.961 (0.945; 0.977) | −0.27 (−0.44; −0.14) | 3.88 (2.32; 5.60) |

*Note*: Multipollutant parametric g-computation; kriging exposures; 250 bootstrap replicates. Excel exports: `02_Output/G-Form-Multi/Summary_results/*_pct20_point_estimates.xlsx`.

### Figure 4. Cumulative Preterm Birth Risk (20% Reduction, Multipollutant G-Formula)

![](/02_Output/G-Form-Multi/Figures/Figure_cumulative_risk_interventions_percent.png)

*Note*: Natural course vs 20% proportional reduction applied to all gestational weeks for PM₂.₅, NO₂, and O₃ (panels A–C). Follow-up weeks 28–36.

### Figure 5. Critical-Window Heatmap (Single-Week 20% Reduction)

![](/02_Output/G-Form-Multi/Figures/Figure_heatmap_rd_interventions_percent.png)

*Note*: Risk difference by intervention week (columns) and follow-up week (rows). Multipollutant models; point estimates; 20% reduction only.

### Supplementary material (reported in `Gformula_Santiago_PTB_29092026_Supp.docx`)

- **Methods S1**, **Table S1**: g-formula estimators and bootstrap.
- **Table S2**: municipality exposure descriptives (from `7.1`).
- **Figures S1–S5**: IDW distributions, temporal trends, kriging vs IDW, pollutant–covariate correlations.
- **Tables S3–S4**: weekly IQRs and risk-window descriptives.
- **Tables S5–S6**, **S7–S10**, **Figure S6**: period-average and single-pollutant / IDW DLM results.

### Repository-only descriptive figures (not in the main PDF)

These support exploration and the supplement but are **not** numbered manuscript figures:

- Preterm trends: `Preterm_trends_2010_2020.png` (`7.0`)
- Spatial means and 20% density overlays: `Map_Exposure_daily_mean_and_density_KRG.png` (`7.1`)
- Regional daily time series: `Time_distribution_pm25_no2_o3.png` (`7.1`)

---

## :file_folder: Data Availability

### Input Data Sources

1. **Birth Records**: Chilean Ministry of Health (DEIS) vital statistics (2010–2020)
   - Location: `01_Data/Input/Nacimientos/`
   - Variables: Gestational age, birth weight, parental characteristics, municipality of residence

2. **Air Pollution Data**: National air-quality monitoring network
   - Location: `01_Data/Input/Clime_series/`
   - Pollutants: PM₂.₅ (beta-attenuation), NO₂ (chemiluminescence), O₃ (UV photometry)
   - Interpolation: Ordinary kriging (primary) and IDW (sensitivity)
   - Spatial unit: Municipal administrative centers (33 urban comunas in the analytic cohort)

3. **Temperature Data**: CR2MET gridded climate product
   - Processed in: `00_Code/4.0 Climate_data_generate.R`
   - Variable: Daily mean ambient temperature (TAD)

4. **NDVI**: MODIS MOD13Q1 (250 m, 16-day composite)
   - Extraction: `00_Code/3.0 NDVI_EarthEngine_commune_extraction.py`
   - Gaps imputed with Kalman smoother

5. **Socioeconomic Vulnerability Index (SOVI)**
   - Location: `01_Data/Input/SOVI/`
   - Categories: Low, medium-low, medium-high

6. **Municipal Boundaries**
   - Location: `01_Data/Input/district_geo/`

### Processed Datasets

Main analytical datasets are stored in `01_Data/Output/`:

- `births_2010_2020.RData`: cleaned birth records
- `Contamination_Climate_Data_2010_2020.RData`: merged pollution and climate series
- `births_2010_2020_exposure_weeks.RData`: weekly gestational exposure histories
- `births_2010_2020_exposure_weeks_lagged.RData`: weekly data with DLM lag terms

Descriptive outputs from `7.1` are stored in `02_Output/Descriptives/`:

- `Table_exposure_commune_PM25_O3_summary.xlsx`: municipality min/mean/max (Table S2)
- `Table_Annual_Summary_Contaminants_Estimators.xlsx`: annual pollutant summaries
- `Table_Municipality_Summary_Contaminants_Estimators.xlsx`: municipality summaries across estimators
- `Map_Exposure_daily_mean_and_density_KRG.png` / `..._IDW.png`: repository spatial panels (`7.1`)
- `Time_distribution_pm25_no2_o3.png`: daily regional means (`7.1`)
- `Histogram_*`: Figure 2 and Supplementary Figure S1

**Manuscript g-computation** (`14.x`): `02_Output/G-Form-Multi/` (`Summary_results/`, `Figures/`, `Heatmap/`, `Bootstrap/`).

**Extended g-computation** (`10.x`): `02_Output/G-Form/` (17 scenarios, single-pollutant weekly models):

- `Summary_results/`: point estimates and bootstrap CIs (`{stub}_point_estimates.xlsx`)
- `Interventions/`: counterfactual exposure histories (RDS)
- `WeeklyEffects/`, `PopulationEffects/`: detailed effect objects
- `Bootstrap/{stub}/`, `Heatmap/{stub}/`: replicates and single-week maps

Output stubs follow pollutant and scenario: e.g. `pm25_pct10`, `pm25_pct20`, `pm25_pct30`, `pm25_lt20`, `no2_pct30`, `o3_pct10`.

**Note**: Individual-level birth records cannot be publicly shared due to Chilean data protection regulations. Aggregated results and analysis code are available in this repository.

---

## :computer: Reproducibility

### System Requirements

- R ≥ 4.0.0
- Python 3 (for NDVI extraction via Google Earth Engine)
- Recommended: ≥ 16 GB RAM; Linux server for parallel g-computation

### Required R Packages

Automatically installed via `00_Code/0.2 Packages.R` and `00_Code/0.2 Packages_gform.R`:

- **Data manipulation**: `tidyverse`, `data.table`, `janitor`, `rio`
- **Spatial analysis**: `chilemapas`, `sf`, `rnaturalearth`, `maptiles`, `tidyterra`
- **Survival analysis**: `survival`, `flexsurv`, `survminer`
- **Distributed lag / splines**: `dlnm`, `splines`, `mgcv`
- **Parallel computing**: `future`, `furrr`, `doParallel`
- **Visualization**: `ggplot2`, `patchwork`, `ggpubr`, `ggspatial`, `RColorBrewer`, `ragg`, `scales`
- **Imputation**: `imputeTS`, `zoo`

### Running the Analysis

1. **Setup**:
   ```r
   source("00_Code/0.1 Settings.R")
   source("00_Code/0.2 Packages.R")
   source("00_Code/0.3 Functions.R")
   ```

2. **Data processing** (run in order):
   ```r
   source("00_Code/1.0 Pollution_process_data.R")
   source("00_Code/2.0 Births_process_data.R")
   source("00_Code/3.1 Temp_NDVI_data.R")
   source("00_Code/4.0 Climate_data_generate.R")
   source("00_Code/5.0 Exposure_data_births.R")
   source("00_Code/6.0 Join_full_data.R")
   ```

3. **Descriptive analysis**:
   ```r
   source("00_Code/7.0 Descriptive_births.R")
   source("00_Code/7.1 Descriptive_exposition.R")
   source("00_Code/8.0 Correlation_pollulants.R")
   ```

4. **Distributed-lag models (manuscript Figure 3)**:
   ```r
   source("00_Code/9.0 DLM_pollution.R")   # prerequisite lags
   source("00_Code/13.0 DLM_multi_pollution.R")
   source("00_Code/13.1 DLM_multi_plots.R")
   ```

5. **G-computation (manuscript Table 2, Figures 4–5)**:
   ```r
   source("00_Code/14.1 G-Form_multi_build_interventions.R")
   source("00_Code/14.2 G-Form_multi_models.R")   # set GFORM_BOOT_ITER=250 to match the paper
   source("00_Code/14.3 G-Form_multi_plots.R")
   source("00_Code/14.4 G-Form_multi_table.R")
   ```

6. **Extended g-computation (repository only, optional)**:
   ```r
   source("00_Code/10.0 G-Form_build_interventions.R")
   source("00_Code/10.2 G-Form_models.R")
   source("00_Code/10.3 G-Form_plots.R")
   source("00_Code/10.4 G-Form_table.R")
   ```

   For server/parallel execution:
   ```bash
   GFORM_EXEC_MODE=server Rscript "00_Code/10.2 G-Form_models.R"
   ```

   Run a subset of interventions (by number 1–17):
   ```bash
   GFORM_INTERVENTIONS=12,13,16 Rscript "00_Code/10.2 G-Form_models.R"
   ```

### Notes on Computation Time

- **Birth data processing** (`2.0`): moderate (depends on raw file size)
- **Weekly exposure expansion** (`5.0`, `6.0`): several hours (large longitudinal dataset)
- **DLM Cox models** (`9.0`): ~20–30 minutes per pollutant/method
- **G-computation bootstrap** (`10.2`): several hours to days (250–500 bootstrap replicates; parallelized on server)
- **Total pipeline**: plan for multi-hour to overnight runs on a modern workstation or Linux server

Detailed methodological notes: `02_Output/Notas_G-Formula_resultados.md`

---

## :open_book: Codebook

### Birth Variables

- `id`: Unique birth identifier
- `com`: Municipality code
- `name_com`: Municipality name
- `weeks`: Gestational age at delivery (weeks)
- `date_nac`: Date of birth
- `sex`: Infant sex (Boy/Girl)
- `tbw`: Birth weight (grams)
- `birth_preterm`: Preterm birth indicator (<37 weeks)
- `birth_very_preterm`: Very preterm (28–31 weeks)
- `birth_moderately_preterm`: Moderate preterm (32–33 weeks)
- `birth_late_preterm`: Late preterm (34–36 weeks)

### Parental and Context Variables

- `age_group_mom`, `educ_group_mom`, `job_group_mom`: Maternal age, education, employment
- `age_group_dad`, `educ_group_dad`, `job_group_dad`: Paternal age, education, employment
- `month_week1`, `year_week1`: Month and year of last menstrual period
- `covid`: COVID-19 period indicator
- `vulnerability`: SOVI category (Low, Medium-low, Medium-high)

### Exposure Variables

- `pm25_krg`, `no2_krg`, `o3_krg`: Weekly kriging-interpolated concentrations
- `pm25_idw`, `no2_idw`, `o3_idw`: Weekly IDW-interpolated concentrations (sensitivity)
- `tad`: Weekly mean ambient temperature
- `ndvi_full`: Municipality-level NDVI (full pregnancy average)
- Lag term (`Liw`): Time-weighted cumulative lag through prior gestational weeks

---

## :microscope: Methods Detail

### Distributed-Lag Exposure

For gestational week \(w \geq 2\):

\[
L_{iw} = \sum_{s=1}^{w-1} \frac{X_{is}}{w - s}
\]

### Exclusion Criteria

Births were excluded if:

- Outside urban Metropolitan Santiago (33 comunas: Province of Santiago plus Puente Alto)
- Missing date of birth, gestational age, or municipality
- Maternal age <12 or >50 years
- Gestational age <28 weeks
- Multiple births
- Missing covariates
- Implausible birthweight-for-gestational-age (Alexander et al., 1996)
- Fixed-cohort bias: gestational window not fully observed within 2010–2020

### G-Computation Interventions

**Main manuscript:** proportional reduction with \(\pi = 0.20\) on all gestational weeks of the target pollutant, multipollutant natural-course models, co-pollutants fixed at observed values:

\[
X'_{iw} = X_{iw} \times (1 - \pi), \quad \pi = 0.20
\]

**Single-week reduction (Figure 5):** the same 20% rule applied only in week \(j\); all other weeks remain observed; lag terms recomputed.

**Repository extensions (`10.x`):** additional \(\pi \in \{0.10, 0.30\}\) and cap rules \(X'_{iw} = \min(X_{iw}, c)\) with \(c \in \{5, 10, 15, 20\}\) µg/m³ (PM₂.₅) or ppbv (NO₂), fitted under single-pollutant weekly Cox models.

Population metrics at week 36: prevalence, expected cases, RR, RD, AR, and PAF (definitions in supplementary Table S1).

---

## :file_cabinet: Repository Structure

```
Contamination_PTB_G-Formula/
├── 00_Code/                        # Analysis scripts
│   ├── 0.1–0.3                     # Settings, packages, functions
│   ├── 1.0–8.0                     # Data processing and descriptives
│   ├── 9.0–9.1, 11.x, 13.0–13.1, 15.x  # DLM and period Cox models
│   ├── 10.0–10.4                   # Extended g-computation (17 scenarios)
│   ├── 12.x, 14.x, 16.x            # Period and multipollutant g-computation
│   └── old_code/                   # Archived scripts
├── 01_Data/
│   ├── Input/                      # Raw data (not publicly available)
│   └── Output/                     # Processed analytical datasets
├── 02_Output/
│   ├── Descriptives/               # Tables, histograms, maps, time series (see 7.1)
│   │   └── assets/                 # Map figure assets (e.g. warning icon)
│   ├── Models/                     # DLM results and figures
│   ├── G-Form/                     # Extended g-computation (10.x)
│   ├── G-Form-Multi/               # Manuscript g-computation (14.x)
│   ├── G-Form-Period/              # Period-exposure g-computation (12.x)
│   ├── G-Form-Period-Multi/        # Period multipollutant g-computation (16.x)
│   ├── Exposure_PO/                # Period-average Cox (11.x)
│   └── Exposure_PO_Multi/          # Multipollutant period Cox (15.x)
├── 03_Paper/                       # Manuscript and supplementary material
│   ├── Gformula_Santiago_PTB_29092026.docx
│   └── Gformula_Santiago_PTB_29092026_Supp.docx
├── 04_Conference/                  # Conference abstracts
└── README.md
```

---

## :warning: Important Notes

### Data Privacy

Individual-level birth records are confidential and cannot be shared publicly. Researchers interested in data access should contact the Chilean Ministry of Health (DEIS).

### Air Quality and Climate Data

- National air-quality network: Chilean Ministry of Environment
- CR2MET: [Center for Climate and Resilience Research (CR²)](http://www.cr2.cl/datos-productos-grillados/)

### Citation

If you use this code or methodology, please cite:

> Blanco, E., Conejeros, J.D., Bravo, I., Cornejo, F., Osses, A., & Benmarhnia, T. Estimating the perinatal health benefits of hypothetical pollution interventions in Santiago, Chile using parametric g-computation. *Under Review*. 2026.

---

## :email: Contact

For questions about the code or methodology:

- **Estela Blanco**: <estela.blanco@uc.cl>
- **José Daniel Conejeros**: <jdconejeros@uc.cl>

For data access inquiries:

- Chilean Ministry of Health: [https://www.minsal.cl](https://www.minsal.cl)

---

## :page_facing_up: License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

---

## :handshake: Acknowledgments

This research was supported by FONDECYT de Iniciación en Investigación Nº 11240322 and the Center for Climate and Resilience Research (CR²), FONDAP/ANID 1523A0002. We thank the Chilean Ministry of Health (DEIS) for access to birth records, the national air-quality monitoring network for pollution data, and CR² for climate data.

**Data sources**:

- Birth records: DEIS, Chilean Ministry of Health
- Air pollution: National air-quality monitoring network (SINCA)
- Temperature: CR2MET v2.5, Center for Climate and Resilience Research, Universidad de Chile
- NDVI: MODIS MOD13Q1 via Google Earth Engine
