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

:mailbox_with_mail: **Estela Blanco** (<estela.blanco@uc.cl>) — **Principal Investigator / Corresponding Author**

:mailbox_with_mail: **José Daniel Conejeros** (<jdconejeros@uc.cl>) — **Research Assistant / Repository Manager**

**Research Collaborators**: Ismael Bravo, Felipe Cornejo, Axel Osses, & Tarik Benmarhnia

## :pushpin: Publication

*Work in progess*. 

---

## :dart: Project Overview

### Background

Ambient air pollution is a major environmental risk factor for adverse perinatal outcomes. Preterm birth (delivery before 37 completed weeks of gestation) remains a leading cause of neonatal mortality and long-term morbidity. While epidemiological evidence links PM₂.₅, NO₂, and O₃ to preterm birth, fewer studies have translated observed exposure–response associations into population-level estimates of the health benefits of realistic pollution-reduction scenarios, particularly in Latin American urban settings with heterogeneous exposure patterns.

### Objective

To estimate the population-level impact of hypothetical reductions in PM₂.₅, NO₂, and O₃ on the cumulative risk of preterm birth among singleton births in urban Santiago, Chile (2010–2020), using distributed-lag Cox models combined with parametric g-computation.

### Methods

We conducted a population-based retrospective cohort study using:

- **Study population**: Singleton live births to women residing in the urban conurbation of Santiago (33 municipalities of the Province of Santiago plus Puente Alto), 2010–2020
- **Sample size**: 713,918 births after exclusion criteria (51,081 preterm; 7.2%)
- **Exposure**: Weekly municipality-level PM₂.₅, NO₂, and O₃ from the national air-quality monitoring network, spatially interpolated to municipal centroids using **ordinary kriging** (primary) and **inverse distance weighting (IDW)** (sensitivity)
- **Additional exposures**: Daily mean temperature (CR2MET), NDVI (MODIS MOD13Q1 via Google Earth Engine), and municipal social vulnerability index (SOVI)
- **Exposure windows**: Weekly gestational exposure series (weeks 1–44); period-specific averages for descriptives (full pregnancy, trimesters, last 30 and last 4 days)
- **Outcome**: Preterm birth (<37 weeks), with subcategories:
  - Very preterm (28–31 weeks)
  - Moderate preterm (32–33 weeks)
  - Late preterm (34–36 weeks)
- **Distributed-lag models (DLM)**: Week-specific Cox proportional hazards models with time-weighted lagged exposure and delayed entry at week 28
- **G-computation**: Parametric g-formula applied to natural-course Cox models (weeks 28–36) under counterfactual exposure histories
- **Intervention scenarios**:
  - **Threshold (cap)**: PM₂.₅ capped at 5, 10, 15, 20 µg/m³; NO₂ capped at 5, 10, 15, 20 ppbv
  - **Proportional reduction**: 20% reduction applied to PM₂.₅, NO₂, or O₃ across all gestational weeks
  - **Single-week interventions**: 20% reduction in one gestational week at a time (critical-window map)
- **Covariates**: Newborn sex; maternal and paternal age, education, and occupation; month and year of last menstrual period; COVID-19 period; SOVI; weekly temperature; NDVI
- **Inference**: Parametric bootstrap (250 replicates in manuscript; configurable in code via `GFORM_BOOT_ITER`)

### Key Findings

- **NO₂ interventions** showed the most consistent protective effects on cumulative preterm birth risk. A 20% reduction across gestation was associated with a risk difference of −0.12 percentage points (PAF −1.77%). Stronger benefits were observed under lower caps (e.g., NO₂ < 5 ppbv: RD −0.48 pp; PAF −6.98%).
- **O₃**: A uniform 20% reduction was associated with a risk difference of −0.26 percentage points (PAF 3.80%), suggesting a meaningful share of observed preterm risk may be attributable to ozone exposure under the natural course.
- **PM₂.₅**: Threshold and proportional-reduction scenarios showed small and statistically uncertain effects, with confidence intervals spanning both protective and adverse directions.
- **Timing matters**: Single-week intervention heatmaps (Figure 5) identified gestational windows in which exposure reductions would yield the largest expected impact on cumulative preterm risk.

---

## ![R](https://skillicons.dev/icons?i=r) Code Structure

### Setup Scripts

- `00_Code/0.1 Settings.R` — Global settings and locale
- `00_Code/0.2 Packages.R` — Package installation and loading (main pipeline)
- `00_Code/0.2 Packages_gform.R` — Packages for g-computation pipeline
- `00_Code/0.3 Functions.R` — Custom helper functions

### Data Processing Scripts

- `00_Code/1.0 Pollution_process_data.R` — Load and clean interpolated PM₂.₅, NO₂, O₃ series
- `00_Code/2.0 Births_process_data.R` — Birth data cleaning, cohort definition, exclusions
- `00_Code/3.0 NDVI_EarthEngine_commune_extraction.py` — NDVI extraction (Google Earth Engine)
- `00_Code/3.1 Temp_NDVI_data.R` — Temperature and NDVI processing
- `00_Code/4.0 Climate_data_generate.R` — Climate data generation
- `00_Code/5.0 Exposure_data_births.R` — Weekly gestational exposure histories
- `00_Code/6.0 Join_full_data.R` — Merge pollution, climate, and birth data
- `00_Code/8.0 Correlation_pollulants.R` — Pollutant correlation analysis

### Descriptive Analysis

- `00_Code/7.0 Descriptive_births.R` — Birth and preterm trends
- `00_Code/7.1 Descriptive_exposition.R` — Exposure descriptives

### Statistical Models

- `00_Code/9.0 DLM_pollution.R` — Distributed-lag Cox models (PM₂.₅, NO₂, O₃)
- `00_Code/9.1 DLM_plots.R` — DLM visualization

### G-Computation Pipeline

- `00_Code/10.0 G-Form_build_interventions.R` — Build counterfactual exposure histories (Stage 1)
- `00_Code/10.1 G-Form_functions.R` — Core g-formula functions
- `00_Code/10.2 G-Form_models.R` — Run interventions, bootstrap, and heatmaps (Stage 2)
- `00_Code/10.3 G-Form_plots.R` — Publication figures (cumulative risk, heatmaps)
- `00_Code/10.4 G-Form_table.R` — Summary tables (Table 3)

---

## :chart_with_upwards_trend: Principal Findings

### Figure 1. Flowchart of Analytical Sample Construction (2010–2020)

Starting from 2,557,140 singleton births in Chile (2010–2020), sequential exclusions yielded a final analytic sample of **713,918 births** in urban Metropolitan Santiago.

*Note*: Flowchart included in manuscript (`03_Paper/CTA_PTB_Gf_Manuscript.docx`).

### Figure 2. Distribution of Daily PM₂.₅, NO₂, and O₃ (Kriging)

![](/02_Output/Descriptives/Histogram_KRG_panel_compiled.png)

*Note*: Municipality-day pollutant concentrations interpolated by ordinary kriging. Histograms summarize overall, seasonal, and pollutant-specific distributions across the study period.

### Preterm Birth Trends in Santiago (2010–2020)

![](/02_Output/Descriptives/Preterm_trends_2010_2020.png)

*Note*: Annual prevalence of preterm birth and subcategories per 1,000 births. Analytic cohort N = 713,918.

### Figure 3. Hazard Ratios for Preterm Birth by Gestational Week (DLM)

![](/02_Output/Models/DLM_models_krg.png)

*Note*: Week-specific hazard ratios (HR) and 95% CIs from distributed-lag Cox models. Each point represents the acute effect of weekly exposure at gestational week *t*, conditional on time-weighted lagged exposure through week *t*−1. Models adjusted for sex, parental characteristics, calendar time, SOVI, temperature, and NDVI. Kriging-based exposures; N = 713,918.

### Figure 4. Cumulative Preterm Birth Risk Under 20% Pollution Reduction Scenarios

![](/02_Output/G-Form/Figures/Figure_cumulative_risk_interventions.png)

*Note*: Cumulative risk of preterm birth (weeks 28–36) under the natural course (observed exposure) vs. a uniform 20% reduction in weekly PM₂.₅, NO₂, or O₃. Estimates from parametric g-computation with distributed-lag Cox models.

### Figure 5. Critical-Window Heatmap of Risk Differences (Single-Week 20% Reduction)

![](/02_Output/G-Form/Figures/Figure_heatmap_rd_interventions.png)

*Note*: Risk difference (RD) in cumulative preterm birth risk for a 20% reduction applied in a single gestational week (columns) and evaluated through follow-up weeks 28–36 (rows). Marginal single-week interventions; not directly comparable to simultaneous full-history interventions.

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
   - Spatial unit: Municipality centroids (34 comunas)

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

- `births_2010_2020.RData` — Cleaned birth records
- `Contamination_Climate_Data_2010_2020.RData` — Merged pollution and climate series
- `births_2010_2020_exposure_weeks.RData` — Weekly gestational exposure histories
- `births_2010_2020_exposure_weeks_lagged.RData` — Weekly data with DLM lag terms

G-computation outputs are stored in `02_Output/G-Form/`:

- `Summary_results/` — Point estimates and bootstrap CIs by scenario
- `Interventions/` — Counterfactual exposure histories (RDS)
- `WeeklyEffects/`, `PopulationEffects/` — Detailed effect objects

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
- **Spatial analysis**: `chilemapas`, `sf`, `rnaturalearth`
- **Survival analysis**: `survival`, `flexsurv`, `survminer`
- **Distributed lag / splines**: `dlnm`, `splines`, `mgcv`
- **Parallel computing**: `future`, `furrr`, `doParallel`
- **Visualization**: `ggplot2`, `patchwork`, `ggpubr`, `RColorBrewer`
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

4. **Distributed-lag models**:
   ```r
   source("00_Code/9.0 DLM_pollution.R")
   source("00_Code/9.1 DLM_plots.R")
   ```

5. **G-computation** (two stages):
   ```r
   # Stage 1: build intervention objects (run once)
   source("00_Code/10.0 G-Form_build_interventions.R")

   # Stage 2: run models, bootstrap, and heatmaps
   source("00_Code/10.2 G-Form_models.R")
   source("00_Code/10.3 G-Form_plots.R")
   source("00_Code/10.4 G-Form_table.R")
   ```

   For server/parallel execution:
   ```bash
   GFORM_EXEC_MODE=server Rscript "00_Code/10.2 G-Form_models.R"
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

- Outside urban Metropolitan Santiago (33 comunas + Puente Alto)
- Missing date of birth, gestational age, or municipality
- Maternal age <12 or >50 years
- Gestational age <28 weeks
- Multiple births
- Missing covariates
- Implausible birthweight-for-gestational-age (Alexander et al., 1996)
- Fixed-cohort bias: gestational window not fully observed within 2010–2020

### G-Computation Interventions

**Threshold (cap)**:

\[
X'_{iw} = \min(X_{iw}, c)
\]

**Proportional reduction (20%)**:

\[
X'_{iw} = X_{iw} \times (1 - 0.20)
\]

**Single-week reduction**: Apply 20% reduction only in week \(j\); all other weeks remain at observed values.

Population metrics at week 36: prevalence, expected cases, risk ratio (RR), risk difference (RD), attributable risk (AR), and population attributable fraction (PAF).

---

## :file_cabinet: Repository Structure

```
Contamination_PTB_G-Formula/
├── 00_Code/                        # Analysis scripts
│   ├── 0.1–0.3                     # Settings, packages, functions
│   ├── 1.0–8.0                     # Data processing and descriptives
│   ├── 9.0–9.1                     # Distributed-lag Cox models
│   ├── 10.0–10.4                   # G-computation pipeline
│   └── old_code/                   # Archived scripts
├── 01_Data/
│   ├── Input/                      # Raw data (not publicly available)
│   └── Output/                     # Processed analytical datasets
├── 02_Output/
│   ├── Descriptives/               # Tables and descriptive plots
│   ├── Models/                     # DLM results and figures
│   ├── G-Form/                     # G-computation outputs
│   └── idw_vs_kriging/             # Interpolation comparison
├── 03_Paper/                       # Manuscript and supplementary material
│   ├── CTA_PTB_Gf_Manuscript.docx
│   └── CTA_PTB_Gf_Supplementary_Material.docx
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

> Blanco, E., Conejeros, J.D., Bravo, I., Cornejo, F., Osses, A., & Benmarhnia, T. Estimating the perinatal health benefits of hypothetical pollution interventions in Santiago, Chile using parametric g-computation. *Under Review*. 2025.

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
