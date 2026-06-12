# Notas interpretativas — resultados G-Formula (parto pretérmino)

Documentación del pipeline `10.0` / `10.1` / `10.2`. Implementación en `00_Code/10.1 G-Form_functions.R`.

---

## 1. Contexto

| Elemento | Definición |
|----------|------------|
| **Evento** | Parto pretérmino: `birth_preterm = 1` si edad gestacional al parto `< 37` semanas |
| **Cohortes** | Nacimientos con `weeks ≥ 28` (filtro en procesamiento de datos) |
| **Modelo base** | Cox proporcional por semana de exposición, alineado con `9.0 DLM_pollution.R` |
| **Entrada al riesgo** | Semana 28: `Surv(tstart = 28, weeks, birth_preterm)` |
| **Semanas de riesgo** | `t = 1, …, 36` (`risk_weeks`); antes de la semana 28 no hay evento |
| **Contraste** | Curso natural (exposición observada) vs. escenario de intervención contrafactual |

**Outputs principales**

- Excel (estimaciones puntuales + IC): `02_Output/G-Form/Summary_results/{stub}_point_estimates.xlsx`
- Tablas Fig. 3 y 4: `02_Output/G-Form/Other/{stub}_figure3.rds`, `{stub}_figure4.rds`
- Efectos semanales / globales (RDS): `WeeklyEffects/`, `PopulationEffects/`

---

## 2. Exposición, lag e intervenciones

### 2.1 Exposición semanal

Para la persona \(i\) y la semana gestacional \(w\), denotamos \(X_{iw}\) la concentración del contaminante modelado (p. ej. PM\(_{2.5}\) en µg/m³).

### 2.2 Lag ponderado (DLM)

Para semanas \(w \geq 2\):

\[
L_{iw} = \sum_{s=1}^{w-1} \frac{X_{is}}{w - s}
\]

Solo se usan semanas previas a \(w\) (misma definición que en `9.0`).

### 2.3 Intervenciones globales (Etapa 1, script 10.0)

Se transforma la historia completa de exposición y se recalcula el lag.

**Techo (cap):**

\[
X'_{iw} = \min(X_{iw}, c)
\]

**Reducción porcentual:**

\[
X'_{iw} = X_{iw} \times (1 - \pi), \quad \pi = 0{,}20 \text{ para reducción del 20 %}
\]

**Curso natural:** \(X'_{iw} = X_{iw}\).

El lag intervenido \(L'_{iw}\) se obtiene sustituyendo \(X_{is}\) por \(X'_{is}\) en la fórmula anterior.

### 2.4 Intervención en una sola semana (Figura 4)

Para la semana de intervención \(j\):

\[
X'_{iw} =
\begin{cases}
X_{iw}(1-\pi) & \text{si } w = j \\
X_{iw} & \text{si } w \neq j
\end{cases}
\]

---

## 3. Ajuste Cox (una sola vez, curso natural)

Para cada semana de riesgo \(t \geq 28\) se ajusta un modelo independiente usando la exposición **observada**:

\[
h_i(u \mid t) = h_{0t}(u) \exp\left(\eta_{it}\right)
\]

\[
\eta_{it} = \beta_{1t}\, X_{it} + \beta_{2t}\, L_{it} + \beta_{3t}\, \mathrm{TAD}_{it} + \boldsymbol{\gamma}_t^\top \mathbf{Z}_i
\]

donde:

- \(u\) = edad gestacional al parto (`weeks`)
- Entrada retardada desde la semana 28
- \(\mathbf{Z}_i\) = covariables fijas (sexo, edad/educación/ocupación de padres, mes/año inicio embarazo, COVID, vulnerabilidad, TAD basal, NDVI)
- Los coeficientes \(\hat{\boldsymbol{\beta}}_t\) **no se re-estiman** bajo intervención; solo cambian las covariables de exposición en la predicción (g-computación estándar)

El predictor lineal estimado para la persona \(i\) en la semana \(t\) bajo un escenario \(s \in \{\mathrm{nat}, \mathrm{int}\}\):

\[
\widehat{\mathrm{lp}}_{it}^{(s)} = \hat{\eta}_{it}^{(s)}
\]

evaluado con \(X_{it}^{(s)}\), \(L_{it}^{(s)}\) del escenario correspondiente y el mismo \(\hat{\boldsymbol{\beta}}_t\).

---

## 4. De Cox a probabilidad semanal (g-computación discreta)

### 4.1 Incremento del hazard baselina

Del modelo Cox de la semana \(t\), se extrae la función de riesgo acumulada baselina \( \hat{H}_{0t}(u) \). El incremento entre \(t-1\) y \(t\):

\[
\Delta \hat{H}_{0t}(t) = \hat{H}_{0t}(t) - \hat{H}_{0t}(t-1)
\]

(con valor mínimo 0 si no hay incremento).

### 4.2 Probabilidad de evento en la semana \(t\)

\[
\hat{p}_i(t \mid s) = 1 - \exp\left(-\Delta \hat{H}_{0t}(t) \cdot \exp\left(\widehat{\mathrm{lp}}_{it}^{(s)}\right)\right)
\]

### 4.3 Probabilidad de no evento

\[
\hat{q}_i(t \mid s) = 1 - \hat{p}_i(t \mid s)
\]

**Semanas \(t < 28\):** \(\hat{p}_i(t \mid s) = 0\), \(\hat{q}_i(t \mid s) = 1\) (aún no entra al riesgo).

---

## 5. Supervivencia y riesgo acumulado individual

Para cada persona \(i\) y escenario \(s\):

**Supervivencia (no haber parido pretérmino hasta \(t\)):**

\[
\hat{S}_i(t \mid s) = \prod_{k=1}^{t} \hat{q}_i(k \mid s)
\]

**Riesgo acumulado de PTB hasta la semana \(t\):**

\[
\hat{R}_i(t \mid s) = 1 - \hat{S}_i(t \mid s)
\]

Interpretación: \(\hat{R}_i(t \mid s)\) es la probabilidad simulada de que la persona \(i\) haya tenido un parto pretérmino en o antes de la semana \(t\) bajo el escenario \(s\).

---

## 6. Estimaciones semanales poblacionales (`weekly_effects`)

En cada semana de seguimiento \(t\), sobre el conjunto de personas aún en riesgo en \(t\) (denominador \(n_t\)):

**Riesgo acumulado medio — curso natural:**

\[
\overline{R}_{\mathrm{nat}}(t) = \frac{1}{n_t} \sum_{i \in \mathcal{R}_t} \hat{R}_i(t \mid \mathrm{nat})
\]

**Riesgo acumulado medio — intervención:**

\[
\overline{R}_{\mathrm{int}}(t) = \frac{1}{n_t} \sum_{i \in \mathcal{R}_t} \hat{R}_i(t \mid \mathrm{int})
\]

**Contraste semanal:**

\[
\mathrm{RR}(t) = \frac{\overline{R}_{\mathrm{int}}(t)}{\overline{R}_{\mathrm{nat}}(t)}, \qquad
\mathrm{RD}(t) = \overline{R}_{\mathrm{int}}(t) - \overline{R}_{\mathrm{nat}}(t)
\]

| Columna Excel | Fórmula / contenido |
|---------------|---------------------|
| `week` | Semana de seguimiento \(t\) |
| `risk_natural` | \(\overline{R}_{\mathrm{nat}}(t)\) |
| `risk_intervention` | \(\overline{R}_{\mathrm{int}}(t)\) |
| `risk_ratio` | \(\mathrm{RR}(t)\); NA si \(\overline{R}_{\mathrm{nat}}(t) = 0\) |
| `risk_difference` | \(\mathrm{RD}(t)\) |
| `*_lcl`, `*_ucl` | Percentiles 2,5 y 97,5 del bootstrap (200 réplicas) |

**Interpretación sustantiva:** en la semana \(t\), qué fracción de la cohorte (entre quienes siguen embarazadas en \(t\)) habría acumulado un parto pretérmino hasta entonces bajo cada escenario. Valores en escala 0–1 (0,032 = 3,2 %).

**Figura 3:** graficar \(\overline{R}_{\mathrm{nat}}(t)\) y \(\overline{R}_{\mathrm{int}}(t)\) para \(t = 28, …, 36\). La hoja `figure3` y el RDS en `Other/` son el mismo contenido con columnas renombradas (`cumulative_risk_observed`, `cumulative_risk_intervention`).

**Nota:** \(\hat{R}_i(t)\) es monótono creciente en \(t\) a nivel individual, pero \(\overline{R}(t)\) puede no serlo porque cambia la composición de \(\mathcal{R}_t\) (quienes ya parieron salen del promedio).

---

## 7. Estimadores globales (`population_effects`, semana 36)

Sea \(T = 36\) la semana objetivo (`population_week`).

**Prevalencia simulada por escenario:**

\[
\hat{P}^{(s)} = \frac{1}{N} \sum_{i=1}^{N} \hat{R}_i(T \mid s)
\]

**Casos esperados:**

\[
\hat{C}^{(s)} = \hat{P}^{(s)} \times N
\]

**Contraste global** (fila `intervention` respecto de `observed`):

\[
\mathrm{RR}_{\mathrm{global}} = \frac{\hat{P}^{(\mathrm{int})}}{\hat{P}^{(\mathrm{nat})}}
\]

\[
\mathrm{RD}_{\mathrm{global}} = \hat{P}^{(\mathrm{int})} - \hat{P}^{(\mathrm{nat})}
\]

\[
\mathrm{AR} = \hat{P}^{(\mathrm{nat})} - \hat{P}^{(\mathrm{int})} = -\mathrm{RD}_{\mathrm{global}}
\]

| Columna | Fórmula | Interpretación |
|---------|---------|----------------|
| `scenario` | `observed` / `intervention` | Escenario evaluado |
| `prevalence` | \(\hat{P}^{(s)}\) | Probabilidad poblacional simulada de PTB acumulado hasta la semana 36 |
| `cases` | \(\hat{C}^{(s)}\) | Número esperado de casos de PTB en la cohorte |
| `risk_ratio` | \(\mathrm{RR}_{\mathrm{global}}\) | Cuántas veces mayor/menor es el riesgo acumulado con intervención |
| `risk_difference` | \(\mathrm{RD}_{\mathrm{global}}\) | Cambio absoluto en probabilidad acumulada (puntos proporcionales) |
| `attributable_risk` | \(\mathrm{AR}\) | Riesgo proporcional evitado por la intervención (positivo = beneficio) |

**Ejemplo:** \(\mathrm{RD}_{\mathrm{global}} = -0{,}006\) → la intervención reduciría en 0,6 puntos porcentuales la probabilidad acumulada de PTB a la semana 36. Con \(N = 7\,000\), eso equivale a unos 42 casos evitados (\(0{,}006 \times N\)).

---

## 8. Figura 4 — heatmap (`figure4_long`, `figure4_wide`)

Para cada semana de intervención \(j \in \{1, …, 44\}\) y semana de seguimiento \(t \in \{28, …, 36\}\):

1. Construir exposición con reducción del 20 % **solo** en la semana \(j\).
2. Calcular \(\overline{R}_{\mathrm{nat}}(t)\) (referencia común) y \(\overline{R}_{\mathrm{int}, j}(t)\) bajo esa intervención puntual.
3. Diferencia celda del heatmap:

\[
\mathrm{RD}_{j,t} = \overline{R}_{\mathrm{int}, j}(t) - \overline{R}_{\mathrm{nat}}(t)
\]

| Columna | Interpretación |
|---------|----------------|
| `intervention_week` | Semana \(j\) en que se interviene (eje X) |
| `follow_up_week` | Semana \(t\) de evaluación del riesgo acumulado (eje Y) |
| `risk_difference` | \(\mathrm{RD}_{j,t}\) |
| `risk_natural` | \(\overline{R}_{\mathrm{nat}}(t)\) |
| `risk_intervention` | \(\overline{R}_{\mathrm{int}, j}(t)\) |

`figure4_wide` pivota \(\mathrm{RD}_{j,t}\) con filas = `follow_up_week`, columnas = `intervention_week`.

**Distinción importante:** la Figura 4 evalúa intervenciones **marginales** (una semana a la vez). No es lo mismo que una intervención **simultánea en todas las semanas** (Figura 3 / `population_effects` con escenario `pct20` global).

---

## 9. Bootstrap e intervalos de confianza

Para cada réplica \(b = 1, …, B\) (\(B = 200\)):

1. Simular coeficientes Cox del modelo de la semana \(t\):
   \[
   \tilde{\boldsymbol{\beta}}_t^{(b)} \sim \mathcal{N}\left(\hat{\boldsymbol{\beta}}_t,\, \widehat{\mathrm{Var}}(\hat{\boldsymbol{\beta}}_t)\right)
   \]
2. Repetir pasos 4–7 con \(\tilde{\boldsymbol{\beta}}_t^{(b)}\) fijos dentro de la réplica.
3. IC 95 % para cualquier estimador \(\hat{\theta}\):
   \[
   \left[\, Q_{0{,}025}(\hat{\theta}^{(1)}, …, \hat{\theta}^{(B)}),\; Q_{0{,}975}(\hat{\theta}^{(1)}, …, \hat{\theta}^{(B)}) \,\right]
   \]

Columnas `*_lcl` y `*_ucl` en Excel y RDS corresponden a estos percentiles.

---

## 10. De estimaciones semanales a globales (resumen del flujo)

```
Exposición observada / intervenida
        ↓
Ajuste Cox por semana t ≥ 28 (solo natural, una vez)
        ↓
Predicción lp → p_i(t) → q_i(t) = 1 - p_i(t)
        ↓
S_i(t) = ∏ q_i(k)  →  R_i(t) = 1 - S_i(t)
        ↓
Promedio semanal R̄(t)  →  weekly_effects / Figura 3
        ↓
Promedio R_i(36)  →  prevalence, cases, RR, RD, AR  →  population_effects
        ↓
Loop intervención puntual j  →  RD_{j,t}  →  Figura 4
        ↓
Bootstrap (200×)  →  IC 95 %
```

---

## 11. Referencia rápida de signos

| Condición | Interpretación |
|-----------|----------------|
| \(\mathrm{RD} < 0\) | Menos PTB acumulado con intervención (efecto protector) |
| \(\mathrm{RD} > 0\) | Más PTB acumulado con intervención |
| \(\mathrm{RR} < 1\) | Riesgo relativo reducido |
| \(\mathrm{RR} > 1\) | Riesgo relativo aumentado |
| \(\mathrm{AR} > 0\) | Proporción de riesgo evitada respecto del curso natural |

---

## 12. Mapa de archivos de salida

| Ruta | Contenido |
|------|-----------|
| `Summary_results/{stub}_point_estimates.xlsx` | Hojas: `weekly_effects`, `population_effects`, `figure3`, `figure4_long`, `figure4_wide` |
| `Other/{stub}_figure3.rds` | Tabla Figura 3 (punto) |
| `Other/{stub}_figure4.rds` | Lista `long` + matriz `wide` Figura 4 |
| `WeeklyEffects/{stub}_weekly_effects.rds` | Punto + bootstrap semanal completo |
| `PopulationEffects/{stub}_population_effects.rds` | Punto + bootstrap global |
| `Interventions/{id}.rds` | Historias de exposición intervenidas (script 10.0) |
| `Timing/` | Registro de tiempos de ejecución |
