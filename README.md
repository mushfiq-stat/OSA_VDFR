# OSA analysis: AHI and sleep-bout durations

R code for a variable-domain functional regression (VDFR) analysis of the apnea-hypopnea index (AHI) and log(AHI) on the durations of sleep-stage bouts. A bout is a run of consecutive 30-second epochs in the same stage.

This code accompanies the paper:

> Rahman MM, Guha Niyogi P. Association between stage-specific sleep bout durations and obstructive sleep apnea severity: A variable-domain functional regression approach. *medRxiv* (preprint), 2026. DOI: [10.64898/2026.07.14.26358060](https://doi.org/10.64898/2026.07.14.26358060)

If you use this code, please cite the paper (see [Citation](#citation)).

The project contains two scripts:

| Script | What it does |
|---|---|
| `R/01_simulate_data.R` | Writes a simulated data set in the same format as the real one, so the analysis can be run and shared without the original data. |
| `R/02_analysis.R` | Runs the descriptive table, the participant plots, the linear models and the four VDFR analyses. It writes Excel tables and PDF figures. |

---

## Project layout

```
<project>/
├── README.md
├── R/
│   ├── 01_simulate_data.R
│   └── 02_analysis.R
├── data/                         written by 01_simulate_data.R
│   ├── demographics.csv
│   ├── data_5stage.rda
│   ├── data_3stage.rda
│   ├── data_2stage.rda
│   ├── SIMULATED_DATA.txt
│   └── epoch/
│       ├── ID0001.csv
│       ├── ID0002.csv
│       └── ...
└── result/                       written by 02_analysis.R
    ├── tables/
    └── plot/
        ├── LM_diagnostics/
        ├── 5stage/
        ├── 3stage/
        └── 2stage/
```

Both scripts locate the project folder by themselves. You can run them from the project folder or from `R/`.

---

## Requirements

- R 4.1 or later
- Packages: `pacman`, `tidyverse`, `refund`, `ggfortify`, `openxlsx`

`pacman::p_load()` at the top of each script installs any package that is missing. To install them yourself instead:

```r
install.packages(c("pacman", "tidyverse", "refund", "ggfortify", "openxlsx"))
```

---

## How to run

```r
setwd("path/to/<project>")
source("R/01_simulate_data.R")   # about 30 seconds: writes data/
source("R/02_analysis.R")        # about 1 minute:   writes result/
```

Each script starts with `rm(list = ls())`, so it runs in a clean workspace.

---

## 1. `R/01_simulate_data.R` — simulated data

### What it writes

| File | Contents |
|---|---|
| `data/epoch/<ID>.csv` | One file per participant (`ID0001`, `ID0002`, …) covering the whole recording, one row per bout. |
| `data/demographics.csv` | One row per participant: the model covariates and the outcome. |
| `data/data_5stage.rda` | `final_5stage`: bout-duration matrices for W, N1, N2, N3 and REM. |
| `data/data_3stage.rda` | `final_3stage`: bout-duration matrices for W, N (N1–N3 merged) and REM. |
| `data/data_2stage.rda` | `final_2stage`: bout-duration matrices for wake and non-wake (all sleep stages merged). |
| `data/SIMULATED_DATA.txt` | Marker file holding the scenario, seed, sample size and creation time. |

#### `data/epoch/<ID>.csv`

| Column | Meaning |
|---|---|
| `stage` | `W`, `N1`, `N2`, `N3` or `REM` |
| `start_time` | Clock time the bout starts (`hh:mm:ss`) |
| `end_time` | Clock time the bout ends (`hh:mm:ss`) |
| `time` | Bout length in seconds |

The rows are in the order of the night. Recordings start between 22:00 and 24:00, so the clock times can cross midnight.

#### `data/demographics.csv`

| Column | Meaning |
|---|---|
| `ID` | Participant ID; matches the epoch file name |
| `gender` | Female / Male |
| `age` | Years |
| `bmiblquan` | BMI |
| `ethnicity` | Other / White / Black / Asian |
| `maritalstatushp` | Other / Single / Married |
| `highestgradehp` | Years of education (about 2% missing) |
| `ahi` | Apnea-hypopnea index (the outcome) |

#### Bout-duration matrices (`.rda` files)

Each object is a participant-level tibble holding `ID`, the covariates, `ahi` and one matrix column per stage. In each matrix:

- a row is a participant and a column is the bout order;
- the values are bout durations in **minutes**;
- each row is filled from the left and padded with `NA` after the participant's last bout.

The matrices are built the same way as the real ones. **The first bout is dropped if it is W, and so is the last bout**, which removes the wake before sleep onset and after final awakening. For the 3-stage and 2-stage groupings, consecutive bouts that fall in the same merged stage are joined into one bout.

| Object | Elements / matrices | Who is included |
|---|---|---|
| `final_5stage` | a list with elements `W`, `N1`, `N2`, `N3`, `REM`, each holding all five matrices | element `s` = participants with at least one bout of stage `s`; `$W` = everyone |
| `final_3stage` | matrices `N`, `REM`, `W` | participants with all three stages |
| `final_2stage` | matrices `wake`, `non_wake` | participants with both |

### How the data are generated

1. **Covariates** are drawn to match the real sample: 65.5% male, age 51.6 (SD 12.2), BMI 32.2 (7.1), years of education 15.5 (2.6), and the observed marital-status and ethnicity distributions.
2. **One night per participant** comes from a semi-Markov model on 30-second epochs:
   - wake before sleep onset (about 15 min), then alternating sleep-stage bouts, then a final wake (about 10 min); recordings last 7.5–9 h;
   - the next stage is drawn from transition weights that change during the night, with N3 more likely early and REM more likely late;
   - bout lengths are lognormal;
   - about half the participants have no N3, as in the real data;
   - a participant-level fragmentation score shortens every bout.
3. **The outcome** is generated as log(AHI) = α + Z′γ + sleep term + error, and AHI = exp(log(AHI)). The covariate effects are close to those estimated on the real data, and the mean and SD of log(AHI) match the real data (3.50 and 0.64).

### Settings (section S0)

| Setting | Default | Meaning |
|---|---|---|
| `scenario` | `"fragmentation"` | `"fragmentation"`: more fragmented sleepers have higher AHI, so the bout-duration models have a real signal to find. `"null"`: AHI is unrelated to sleep. |
| `n_sim` | `1103` | Number of participants |
| `seed` | `2026` | Random seed |
| `signal_share` | `0.10` | Share of the variance of log(AHI) explained by the sleep term |
| `bout_median`, `bout_sdlog` | per stage | Lognormal bout-length parameters, in epochs |
| `exit_w`, `n3_time`, `rem_time` | — | Transition weights and how they change through the night |

At the end of a run the script prints a check against the real data: sample sizes, median number of bouts per stage, and the mean and SD of AHI and log(AHI).

### Safety

The script never writes into a `data/` folder that it did not create. A folder counts as its own only if it contains `SIMULATED_DATA.txt`. On a re-run it replaces only its own earlier files.

---

## 2. `R/02_analysis.R` — analysis

### Model

All models adjust for gender, age, BMI, ethnicity, marital status and years of education. Every analysis is run twice, once for **AHI** and once for **log(AHI)**.

The VDFR models are fitted with `refund::pfr()`. Each stage enters as a term `lf.vd(X, transform = "standardized")`:

```
Y_i = alpha + Z_i' gamma + sum over stages of  integral beta(t, T_i) X_i(t) dt  + e_i
```

- X_i(t) is participant i's bout-duration curve for that stage.
- T_i is the participant's number of bouts of that stage.
- β(t, T_i) is a coefficient surface that may change with T_i.

### Sections and parts

| Part | Analysis | Data |
|---|---|---|
| Section 3 | Hypnogram and bout-duration plot for one participant (`plot_id`) | `data/epoch/<ID>.csv` |
| Section 4 | Table 1: descriptive statistics | `final_5stage$W` |
| **A** | Linear regression on **mean** bout duration: one model per stage plus one with all five means, with diagnostic plots | `final_5stage$W` |
| **B** | **5 stages, separate**: one VDFR per stage (W, N1, N2, N3, REM) | `final_5stage[[stage]]` |
| **C** | **5 stages, joint**: W + N1 + N2 + N3 + REM in one VDFR | participants with all five stages |
| **D** | **3 stages, separate**: one VDFR per stage (W, N, REM) | `final_3stage` |
| **E** | **2 stages, separate**: one VDFR per stage (wake, non-wake) | `final_2stage` |

### Outputs

#### Tables: `result/tables/`

| File | Sheets |
|---|---|
| `Table1_demographics.xlsx` | Table 1 |
| `LM_mean_bout_duration.xlsx` | `AHI`, `log(AHI)`, `estimates` |
| `5stage_separate_VDFR.xlsx` | `AHI`, `log(AHI)`, `estimates` |
| `5stage_joint_VDFR.xlsx` | `AHI`, `log(AHI)`, `estimates`, `sample_size` |
| `3stage_separate_VDFR.xlsx` | `AHI`, `log(AHI)`, `estimates` |
| `2stage_separate_VDFR.xlsx` | `AHI`, `log(AHI)`, `estimates` |

- **`AHI` and `log(AHI)` sheets** are formatted tables. Rows are the model terms and columns are the stages (or the models). Each cell reads `estimate (95% CI) stars`, where * means p < 0.05, ** p < 0.01 and *** p < 0.001. The VDFR tables also have **EDF** rows: the effective degrees of freedom of each coefficient surface, starred by the p-value for the test of no effect of that stage. Each table ends with R², deviance explained (VDFR) or adjusted R² (linear models), and n.
- **`estimates` sheet** holds the same results as unformatted numbers (estimate, SE, t, p).
- **`sample_size` sheet** (joint model only) shows how many participants have each stage and how many have all five.

The 95% CI is a Wald interval (estimate ± 1.96 SE) for the VDFR models and comes from `confint()` (t-based) for the linear models.

#### Plots: `result/plot/`

| Location | Files |
|---|---|
| `result/plot/` | `Hypnogram <ID>.pdf`, `Duration <ID>.pdf` |
| `result/plot/LM_diagnostics/` | `diag_<AHI or logAHI>_<W, N1, N2, N3, REM or Combined>.pdf`: residuals vs fitted, normal Q-Q and scale-location |
| `result/plot/5stage/` | `5stage_separate_<stage>[_logAHI].pdf`, `5stage_joint_<stage>[_logAHI].pdf` |
| `result/plot/3stage/` | `3stage_separate_<stage>[_logAHI].pdf` |
| `result/plot/2stage/` | `2stage_separate_<stage>[_logAHI].pdf` |

The coefficient-surface heat maps show β(t, T_i) on the triangle t ≤ T_i, where it is defined. The axes are on the original bout scale: **t is the bout number** (occurrence within the stage) and **T_i is the participant's number of bouts** of that stage. File names without `_logAHI` are for AHI.

### Settings (Section 0)

| Setting | Default | Meaning |
|---|---|---|
| `run_participant_plot`, `run_table1` | `TRUE` | Switch Sections 3 and 4 on or off |
| `run_part_A` … `run_part_E` | `TRUE` | Switch each part on or off |
| `plot_id` | `"ID0063"` | Participant shown in Section 3 |
| `covars` | six covariates | Covariates in every model |
| `outcomes` | AHI, log(AHI) | Outcomes, their table labels and file-name suffixes |
| `complete_only` (Part A) | `FALSE` | `TRUE` restricts every linear model to participants with all five stages |

---

## Using the real data

The paper uses baseline polysomnography from the Apnea Positive Pressure Long-term Efficacy Study (APPLES, n = 1,103). The de-identified data are available from the National Sleep Research Resource (NSRR) at <https://sleepdata.org> after a data-access request; they are not included in this project.

`R/02_analysis.R` runs on real data without changes if `data/` holds files in the same format:

- `data_5stage.rda`, `data_3stage.rda` and `data_2stage.rda`, with the object names and matrix columns described above;
- `epoch/<ID>.csv` with columns `stage`, `start_time`, `end_time` and `time`.

`demographics.csv` is not read by the analysis. The covariates and AHI come from the `.rda` files.

## Notes

- Sleep is scored in 30-second epochs, so every bout length is a multiple of 30 seconds.
- In the simulated data every recording starts and ends with W. Dropping the first and last W bout therefore removes exactly the wake before sleep onset and after final awakening.
- Participants missing years of education are dropped from the models (about 2%), so a model's n can be smaller than its data set.

---

## Citation

If you use this code, please cite:

> Rahman MM, Guha Niyogi P. Association between stage-specific sleep bout durations and obstructive sleep apnea severity: A variable-domain functional regression approach. *medRxiv* (preprint), 2026.
> DOI: [10.64898/2026.07.14.26358060](https://doi.org/10.64898/2026.07.14.26358060)

BibTeX:

```bibtex
@article{rahman2026vdfr,
  title   = {Association between stage-specific sleep bout durations and obstructive sleep apnea severity: A variable-domain functional regression approach},
  author  = {Rahman, Md Mushfiqur and Guha Niyogi, Pratim},
  journal = {medRxiv},
  year    = {2026},
  note    = {Preprint},
  doi     = {10.64898/2026.07.14.26358060},
  url     = {https://doi.org/10.64898/2026.07.14.26358060}
}
```

The VDFR method implemented by `refund::lf.vd()` is described in:

> Gellar JE, Colantuoni E, Needham DM, Crainiceanu CM. Variable-domain functional regression for modeling ICU data. *Journal of the American Statistical Association*. 2014;109(508):1425–1439. DOI: [10.1080/01621459.2014.940044](https://doi.org/10.1080/01621459.2014.940044)

Users of the APPLES data should also follow the NSRR citation and acknowledgement requirements listed at <https://sleepdata.org>.
