################################################################################
#  02_analysis.R
#  OSA ANALYSIS: AHI and log(AHI) vs. sleep-bout durations
#  ---------------------------------------------------------------------------
#  Outcomes   : AHI (apnea-hypopnea index) and log(AHI)
#  Exposures  : per-participant sequences of sleep-bout durations by stage.
#               Each participant has a different number of bouts, so each
#               stage is a curve on its own domain [0, T_i] -> variable-domain
#               functional regression (VDFR), refund::pfr() + lf.vd().
#  Covariates : gender, age, BMI, ethnicity, marital status, years of education
#
#  Section 3 : one participant's hypnogram + bout-duration plot
#  Section 4 : descriptive table (Table 1)
#
#  Models (each for AHI and for log(AHI)):
#    PART A : linear regression on mean bout duration (5 stages): one model
#             per stage + one with all five means, with diagnostic plots
#    PART B : 5 stages, separate  - one VDFR per stage: W, N1, N2, N3, REM
#    PART C : 5 stages, joint     - W + N1 + N2 + N3 + REM in one VDFR
#                                   (participants with all five stages)
#    PART D : 3 stages, separate  - one VDFR per stage: W, N (N1-N3), REM
#    PART E : 2 stages, separate  - one VDFR per stage: wake, non-wake
#
#  Project layout (run from the project folder, or from R/):
#      <project>/
#        R/              01_simulate_data.R, 02_analysis.R
#        data/           data_5stage.rda, data_3stage.rda, data_2stage.rda,
#                        demographics.csv
#          epoch/        <ID>.csv: stage, start_time, end_time, time (s)
#        result/
#          tables/       Excel files with the coefficient tables (written here)
#          plot/         participant plots (written here), and
#            LM_diagnostics/  residual plots of the PART A linear models
#            5stage/     coefficient surfaces, PARTS B and C
#            3stage/     coefficient surfaces, PART D
#            2stage/     coefficient surfaces, PART E
#
#  Run order: 0 Setup -> 1 Helpers -> 2 Data -> 3 Participant plot ->
#             4 Table 1 -> PART A -> B -> C -> D -> E.
#  Sections 3-4 and each PART can be switched off in Section 0.
################################################################################

rm(list = ls())


################################################################################
# 0. SETUP: packages, paths, switches, model settings
################################################################################

pacman::p_load(
  tidyverse,   # dplyr / tidyr / purrr / ggplot2 / readr / stringr
  refund,      # pfr() and lf.vd() for variable-domain functional regression
  ggfortify,   # autoplot() for the linear-model diagnostic plots
  openxlsx     # Excel output
)

# ---- Paths -------------------------------------------------------------------
# The project folder holds R/, data/ and result/. Run from it or from R/.
project_dir <- getwd()
if (basename(project_dir) == "R") project_dir <- dirname(project_dir)
data_dir  <- file.path(project_dir, "data")
epoch_dir <- file.path(data_dir, "epoch")
table_dir <- file.path(project_dir, "result", "tables")
plot_dir  <- file.path(project_dir, "result", "plot")

if (!file.exists(file.path(data_dir, "data_5stage.rda"))) {
  stop("data/data_5stage.rda not found under ", project_dir,
       ". Run from the project folder (or R/), after R/01_simulate_data.R.", call. = FALSE)
}
plot_sub <- c("LM_diagnostics", "5stage", "3stage", "2stage")
for (d in c(table_dir, plot_dir, file.path(plot_dir, plot_sub))) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ---- Which parts to run ------------------------------------------------------
run_participant_plot <- TRUE   # Section 3: one participant's hypnogram + duration plot
run_table1           <- TRUE   # Section 4: descriptive table
run_part_A           <- TRUE   # linear regression on mean bout duration + diagnostics
run_part_B           <- TRUE   # 5 stages, separate (VDFR)
run_part_C           <- TRUE   # 5 stages, joint (VDFR)
run_part_D           <- TRUE   # 3 stages, separate (VDFR)
run_part_E           <- TRUE   # 2 stages, separate (VDFR)

# Participant shown in Section 3
plot_id <- "ID0063"

# ---- Model settings shared by every part ---------------------------------------
covars <- c("gender", "age", "bmiblquan", "ethnicity",
            "maritalstatushp", "highestgradehp")

# Outcomes: label used in tables, column in the data, suffix for file names
outcomes <- tibble(
  label  = c("AHI",  "log(AHI)"),
  column = c("ahi",  "log_ahi"),
  suffix = c("",     "_logAHI")
)

# Colours for the five stages (Section 3)
stage_cols <- c(W = "#F26419", N1 = "#8ECAE6", N2 = "#219EBC",
                N3 = "#023047", REM = "#F6AE2D")

# PDF device: cairo_pdf renders the plot-math axis labels (T[i]) best
pdf_device <- if (capabilities("cairo")) cairo_pdf else "pdf"


################################################################################
# 1. HELPER FUNCTIONS
################################################################################

# ---- 1.1 Add log(AHI), refusing silently-bad values ----------------------------
add_log_ahi <- function(df) {
  n_bad <- sum(df$ahi <= 0, na.rm = TRUE)
  if (n_bad > 0) warning(n_bad, " participant(s) have AHI <= 0; log(AHI) is undefined for them.")
  df %>% mutate(log_ahi = log(ahi))
}

# ---- 1.2 Fit a VDFR model ------------------------------------------------------
# One lf.vd() term per stage, plus the scalar covariates.
#   lf.vd(X, transform = "standardized"): each curve is rescaled to a common
#   [0, 1] domain, and the coefficient surface beta(t, T_i) may vary with the
#   participant's number of bouts T_i.
fit_vdfr <- function(data, y, stages) {
  vd_terms <- sprintf("lf.vd(%s, transform = 'standardized')", stages)
  f <- as.formula(paste(y, "~", paste(c(vd_terms, covars), collapse = " + ")))
  pfr(f, family = gaussian(), data = data)
}

# ---- 1.3 Tidy a VDFR fit into one table --------------------------------------
# Parametric rows : estimate, SE, t, p for intercept and covariates
# Functional rows : EDF of each coefficient surface and its p-value (test of
#                   "no effect of this stage"), plus reference df
# Fit rows        : adjusted R-squared, deviance explained, n
tidy_pfr <- function(fit, stages) {
  s  <- summary(fit)
  st <- s$s.table
  k  <- length(s$p.coeff)
  edf_names <- if (length(stages) == 1) "edf"     else paste0("edf_", stages)
  ref_names <- if (length(stages) == 1) "ref_edf" else paste0("ref_edf_", stages)
  
  bind_rows(
    tibble(term = names(s$p.coeff), estimate = unname(s$p.coeff),
           se = unname(s$se[seq_len(k)]), t = unname(s$p.t), pval = unname(s$p.pv)),
    tibble(term = c(rbind(edf_names, ref_names)),
           estimate = c(rbind(st[, "edf"], st[, "Ref.df"])),
           se = NA_real_, t = NA_real_, pval = c(rbind(st[, "p-value"], NA_real_))),
    tibble(term = c("r_square", "dev_explained", "n"),
           estimate = c(s$r.sq, s$dev.expl, s$n),
           se = NA_real_, t = NA_real_, pval = NA_real_)
  )
}

# ---- 1.4 Coefficient table -> readable wide table ----------------------------
# One wide table per response: rows = terms, columns = stages (or "Joint"),
# cells = "estimate (95% CI) stars" (Wald CI, estimate +/- 1.96 SE);
# EDF cells = "edf stars" (p-value of the functional term).
term_label <- function(term) {
  case_when(
    term == "edf"                        ~ "EDF (functional term)",
    term == "ref_edf"                    ~ "Reference df",
    str_detect(term, "^edf_")            ~ str_replace(term, "^edf_(.*)$", "EDF (\\1)"),
    str_detect(term, "^ref_edf_")        ~ str_replace(term, "^ref_edf_(.*)$", "Reference df (\\1)"),
    term == "r_square"                   ~ "R-squared (adj.)",
    term == "adj_r_squared"              ~ "R-squared (adj.)",
    term == "r_squared"                  ~ "R-squared",
    str_detect(term, "^avg_")            ~ paste("Mean bout duration:",
                                                 toupper(str_remove(term, "^avg_"))),
    term == "dev_explained"              ~ "Deviance explained",
    term == "n"                          ~ "n",
    term == "age"                        ~ "Age",
    term == "bmiblquan"                  ~ "BMI",
    term == "highestgradehp"             ~ "Years of education",
    str_detect(term, "^gender")          ~ str_replace(term, "^gender", "Gender: "),
    str_detect(term, "^ethnicity")       ~ str_replace(term, "^ethnicity", "Ethnicity: "),
    str_detect(term, "^maritalstatushp") ~ str_replace(term, "^maritalstatushp", "Marital status: "),
    term == "(Intercept)"                ~ "Intercept",
    TRUE                                 ~ term
  )
}

format_coef <- function(df) {
  # Wald CI unless the table already has one (lower/upper, e.g. from confint())
  if (!all(c("lower", "upper") %in% names(df))) {
    df <- df %>% mutate(lower = estimate - 1.96 * se, upper = estimate + 1.96 * se)
  }
  df %>%
    mutate(
      stars = case_when(is.na(pval) ~ "", pval < 0.001 ~ "***",
                        pval < 0.01 ~ "**", pval < 0.05 ~ "*", TRUE ~ ""),
      value = case_when(
        term == "n"               ~ sprintf("%.0f", estimate),
        str_detect(term, "^edf")  ~ trimws(sprintf("%.3f %s", estimate, stars)),
        is.na(se)                 ~ sprintf("%.3f", estimate),
        TRUE ~ trimws(sprintf("%.3f (%.3f, %.3f) %s", estimate, lower, upper, stars))
      ),
      Term = term_label(term)
    ) %>%
    select(response, stage, Term, value) %>%
    pivot_wider(names_from = stage, values_from = value) %>%
    # row order: intercept, exposures, covariates, functional terms, fit statistics
    mutate(.order = case_when(
      Term == "Intercept"                            ~ 1,
      str_detect(Term, "^Mean bout duration")        ~ 2,
      str_detect(Term, "^EDF|^Reference df")         ~ 4,
      str_detect(Term, "^R-squared|^Deviance|^n$")   ~ 5,
      TRUE                                           ~ 3
    )) %>%
    arrange(response, .order) %>%
    select(-.order) %>%
    split(.$response) %>%
    map(~ select(.x, -response))
}

# ---- 1.5 Write a named list of tables to one Excel file (one sheet each) ------
write_tables <- function(sheets, file, note = NULL) {
  wb <- createWorkbook()
  hs <- createStyle(textDecoration = "bold", fgFill = "#E8EEF4", border = "bottom")
  for (nm in names(sheets)) {
    sh <- substr(gsub("[\\[\\]\\*\\?/\\\\:]", "", nm), 1, 31)   # Excel sheet-name rules
    addWorksheet(wb, sh)
    start <- 1
    if (!is.null(note)) { writeData(wb, sh, note, startRow = 1); start <- 3 }
    writeData(wb, sh, sheets[[nm]], startRow = start, headerStyle = hs)
    setColWidths(wb, sh, cols = seq_len(ncol(sheets[[nm]])), widths = "auto")
    freezePane(wb, sh, firstActiveRow = start + 1)
  }
  saveWorkbook(wb, file, overwrite = TRUE)
  cat("  Table:", file.path("result/tables", basename(file)), "\n")
}

stars_note <- paste("Cells: estimate (95% CI); EDF rows: effective df of the coefficient",
                    "surface. * p < 0.05, ** p < 0.01, *** p < 0.001.")

# ---- 1.6 Coefficient surface heat map ------------------------------------------
# beta(t, T_i) on the triangle t <= T_i. With n_bouts = J (columns of the
# stage matrix) both axes are mapped back from lf.vd()'s [0, 1] scale to
# bout counts: t = bout number (occurrence), T_i = number of bouts.
plot_heatmap <- function(fit, select = 1, n_bouts = NULL, title = NULL,
                         n1 = 100, n2 = 100) {
  coef_df <- coef(fit, select = select, n1 = n1, n2 = n2)
  arg_col <- grep("arg", names(coef_df), value = TRUE)[1]
  vd_col  <- grep("vd",  names(coef_df), value = TRUE)[1]
  coef_df <- coef_df[coef_df[[arg_col]] <= coef_df[[vd_col]], ]
  if (!is.null(n_bouts)) {
    coef_df[[arg_col]] <- 1 + coef_df[[arg_col]] * (n_bouts - 1)
    coef_df[[vd_col]]  <- 1 + coef_df[[vd_col]]  * (n_bouts - 1)
  }
  rng  <- range(coef_df$value, na.rm = TRUE)
  brks <- round(seq(rng[1], rng[2], length.out = 6)[2:5], 2)
  
  ggplot(coef_df, aes(x = .data[[arg_col]], y = .data[[vd_col]])) +
    geom_tile(aes(fill = value)) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_gradient2(low = "blue", mid = "yellow", high = "red",
                         midpoint = 0, breaks = brks) +
    labs(x = expression(t), y = expression(T[i]), title = title) +
    theme_bw() +
    theme(panel.grid   = element_blank(),
          plot.title   = element_text(hjust = 0.5, face = "bold"),
          axis.text    = element_text(size = 13, face = "bold"),
          axis.title   = element_text(size = 15, face = "bold"),
          legend.text  = element_text(size = 10, face = "bold"),
          legend.title = element_blank())
}

# `sub` = subfolder of result/plot ("" for result/plot itself)
save_plot <- function(p, name, sub = "", width = 4, height = 3.5) {
  ggsave(file.path(plot_dir, sub, name), p, width = width, height = height,
         units = "in", device = pdf_device)
}

# ---- 1.7 Separate models: one VDFR per stage, both outcomes ------------------------
# data_for(s) returns the analysis data for stage s. Writes one Excel file
# (sheets AHI, log(AHI), estimates) and one heat map per stage and outcome
# into result/plot/<sub>/.
run_separate <- function(stages, data_for, part_name, file_stub, sub) {
  coef_tab <- tibble()
  for (o in seq_len(nrow(outcomes))) {
    cat(sprintf("  %s, %s\n", part_name, outcomes$label[o]))
    for (s in stages) {
      dat <- data_for(s)
      fit <- fit_vdfr(dat, outcomes$column[o], s)
      coef_tab <- bind_rows(coef_tab,
                            tidy_pfr(fit, s) %>%
                              mutate(response = outcomes$label[o], stage = s, .before = 1))
      save_plot(plot_heatmap(fit, n_bouts = ncol(dat[[s]]),
                             title = paste0(outcomes$label[o], ": ", s)),
                paste0(file_stub, "_", s, outcomes$suffix[o], ".pdf"), sub = sub)
    }
  }
  write_tables(c(format_coef(coef_tab), list(estimates = coef_tab)),
               file.path(table_dir, paste0(file_stub, "_VDFR.xlsx")),
               note = paste("Columns: stage (one model each).", stars_note))
  invisible(coef_tab)
}


################################################################################
# 2. LOAD DATA
################################################################################
# final_5stage : list (W, N1, N2, N3, REM); element s = participants with at
#                least one bout of stage s; $W = all participants
# final_3stage : matrix columns N, REM, W
# final_2stage : matrix columns wake, non_wake
# Each is a participant-level tibble: ID, covariates, ahi and the bout-duration
# matrices (row = participant, column = bout order, minutes; NA-padded).

cat("Loading data from", data_dir, "\n")
load(file.path(data_dir, "data_5stage.rda"))   # -> final_5stage
load(file.path(data_dir, "data_3stage.rda"))   # -> final_3stage
load(file.path(data_dir, "data_2stage.rda"))   # -> final_2stage

final_5stage <- map(final_5stage, add_log_ahi)
final_3stage <- add_log_ahi(final_3stage)
final_2stage <- add_log_ahi(final_2stage)


################################################################################
# 3. ONE PARTICIPANT'S HYPNOGRAM AND BOUT-DURATION PLOT
################################################################################
# From data/epoch/<plot_id>.csv (whole recording, one row per bout, in order).
# Outputs: result/plot/Hypnogram <ID>.pdf, result/plot/Duration <ID>.pdf

if (run_participant_plot) {
  cat("\nSection 3: participant plots for", plot_id, "...\n")
  f <- file.path(epoch_dir, paste0(plot_id, ".csv"))
  if (!file.exists(f)) stop("No file for ", plot_id, " in ", epoch_dir, call. = FALSE)
  
  bouts <- suppressMessages(read_csv(f, show_col_types = FALSE)) %>%
    mutate(stage    = factor(recode(stage, "R" = "REM"), levels = names(stage_cols)),
           end_h    = cumsum(time) / 3600,                 # hours since start of recording
           start_h  = end_h - time / 3600,
           duration = time / 60) %>%                       # minutes
    group_by(stage) %>%
    mutate(occurrence = row_number()) %>%                  # file order = night order
    ungroup()
  
  participant_theme <- theme_bw() +
    theme(panel.grid   = element_blank(),
          plot.title   = element_text(hjust = 0.5, face = "bold", size = 14),
          axis.text    = element_text(face = "bold", size = 11, colour = "black"),
          axis.title   = element_text(face = "bold", size = 12),
          legend.title = element_text(face = "bold"),
          legend.text  = element_text(face = "bold"))
  
  p_duration <- ggplot(bouts, aes(occurrence, duration, colour = stage)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    scale_colour_manual(values = stage_cols) +
    labs(title = "Duration (min) Across Stage Occurrences",
         x = "Occurrence Number within Stage", y = NULL) +
    participant_theme +
    theme(legend.position = "none")
  
  p_timeline <- ggplot(bouts) +
    geom_rect(aes(xmin = start_h, xmax = end_h,
                  ymin = as.numeric(stage) - 0.4, ymax = as.numeric(stage) + 0.4,
                  fill = stage)) +
    scale_fill_manual(values = stage_cols, name = "Stage", drop = FALSE) +
    scale_y_reverse(breaks = seq_along(stage_cols), labels = names(stage_cols),
                    limits = c(length(stage_cols) + 0.5, 0.5)) +
    scale_x_continuous(breaks = 0:ceiling(max(bouts$end_h))) +
    labs(x = "Hours Since Start of Recording", y = NULL) +
    participant_theme
  
  save_plot(p_duration, paste0("Duration ", plot_id, ".pdf"), width = 8, height = 4)
  save_plot(p_timeline, paste0("Hypnogram ", plot_id, ".pdf"), width = 7, height = 3.5)
}


################################################################################
# 4. DESCRIPTIVE TABLE (Table 1)
################################################################################
# All participants (final_5stage$W). Continuous: mean (SD). Categorical: n (%)
# of non-missing values. Output: result/tables/Table1_demographics.xlsx

if (run_table1) {
  cat("\nSection 4: descriptive table...\n")
  desc_data <- final_5stage$W
  
  fmt_mean_sd <- function(x) sprintf("%.2f (%.2f)", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
  fmt_n_pct <- function(x, levels) {
    x <- as.character(x)
    n <- sapply(levels, function(l) sum(x == l, na.rm = TRUE))
    sprintf("%d (%.2f%%)", n, 100 * n / sum(!is.na(x)))
  }
  cont_row <- function(label, x) tibble(Characteristic = label, Value = fmt_mean_sd(x))
  cat_rows <- function(label, x, levels) {
    tibble(Characteristic = c(label, paste0("  ", levels)),
           Value          = c("", fmt_n_pct(x, levels)))
  }
  
  table1 <- bind_rows(
    cat_rows("Sex",                desc_data$gender,          c("Male", "Female")),
    cont_row("Age",                desc_data$age),
    cont_row("BMI",                desc_data$bmiblquan),
    cont_row("Years of Education", desc_data$highestgradehp),
    cat_rows("Marital Status",     desc_data$maritalstatushp, c("Single", "Married", "Other")),
    cat_rows("Ethnicity",          desc_data$ethnicity,       c("Asian", "Black", "White", "Other")),
    cont_row("AHI",                desc_data$ahi),
    cont_row("log(AHI)",           desc_data$log_ahi)
  )
  names(table1)[2] <- paste0("n = ", format(nrow(desc_data), big.mark = ","))
  print(table1, n = Inf)
  write_tables(list(`Table 1` = table1), file.path(table_dir, "Table1_demographics.xlsx"),
               note = "Continuous: mean (SD). Categorical: n (%).")
}


################################################################################
# PART A: LINEAR REGRESSION ON MEAN BOUT DURATION (5 stages)
################################################################################
# Scalar benchmark for the VDFR models: each stage's curve is summarised by
# its mean bout duration (minutes).
#   - one model per stage (W, N1, N2, N3, REM)
#   - one combined model with all five means (complete cases)
# 95% CI from confint() (t-based).
# Outputs: result/tables/LM_mean_bout_duration.xlsx
#          result/plot/LM_diagnostics/diag_<outcome>_<model>.pdf
#          (residuals vs fitted, normal Q-Q, scale-location)

if (run_part_A) {
  cat("\nPART A: linear regression on mean bout duration...\n")
  
  df_lm <- final_5stage$W %>%
    mutate(
      avg_w   = rowMeans(unclass(W),   na.rm = TRUE),
      avg_n1  = rowMeans(unclass(N1),  na.rm = TRUE),
      avg_n2  = rowMeans(unclass(N2),  na.rm = TRUE),
      avg_n3  = rowMeans(unclass(N3),  na.rm = TRUE),
      avg_rem = rowMeans(unclass(REM), na.rm = TRUE)
    ) %>%
    # a stage with no bouts gives NaN; turn it into NA
    mutate(across(starts_with("avg_"), ~ ifelse(is.nan(.x), NA_real_, .x))) %>%
    select(-c(W, N1, N2, N3, REM))
  
  # TRUE = everyone must have all five stages (same n in every model);
  # FALSE = each single-stage model uses everyone who has that stage
  # (the combined model is complete-case either way, lm() drops NA)
  complete_only <- FALSE
  if (complete_only) df_lm <- df_lm %>% filter(if_all(starts_with("avg_"), ~ !is.na(.x)))
  
  lm_specs <- list(W = "avg_w", N1 = "avg_n1", N2 = "avg_n2", N3 = "avg_n3", REM = "avg_rem",
                   Combined = c("avg_w", "avg_n1", "avg_n2", "avg_n3", "avg_rem"))
  
  # lm_fits$AHI$W, lm_fits$`log(AHI)`$Combined, ...
  lm_fits <- map(set_names(outcomes$column, outcomes$label), function(y) {
    map(lm_specs, function(preds) {
      lm(as.formula(paste(y, "~", paste(c(preds, covars), collapse = " + "))), data = df_lm)
    })
  })
  
  tidy_lm <- function(fit) {
    s  <- summary(fit)
    cm <- s$coefficients
    ci <- confint(fit)[rownames(cm), , drop = FALSE]
    bind_rows(
      tibble(term = rownames(cm), estimate = cm[, "Estimate"], se = cm[, "Std. Error"],
             t = cm[, "t value"], pval = cm[, "Pr(>|t|)"], lower = ci[, 1], upper = ci[, 2]),
      tibble(term = c("r_squared", "adj_r_squared", "n"),
             estimate = c(s$r.squared, s$adj.r.squared, nobs(fit)))
    )
  }
  lm_results <- map_dfr(lm_fits, ~ map_dfr(.x, tidy_lm, .id = "stage"), .id = "response")
  
  write_tables(c(format_coef(lm_results), list(estimates = lm_results)),
               file.path(table_dir, "LM_mean_bout_duration.xlsx"),
               note = paste("Columns: model (one stage, or all five means combined).",
                            "Cells: estimate (95% CI). * p < 0.05, ** p < 0.01, *** p < 0.001."))
  
  # Diagnostic plots, one file per outcome and model
  for (o in seq_len(nrow(outcomes))) {
    out_name <- gsub("[()]", "", outcomes$label[o])           # "AHI" / "logAHI"
    for (st in names(lm_fits[[o]])) {
      p <- autoplot(lm_fits[[o]][[st]], which = 1:3, ncol = 3) + theme_bw()
      pdf(file.path(plot_dir, "LM_diagnostics", sprintf("diag_%s_%s.pdf", out_name, st)),
          width = 6, height = 2.5)
      print(p)
      dev.off()
    }
  }
}


################################################################################
# PART B: 5 STAGES, SEPARATE MODELS (W, N1, N2, N3, REM)
################################################################################
#   y_i = covariates_i'gamma + integral beta(t, T_i) X_i(t) dt + e_i, one stage
# Data for stage s: final_5stage[[s]] (participants with at least one bout of s)
# Outputs: result/tables/5stage_separate_VDFR.xlsx
#          result/plot/5stage/5stage_separate_<stage><suffix>.pdf

if (run_part_B) {
  cat("\nPART B: 5 stages, separate...\n")
  run_separate(c("W", "N1", "N2", "N3", "REM"), function(s) final_5stage[[s]],
               "5 stages, separate", "5stage_separate", sub = "5stage")
}


################################################################################
# PART C: 5 STAGES, JOINT MODEL (W + N1 + N2 + N3 + REM)
################################################################################
# All five curves in one model, so only participants with at least one bout
# of every stage (complete cases).
# Outputs: result/tables/5stage_joint_VDFR.xlsx (+ sheet sample_size)
#          result/plot/5stage/5stage_joint_<stage><suffix>.pdf

if (run_part_C) {
  cat("\nPART C: 5 stages, joint...\n")
  five <- c("W", "N1", "N2", "N3", "REM")
  all_data <- final_5stage$W
  N <- nrow(all_data)
  
  # a participant is "missing" a stage if that stage's row is entirely NA
  missing_ids <- map(set_names(five), function(s) {
    all_data$ID[apply(all_data[[s]], 1, function(x) all(is.na(x)))]
  })
  dat_joint <- all_data %>% filter(!ID %in% unique(unlist(missing_ids)))
  
  sample_size <- tibble(
    Group = c("All participants", paste("With stage", five), "All five stages (analysed)"),
    n     = c(N, map_int(missing_ids, ~ N - length(.x)), nrow(dat_joint))
  )
  print(sample_size)
  
  joint_tab <- tibble()
  for (o in seq_len(nrow(outcomes))) {
    fit <- fit_vdfr(dat_joint, outcomes$column[o], five)
    joint_tab <- bind_rows(joint_tab,
                           tidy_pfr(fit, five) %>%
                             mutate(response = outcomes$label[o], stage = "Joint", .before = 1))
    for (i in seq_along(five)) {
      save_plot(plot_heatmap(fit, select = i, n_bouts = ncol(dat_joint[[five[i]]]),
                             title = paste0(outcomes$label[o], ": ", five[i], " (joint)")),
                paste0("5stage_joint_", five[i], outcomes$suffix[o], ".pdf"), sub = "5stage")
    }
  }
  write_tables(c(format_coef(joint_tab), list(estimates = joint_tab, sample_size = sample_size)),
               file.path(table_dir, "5stage_joint_VDFR.xlsx"), note = stars_note)
}


################################################################################
# PART D: 3 STAGES, SEPARATE MODELS (W, N, REM)
################################################################################
# N = N1, N2 and N3 merged (consecutive NREM epochs form one bout).
# Outputs: result/tables/3stage_separate_VDFR.xlsx
#          result/plot/3stage/3stage_separate_<stage><suffix>.pdf

if (run_part_D) {
  cat("\nPART D: 3 stages, separate...\n")
  run_separate(c("W", "N", "REM"), function(s) final_3stage,
               "3 stages, separate", "3stage_separate", sub = "3stage")
}


################################################################################
# PART E: 2 STAGES, SEPARATE MODELS (wake, non-wake)
################################################################################
# non-wake = any sleep stage merged (consecutive sleep epochs form one bout).
# Outputs: result/tables/2stage_separate_VDFR.xlsx
#          result/plot/2stage/2stage_separate_<stage><suffix>.pdf

if (run_part_E) {
  cat("\nPART E: 2 stages, separate...\n")
  run_separate(c("wake", "non_wake"), function(s) final_2stage,
               "2 stages, separate", "2stage_separate", sub = "2stage")
}


################################################################################
# DONE
################################################################################
cat("\nAll analyses completed.\n  Tables:", table_dir, "\n  Plots: ", plot_dir, "\n")
sessionInfo()