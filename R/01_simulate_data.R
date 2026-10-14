################################################################################
#  01_simulate_data.R
#  Simulated data for the OSA analysis (AHI vs. sleep-bout durations)
#  ---------------------------------------------------------------------------
#  Project layout (run from the project folder, or from R/):
#      <project>/
#        R/        01_simulate_data.R, 02_analysis.R
#        data/     analysis data sets + demographics.csv (written here)
#          epoch/  one sleep-stage file per participant (written here)
#        result/   tables/ and plot/ (written by 02_analysis.R)
#
#  What this script writes
#    data/epoch/<ID>.csv     the participant's night, one row per bout
#                            (consecutive 30-s epochs of the same stage):
#                              stage       W, N1, N2, N3 or REM
#                              start_time  clock time the bout starts (hh:mm:ss)
#                              end_time    clock time the bout ends   (hh:mm:ss)
#                              time        bout length in seconds
#    data/demographics.csv   one row per participant: ID, the model
#                            covariates (gender, age, bmiblquan, ethnicity,
#                            maritalstatushp, highestgradehp) and the outcome ahi
#    data/data_5stage.rda    final_5stage: list (W, N1, N2, N3, REM); element s =
#                            participants with at least one bout of stage s,
#                            with all five bout-duration matrices
#    data/data_3stage.rda    final_3stage: matrices N (N1-N3 merged), REM, W
#    data/data_2stage.rda    final_2stage: matrices wake, non_wake (all sleep merged)
#    data/SIMULATED_DATA.txt marker; the script never overwrites a data/ folder
#                            without it
#  The bout-duration matrices are built as the real ones were: the first bout
#  is dropped if it is W and the last bout if it is W (wake before sleep onset
#  and after final awakening), durations are in minutes, rows are NA-padded.
#
#  How the data are generated
#  1. Covariates     : drawn to match the real sample (sex, age, BMI,
#                      education, marital status, ethnicity).
#  2. One night each : a semi-Markov model of sleep stages on 30-s epochs.
#                      Wake before sleep onset -> alternating sleep-stage bouts
#                      (next stage from time-of-night-dependent transition
#                      weights; bout length lognormal) -> final wake.
#                      N3 is more likely early in the night, REM late, and
#                      about half the participants have no N3 at all, as in
#                      the real data. A participant-level "fragmentation"
#                      score shortens all bouts.
#  3. Outcome        : log(AHI) = alpha + Z'gamma + [sleep term] + error,
#                      AHI = exp(log(AHI)). The sleep term depends on `scenario`:
#        "fragmentation" log(AHI) rises with the fragmentation score, which
#                        also shortens bouts -> the bout-duration models have
#                        a real signal to find.
#        "null"          no sleep term -> no association with bout durations.
################################################################################

rm(list = ls())


################################################################################
# S0. SETTINGS
################################################################################

pacman::p_load(tidyverse)

# ---- Project folder: the folder that holds R/, data/ and result/ --------------
project_dir <- getwd()
if (basename(project_dir) == "R") project_dir <- dirname(project_dir)
data_dir  <- file.path(project_dir, "data")
epoch_dir <- file.path(data_dir, "epoch")

scenario <- "fragmentation"   # "fragmentation" or "null"
n_sim    <- 150            # participants (same as the real full sample)
seed     <- 2026

# Share of the variance of log(AHI) explained by the sleep term
signal_share <- 0.10

# Stage coding used throughout: 1 = W, 2 = N1, 3 = N2, 4 = N3, 5 = REM
stage_names <- c("W", "N1", "N2", "N3", "REM")

# Bout length (in 30-s epochs): lognormal with these medians and log-SDs
bout_median <- c(W = 3,   N1 = 3,   N2 = 7.5, N3 = 4.5, REM = 12)
bout_sdlog  <- c(W = 1.0, N1 = 0.9, N2 = 1.2, N3 = 1.3, REM = 1.0)

# Transition weights: row = current stage, column = next stage (self = 0).
# N3 and REM columns are further scaled by time of night (see S2).
exit_w <- rbind(
  W   = c(0,    0.70, 0.25, 0.00, 0.05),
  N1  = c(0.35, 0,    0.55, 0.00, 0.10),
  N2  = c(0.30, 0.40, 0,    0.20, 0.10),
  N3  = c(0.20, 0.10, 0.70, 0,    0.00),
  REM = c(0.40, 0.40, 0.20, 0.00, 0)
)
n3_time  <- function(t) 3 * exp(-4 * t) + 0.1     # N3 front-loaded
rem_time <- function(t) 0.3 + 1.5 * t             # REM back-loaded


################################################################################
# S1. COVARIATES (matched to the real sample)
################################################################################

simulate_covariates <- function(n) {
  tibble(
    ID              = sprintf("ID%04d", seq_len(n)),
    gender          = factor(if_else(rbinom(n, 1, 0.655) == 1, "Male", "Female"),
                             levels = c("Female", "Male")),
    age             = round(pmin(pmax(rnorm(n, 51.6, 12.2), 18), 85), 1),
    bmiblquan       = round(pmax(18, rnorm(n, 32.2, 7.1)), 2),
    ethnicity       = factor(sample(c("Asian", "Black", "White", "Other"), n, TRUE,
                                    prob = c(0.055, 0.0925, 0.7616, 0.0909)),
                             levels = c("Other", "White", "Black", "Asian")),
    maritalstatushp = factor(sample(c("Single", "Married", "Other"), n, TRUE,
                                    prob = c(0.2375, 0.5793, 0.1831)),
                             levels = c("Other", "Single", "Married")),
    highestgradehp  = if_else(runif(n) < 0.02, NA_real_,
                              round(pmin(pmax(rnorm(n, 15.5, 2.6), 8), 22)))
  )
}


################################################################################
# S2. ONE NIGHT PER PARTICIPANT (semi-Markov on 30-s epochs)
################################################################################
# Returns an n x max_epochs integer matrix of stage codes (NA after the end
# of each participant's recording) and the recording length in epochs.

simulate_nights <- function(n, frag, has_n3) {
  len <- round(runif(n, 7.5, 9) * 120)                            # 7.5-9 h
  lat <- pmax(2, round(rlnorm(n, log(30), 0.6)))                  # onset latency, ~15 min
  fwk <- pmax(1, round(rlnorm(n, log(20), 0.7)))                  # final wake, ~10 min
  sleep_end <- len - fwk

  E <- matrix(NA_integer_, n, max(len))
  for (i in seq_len(n)) E[i, seq_len(len[i])] <- 1L               # start all-wake

  pos    <- lat                    # last filled epoch (wake up to sleep onset)
  cur    <- rep(1L, n)             # current stage
  active <- pos < sleep_end
  first  <- TRUE

  while (any(active)) {
    a  <- which(active)
    tt <- pos[a] / len[a]                                          # time of night in [0, 1]

    # next stage: first sleep bout is N1 or N2, otherwise transition weights
    if (first) {
      nxt <- ifelse(runif(length(a)) < 0.7, 2L, 3L)
      first <- FALSE
    } else {
      w <- exit_w[cur[a], , drop = FALSE]
      w[, 4] <- w[, 4] * n3_time(tt) * has_n3[a]
      w[, 5] <- w[, 5] * rem_time(tt)
      cp <- w
      for (j in 2:5) cp[, j] <- cp[, j - 1] + w[, j]
      cp  <- cp / cp[, 5]
      nxt <- 1L + rowSums(runif(length(a)) > cp)
    }

    # bout length: lognormal, shorter for more fragmented sleepers
    d <- pmax(1L, as.integer(round(rlnorm(length(a),
                                          log(bout_median[nxt]) - 0.3 * frag[a],
                                          bout_sdlog[nxt]))))
    end <- pmin(pos[a] + d, sleep_end[a])
    for (k in seq_along(a)) {
      i <- a[k]
      if (end[k] > pos[i]) E[i, (pos[i] + 1):end[k]] <- nxt[k]
    }
    pos[a] <- end
    cur[a] <- nxt
    active <- pos < sleep_end
  }
  list(E = E, len = len)
}


################################################################################
# S3. HELPERS: bouts and matrices
################################################################################

# Epochs from sleep onset to final awakening: the first bout dropped if it is
# W and the last bout dropped if it is W
sleep_period <- function(x) {
  sl <- which(x != 1L)
  if (length(sl) == 0) return(integer(0))
  x[min(sl):max(sl)]
}

# Bout durations (minutes) inside the sleep period, merging consecutive
# epochs whose codes are all in `codes`
bouts_in_period <- function(x, codes) {
  x <- sleep_period(x)
  if (length(x) == 0) return(numeric(0))
  r <- rle(x %in% codes)
  r$lengths[r$values] * 0.5
}

# List of duration vectors -> NA-padded matrix (row = participant)
pad_matrix <- function(v) {
  J <- max(1, max(lengths(v)))
  m <- matrix(NA_real_, length(v), J)
  for (i in seq_along(v)) if (length(v[[i]])) m[i, seq_along(v[[i]])] <- v[[i]]
  m
}

clock <- function(sec) {
  sec <- sec %% 86400
  sprintf("%02d:%02d:%02d", sec %/% 3600, (sec %% 3600) %/% 60, sec %% 60)
}


################################################################################
# S4. GENERATE
################################################################################

set.seed(seed)
cat("Simulating", n_sim, "participants, scenario:", scenario, "\n")

# ---- S4.1 Safety: never overwrite a data/ folder this script did not create ------
marker <- file.path(data_dir, "SIMULATED_DATA.txt")
if (dir.exists(data_dir) && length(list.files(data_dir)) > 0 && !file.exists(marker)) {
  stop(data_dir, " already has files that were not made by this script. ",
       "Use an empty project folder.", call. = FALSE)
}
dir.create(epoch_dir, recursive = TRUE, showWarnings = FALSE)
old_files <- list.files(epoch_dir, pattern = "\\.csv$", full.names = TRUE)
if (length(old_files)) file.remove(old_files)          # earlier simulated files only
writeLines(c("Simulated data written by R/01_simulate_data.R",
             paste("scenario:", scenario), paste("seed:", seed), paste("n:", n_sim),
             paste("created:", Sys.time())), marker)

# ---- S4.2 Covariates and nights ----------------------------------------------
cov    <- simulate_covariates(n_sim)
frag   <- rnorm(n_sim)                                         # fragmentation score
has_n3 <- rbinom(n_sim, 1, plogis(0.1 - 0.04 * (cov$age - 51.6) - 0.3 * frag))
night  <- simulate_nights(n_sim, frag, has_n3)
ep <- map(seq_len(n_sim), ~ night$E[.x, seq_len(night$len[.x])])   # one vector per person

# ---- S4.3 Sleep-stage files: data/epoch/<ID>.csv (one row per bout) ---------------
lights_off <- 22 * 3600 + 30 * sample(0:239, n_sim, TRUE)      # 22:00-24:00
for (i in seq_len(n_sim)) {
  r     <- rle(ep[[i]])
  secs  <- 30 * r$lengths
  start <- lights_off[i] + c(0, cumsum(secs)[-length(secs)])
  write_csv(tibble(stage      = stage_names[r$values],
                   start_time = clock(start),
                   end_time   = clock(start + secs),
                   time       = secs),
            file.path(epoch_dir, paste0(cov$ID[i], ".csv")))
}

# ---- S4.4 Outcome ------------------------------------------------------------
# covariate part: effects close to the fitted ones on the real data
gamma <- c(male = 0.29, age = 0.010, bmi = 0.032, educ = -0.020)
lin   <- gamma["male"] * (cov$gender == "Male") + gamma["age"] * cov$age +
         gamma["bmi"] * cov$bmiblquan + gamma["educ"] * coalesce(cov$highestgradehp, 15.5)
total_var <- 0.64^2                                           # SD of log(AHI) in the real data

sleep_term  <- if (scenario == "fragmentation") frag else rep(0, n_sim)
scale_sleep <- if (scenario == "null") 0 else sqrt(signal_share * total_var / var(sleep_term))
sleep_term  <- scale_sleep * sleep_term
err_sd <- sqrt(max(total_var - var(lin) - var(sleep_term), 0.05))

log_ahi <- lin + sleep_term + rnorm(n_sim, 0, err_sd)
log_ahi <- log_ahi + (3.50 - mean(log_ahi))                   # mean log(AHI) as in the real data
cov$ahi <- round(exp(log_ahi), 1)

write_csv(cov, file.path(data_dir, "demographics.csv"))

# ---- S4.5 Bout-duration matrices ----------------------------------------------
mats <- list(
  W        = map(ep, bouts_in_period, codes = 1L),             # wake inside sleep period
  N1       = map(ep, bouts_in_period, codes = 2L),
  N2       = map(ep, bouts_in_period, codes = 3L),
  N3       = map(ep, bouts_in_period, codes = 4L),
  REM      = map(ep, bouts_in_period, codes = 5L),
  N        = map(ep, bouts_in_period, codes = 2:4),            # NREM merged
  wake     = map(ep, bouts_in_period, codes = 1L),
  non_wake = map(ep, bouts_in_period, codes = 2:5)             # any sleep, merged
) %>% map(pad_matrix)
has <- map(mats, ~ rowSums(!is.na(.x)) > 0)

add_mats <- function(df, idx, names) {
  for (nm in names) df[[nm]] <- I(mats[[nm]][idx, , drop = FALSE])
  df
}

# 5 stages: everyone with wake inside the sleep period ($W), with all five
# matrices; element s = the subset with at least one bout of stage s
five     <- c("W", "N1", "N2", "N3", "REM")
base_idx <- which(has$W)
all5     <- add_mats(cov[base_idx, ], base_idx, five)
final_5stage <- map(set_names(five), ~ all5[has[[.x]][base_idx], ])

# 3 stages / 2 stages: participants with every stage of that grouping
idx3 <- which(has$N & has$REM & has$W)
final_3stage <- add_mats(cov[idx3, ], idx3, c("N", "REM", "W"))
idx2 <- which(has$wake & has$non_wake)
final_2stage <- add_mats(cov[idx2, ], idx2, c("wake", "non_wake"))

save(final_5stage, file = file.path(data_dir, "data_5stage.rda"))
save(final_3stage, file = file.path(data_dir, "data_3stage.rda"))
save(final_2stage, file = file.path(data_dir, "data_2stage.rda"))


################################################################################
# S5. SUMMARY: does it look like the real data?
################################################################################

cc <- base_idx[Reduce(`&`, has[five])[base_idx]]
cat("\nSample sizes (real data in brackets):\n")
cat("  5-stage (W):", nrow(final_5stage$W), "[1103]   with N3:", nrow(final_5stage$N3),
    "[545]   with REM:", nrow(final_5stage$REM), "[1093]   all five stages:", length(cc), "[541]\n")
cat("  3-stage:", nrow(final_3stage), "  2-stage:", nrow(final_2stage), "\n")
cat("Median bouts per person with all five stages (real: W 23, N1 34, N2 35, N3 4, REM 7):\n")
print(map_dbl(set_names(five), ~ median(rowSums(!is.na(mats[[.x]][cc, , drop = FALSE])))))
cat("AHI: mean", round(mean(cov$ahi), 2), "SD", round(sd(cov$ahi), 2),
    "[40.09, 25.31]   log(AHI): mean", round(mean(log(cov$ahi)), 2),
    "SD", round(sd(log(cov$ahi)), 2), "[3.50, 0.64]\n")
cat("\nWritten to", data_dir, "\n  ", length(list.files(epoch_dir)), "files in data/epoch\n")
