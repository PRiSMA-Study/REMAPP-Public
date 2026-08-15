#*****************************************************************************
# ReMAPP Aim 2 public analysis code ----
#*****************************************************************************
# Description:
# Runs the primary ReMAPP Aim 2 Hb-outcome analysis.
# Automatically runs all trimesters and the combined analysis.
#
# The user only needs to specify:
#   1. Dataset
#   2. Outcome variable
#   3. Hb variable
#   4. Outcome type: binary or continuous (fatigue only)
#
# Created by:
# Precious Williams
# williams_pj@gwu.edu
#*****************************************************************************


# 1. Load packages ----
# Report generation
library(rmarkdown)
library(knitr)
library(officer)
library(flextable)
library(gt)

# Data wrangling & visualization
library(tidyverse)
library(ggrepel)
library(wesanderson)
library(ggbreak)
library(VennDiagram)

# Modeling
library(rms)

# Data I/O
library(haven)
library(webshot2)

# 2. Set file paths ----

## Data upload date
UploadDate <- "2026-01-30"

## TNT outcome data
path_to_tnt <- file.path(
  "Z:/Outcome Data",
  UploadDate
)

## Main ReMAPP analysis directory
# FOLDER SETUP (READ THIS BEFORE RUNNING!!!)
# Before running this script, choose a location on your computer where the
# ReMAPP analysis files and results should be stored.
#
# The analysis expects the following folder structure:
#
# ReMAPP/
# └── Aim2/
#     ├── iso_code/
#     └── iso_results/
#
# IMPORTANT:
# The "iso_code" folder contains the supporting isotonic regression, spline,
# plotting, and output functions required by this analysis.
#
# Copy the complete "iso_code" folder provided with the public analysis code
# on GitHub into your local Aim2 folder before running this script.
#
# The script will create the "Aim2" and "iso_results" folders automatically
# if they do not already exist. It will also create an empty "iso_code"
# folder if one is not present; however, the required R function files must
# be copied from GitHub into this folder before the analysis can be run.
#
# Update base_dir below to the location where you would like the ReMAPP
# analysis folder to be stored.

#Create your version!

base_dir <- "D:/Users/williams_pj/Documents/Analysis/ReMAPP/"

## Working directory
working_dir <- file.path(base_dir, "Aim2")

## Subfolders
iso_code_path    <- file.path(working_dir, "iso_code")
iso_results_path <- file.path(working_dir, "iso_results")

## Create testing folders if they do not already exist
dir.create(
  working_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

dir.create(
  iso_code_path,
  showWarnings = FALSE,
  recursive = TRUE
)

dir.create(
  iso_results_path,
  showWarnings = FALSE,
  recursive = TRUE
)

# Set working directory
setwd(working_dir)

# 3. Select analysis dataset ----
## Maternal outcomes (load if running maternal outcome)
analysis_data <- read.csv(
  file.path(path_to_tnt, "MAT_REMAPP_AIM2.csv"),
  stringsAsFactors = FALSE
)

## Infant outcomes (load if running infant outcome)
# analysis_data <- read.csv(
#   file.path(path_to_tnt, "INF_REMAPP_AIM2.csv"),
#   stringsAsFactors = FALSE
# )


# 4. Specify analysis variables ----

## Change only these values for each analysis
#This will include the variable, the labels
# Preeclampsia
outcome_var <- "PRECLAMP"
hb_var      <- "PRECLAMP_HB"
outcome_type <- "binary"
outcome_label <- "Preeclampsia"
outcome_file_label <- "preclamp"


## Examples of outcomes ----
# Copy and paste the desired block above
# ─── INFANT OUTCOMES (Binary) ────────

# Preterm birth <37 weeks
# outcome_var <- "PRETERM37"
# hb_var      <- "PRETERM37_HB"
# outcome_type <- "binary"
# outcome_label <- "Preterm birth <37 weeks"
# outcome_file_label <- "preterm37"

# Preterm birth <34 weeks
# outcome_var <- "PRETERM34"
# hb_var      <- "PRETERM34_HB"
# outcome_type <- "binary"
# outcome_label <- "Preterm birth <34 weeks"
# outcome_file_label <- "preterm34"

# Low birth weight <2500g
# outcome_var <- "LBW2500"
# hb_var      <- "LBW2500_HB"
# outcome_type <- "binary"
# outcome_label <- "Low birth weight <2500g"
# outcome_file_label <- "lbw2500"

# Very low birth weight <1500g
# outcome_var <- "LBW1500"
# hb_var      <- "LBW1500_HB"
# outcome_type <- "binary"
# outcome_label <- "Very low birth weight <1500g"
# outcome_file_label <- "lbw1500"

# Small for gestational age <10th percentile
# outcome_var <- "SGA10"
# hb_var      <- "SGA10_HB"
# outcome_type <- "binary"
# outcome_label <- "SGA <10th percentile"
# outcome_file_label <- "sga10"

# Small for gestational age <3rd percentile
# outcome_var <- "SGA3"
# hb_var      <- "SGA3_HB"
# outcome_type <- "binary"
# outcome_label <- "SGA <3rd percentile"
# outcome_file_label <- "sga3"

# Composite outcome (preterm, LBW, or SGA)
# outcome_var <- "COMPOSITE"
# hb_var      <- "COMPOSITE_HB"
# outcome_type <- "binary"
# outcome_label <- "Composite (preterm/LBW/SGA)"
# outcome_file_label <- "composite"

# Possible severe bacterial infection
# outcome_var <- "PSBI"
# hb_var      <- "PSBI_HB"
# outcome_type <- "binary"
# outcome_label <- "Possible severe bacterial infection"
# outcome_file_label <- "psbi"

# Neonatal asphyxia
# outcome_var <- "ASPH"
# hb_var      <- "ASPH_HB"
# outcome_type <- "binary"
# outcome_label <- "Neonatal asphyxia"
# outcome_file_label <- "asph"

# Stillbirth ≥20 weeks
# outcome_var <- "STILLBIRTH20"
# hb_var      <- "STILLBIRTH20_HB"
# outcome_type <- "binary"
# outcome_label <- "Stillbirth ≥20 weeks"
# outcome_file_label <- "stillbirth20"

# Stillbirth ≥28 weeks
# outcome_var <- "STILLBIRTH28"
# hb_var      <- "STILLBIRTH28_HB"
# outcome_type <- "binary"
# outcome_label <- "Stillbirth ≥28 weeks"
# outcome_file_label <- "stillbirth28"

# Neonatal hyperbilirubinemia
# outcome_var <- "HYPERBILI"
# hb_var      <- "HYPERBILI_HB"
# outcome_type <- "binary"
# outcome_label <- "Neonatal hyperbilirubinemia"
# outcome_file_label <- "hyperbili"

# Neonatal mortality (0-28 days)
# outcome_var <- "NEO_MORTALITY"
# hb_var      <- "NEO_MORTALITY_HB"
# outcome_type <- "binary"
# outcome_label <- "Neonatal mortality"
# outcome_file_label <- "neo_mortality"


# ─── MATERNAL OUTCOMES  ────

# Maternal composite outcome
# outcome_var <- "MAT_COMPO"
# hb_var      <- "MAT_COMPO_HB"
# outcome_type <- "binary"
# outcome_label <- "Maternal composite outcome"
# outcome_file_label <- "mat_compo"

# Postpartum hemorrhage
# outcome_var <- "HEM_PPH"
# hb_var      <- "HEM_PPH_HB"
# outcome_type <- "binary"
# outcome_label <- "Postpartum hemorrhage"
# outcome_file_label <- "hem_pph"

# Preterm premature rupture of membranes
# outcome_var <- "PPROM"
# hb_var      <- "PPROM_HB"
# outcome_type <- "binary"
# outcome_label <- "Preterm premature rupture of membranes"
# outcome_file_label <- "pprom"

# Preeclampsia
# outcome_var <- "PRECLAMP"
# hb_var      <- "PRECLAMP_HB"
# outcome_type <- "binary"
# outcome_label <- "Preeclampsia"
# outcome_file_label <- "preclamp"

# Postpartum anemia at PNC6
# outcome_var <- "PPA_PNC6"
# hb_var      <- "PPA_PNC6_HB"
# outcome_type <- "binary"
# outcome_label <- "Postpartum anemia at PNC6"
# outcome_file_label <- "ppa_pnc6"

# Postpartum anemia at PNC26
# outcome_var <- "PPA_PNC26"
# hb_var      <- "PPA_PNC26_HB"
# outcome_type <- "binary"
# outcome_label <- "Postpartum anemia at PNC26"
# outcome_file_label <- "ppa_pnc26"

# Depression (binary)
# outcome_var <- "DEPRESS"
# hb_var      <- "DEPRESS_HB"
# outcome_type <- "binary"
# outcome_label <- "Depression (EPDS ≥13)"
# outcome_file_label <- "depression"

# Fatigue (continuous only)
# outcome_var <- "FATIGUE"
# hb_var      <- "FATIGUE_HB"
# outcome_type <- "continuous"
# outcome_label <- "Fatigue score"
# outcome_file_label <- "fatigue"


# 5. Source analysis functions ----
## Binary outcomes
if (outcome_type == "binary") {
  
  source(file.path(iso_code_path, "spline_binary.R"))
  source(file.path(iso_code_path, "iso_binary.R"))
  source(file.path(iso_code_path, "plot_figure.R"))
  source(file.path(iso_code_path, "out_table.R"))
}

## Continuous outcome: fatigue only
if (outcome_type == "continuous") {
  
  source(file.path(iso_code_path, "spline_continuous.R"))
  source(file.path(iso_code_path, "iso_continuous.R"))
  source(file.path(iso_code_path, "plot_continuous.R"))
  source(file.path(iso_code_path, "out_table.R"))
}

# 6. Run analysis function ----
## Site variable
site_var <- "SITE"
## Hb analysis range
hb_min <- 5
hb_max <- 18

run_analysis <- function(df, label, suffix, outcome_var_override = NULL, hb_var_override = NULL) {
  
  if (is.null(df) || nrow(df) == 0) {
    cat("⚠️ Skipping", label, "- no data\n")
    return(NULL)
  }
  
  # Use override variables if provided (for trimester-specific outcomes)
  use_outcome <- if (!is.null(outcome_var_override)) outcome_var_override else outcome_var
  use_hb <- if (!is.null(hb_var_override)) hb_var_override else hb_var
  
  # Prepare data - filter out missing Hb for this specific variable
  analysis_df <- df %>%
    transmute(
      SITE    = .data[[site_var]],
      outcome = as.numeric(.data[[use_outcome]]),
      hb      = as.numeric(.data[[use_hb]])
    ) %>%
    filter(!is.na(outcome), !is.na(hb), !is.na(SITE)) %>%
    mutate(hb = round(pmin(pmax(hb, hb_min), hb_max), 1))
  
  if (nrow(analysis_df) == 0) {
    cat("⚠️ Skipping", label, "- no complete cases\n")
    return(NULL)
  }
  
  # CHECK: Skip if N < 100 (too few observations for reliable analysis)
  if (nrow(analysis_df) < 100) {
    cat("⚠️ Skipping", label, "- N =", nrow(analysis_df), "(< 100, insufficient sample size)\n")
    return(NULL)
  }
  
  cat("\n", paste(rep("─", 60), collapse = ""), "\n")
  cat("📊 Running:", label, "| N =", nrow(analysis_df), "\n")
  cat(paste(rep("─", 60), collapse = ""), "\n")
  
  # Fit spline
  if (outcome_type == "binary") {
    spline_result <- knot_fun_boot(analysis_df, "hb", "outcome")
  } else {
    spline_result <- knot_fun_boot(analysis_df, "hb", "outcome")
  }
  
  # Fit 0.50 isotonic
  if (outcome_type == "binary") {
    iso_result <- flexstepreg_glmer(analysis_df$outcome, analysis_df$hb, analysis_df$SITE,
                                    covar2 = NULL, random_effect = NULL, alpha = 0.01)
  } else {
    iso_result <- flexstepreg_lmer(analysis_df$outcome, analysis_df$hb, analysis_df$SITE,
                                   covar2 = NULL, random_effect = analysis_df$MOMID, alpha = 0.01)
  }
  
  # Save models
  saveRDS(spline_result, file.path("iso_results", paste0("spline_", suffix, ".rds")))
  saveRDS(iso_result, file.path("iso_results", paste0("iso_", suffix, "_50.rds")))
  
  # Generate plot
  png(file.path("iso_results", paste0("vioplot_", suffix, "_50.png")), width = 800, height = 600)
  
  if (outcome_type == "binary") {
    plot_boot_violin(x = analysis_df$hb, y = analysis_df$outcome,
                     xlab = "Hemoglobin (g/dL)", ylab = outcome_label,
                     rcs_result = spline_result, iso_model = iso_result,
                     outcome_var = analysis_df$outcome,
                     title = paste0(label, " (0.50 Isotonic)"))
  } else {
    plot_boot_continuous(x = analysis_df$hb, y = analysis_df$outcome,
                         xlab = "Hemoglobin (g/dL)", ylab = outcome_label,
                         rcs_result = spline_result, iso_model = iso_result,
                         outcome_var = analysis_df$outcome,
                         title = paste0(label, " (0.50 Isotonic)"))
  }
  dev.off()
  
  # Save output
  if (outcome_type == "binary") {
    iso_output <- outdata(iso_result, paste0(label, " (.50)"))
  } else {
    iso_output <- outdata(iso_result, paste0(label, " (.50)"))
  }
  
  saveRDS(iso_output, file.path("iso_results", paste0("out_", suffix, "_50.rds")))
  
  # Return results
  return(list(
    df = analysis_df,
    spline = spline_result,
    iso = iso_result,
    output = iso_output,
    n = nrow(analysis_df)
  ))
}


# 7. Run all analyses ----
dir.create("iso_results", showWarnings = FALSE, recursive = TRUE)
results <- list()

# Combined (all trimesters)
results$all <- run_analysis(analysis_data, outcome_label, outcome_file_label)

# Trimester 1 ----
# For regular outcomes: outcome stays same, Hb changes
# For depression/fatigue: both outcome and Hb change
hb_trim1 <- paste0(hb_var, "_TRIM1")

# Check if this is depression or fatigue (outcome changes by trimester)
if (outcome_file_label %in% c("depression", "fatigue")) {
  out_trim1 <- paste0(outcome_var, "_TRIM1")
} else {
  out_trim1 <- outcome_var  # outcome stays the same
}

if (!hb_trim1 %in% names(analysis_data)) {
  hb_trim1 <- paste0(sub("_HB$", "", hb_var), "_HB_TRIM1")
}

if (hb_trim1 %in% names(analysis_data)) {
  results$trim1 <- run_analysis(analysis_data, paste0(outcome_label, " (Trim1)"), 
                                paste0(outcome_file_label, "_trim1"),
                                outcome_var_override = out_trim1,
                                hb_var_override = hb_trim1)
}

# Trimester 2 ----
hb_trim2 <- paste0(hb_var, "_TRIM2")

if (outcome_file_label %in% c("depression", "fatigue")) {
  out_trim2 <- paste0(outcome_var, "_TRIM2")
} else {
  out_trim2 <- outcome_var
}

if (!hb_trim2 %in% names(analysis_data)) {
  hb_trim2 <- paste0(sub("_HB$", "", hb_var), "_HB_TRIM2")
}

if (hb_trim2 %in% names(analysis_data)) {
  results$trim2 <- run_analysis(analysis_data, paste0(outcome_label, " (Trim2)"), 
                                paste0(outcome_file_label, "_trim2"),
                                outcome_var_override = out_trim2,
                                hb_var_override = hb_trim2)
}

# Trimester 3 ----
hb_trim3 <- paste0(hb_var, "_TRIM3")

if (outcome_file_label %in% c("depression", "fatigue")) {
  out_trim3 <- paste0(outcome_var, "_TRIM3")
} else {
  out_trim3 <- outcome_var
}

if (!hb_trim3 %in% names(analysis_data)) {
  hb_trim3 <- paste0(sub("_HB$", "", hb_var), "_HB_TRIM3")
}

if (hb_trim3 %in% names(analysis_data)) {
  results$trim3 <- run_analysis(analysis_data, paste0(outcome_label, " (Trim3)"), 
                                paste0(outcome_file_label, "_trim3"),
                                outcome_var_override = out_trim3,
                                hb_var_override = hb_trim3)
}

# PNC6 (for relevant outcomes: PPA_PNC26, DEPRESS, FATIGUE) ----
hb_pnc6 <- paste0(hb_var, "_PNC6")

if (outcome_file_label %in% c("depression", "fatigue")) {
  out_pnc6 <- paste0(outcome_var, "_PNC6")
} else {
  out_pnc6 <- outcome_var
}

if (!hb_pnc6 %in% names(analysis_data)) {
  hb_pnc6 <- paste0(sub("_HB$", "", hb_var), "_HB_PNC6")
}

if (hb_pnc6 %in% names(analysis_data)) {
  results$pnc6 <- run_analysis(analysis_data, paste0(outcome_label, " (PNC6)"), 
                               paste0(outcome_file_label, "_pnc6"),
                               outcome_var_override = out_pnc6,
                               hb_var_override = hb_pnc6)
}
# 8. Summary results ----
cat("\n", paste(rep("═", 80), collapse = ""), "\n")
cat("📊 SUMMARY RESULTS\n")
cat(paste(rep("═", 80), collapse = ""), "\n")

# Extract and display outdata results in the original outdata format
for (name in names(results)) {
  if (!is.null(results[[name]])) {
    label <- switch(name,
                    all = "Combined",
                    trim1 = "Trimester 1",
                    trim2 = "Trimester 2",
                    trim3 = "Trimester 3",
                    pnc6 = "PNC6",
                    name)
    
    cat("\n", paste(rep("─", 60), collapse = ""), "\n")
    cat("📋", label, "| N =", results[[name]]$n, "\n")
    cat(paste(rep("─", 60), collapse = ""), "\n")
    
    # Print the outdata table in its original format
    print(results[[name]]$output)
  }
}

# Sample size summary (only including analyses with N >= 100)
n_summary <- tibble(
  Analysis = names(results),
  N = map_dbl(results, ~ ifelse(is.null(.x), NA, .x$n))
) %>% 
  filter(!is.na(N)) %>%
  arrange(desc(N))

cat("\n", paste(rep("═", 80), collapse = ""), "\n")
cat("📋 SAMPLE SIZE SUMMARY (Analyses with N >= 100)\n")
cat(paste(rep("─", 80), collapse = ""), "\n")
print(n_summary)

# List any skipped analyses (N < 100)
skipped <- names(results)[sapply(results, is.null)]
if (length(skipped) > 0) {
  cat("\n⚠️ SKIPPED ANALYSES (N < 100 or no data):\n")
  for (s in skipped) {
    label <- switch(s,
                    all = "Combined",
                    trim1 = "Trimester 1",
                    trim2 = "Trimester 2",
                    trim3 = "Trimester 3",
                    pnc6 = "PNC6",
                    s)
    cat("  -", label, "\n")
  }
}

# SHOW FULL PATHS OF SAVED FILES
cat("\n", paste(rep("═", 80), collapse = ""), "\n")
cat("📁 FILES SAVED\n")
cat(paste(rep("─", 80), collapse = ""), "\n")
cat("Location:", iso_results_path, "\n\n")

saved_files <- list.files(iso_results_path, pattern = outcome_file_label, full.names = TRUE)
if (length(saved_files) > 0) {
  for (f in saved_files) {
    cat("  -", f, "\n")
  }
} else {
  cat("  No files found with pattern:", outcome_file_label, "\n")
}

cat("\n", paste(rep("═", 80), collapse = ""), "\n")
cat("✅ ANALYSIS COMPLETE\n")
cat(paste(rep("═", 80), collapse = ""), "\n")