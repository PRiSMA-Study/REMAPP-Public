# ───────────────────────────────────────
# 📁 ReMAPP Aim 2 – PSBI Analysis Pipeline
# 🔧 Modified by: Williams Precious
# 📅 Last updated: 2025-04-18
# 📄 Description: Fits spline and isotonic models for PSBI, generates predictions and plots.
# ───────────────────────────────────────
rm(list = ls())

# ───────────────────────────────────────
# 📦 Load Libraries
# ───────────────────────────────────────
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


# ───────────────────────────────────────
# Setup Paths & Load Data ----
# ───────────────────────────────────────
UploadDate <- "2026-01-30"

base_dir <- file.path("D:/Users/williams_pj/Documents/Analysis/ReMAPP/Aim2", UploadDate)
setwd(base_dir)

load("derived_data/df_inf_psbi.rda")
load("derived_data/df_inf_psbi_trim1.rda")
load("derived_data/df_inf_psbi_trim2.rda")
load("derived_data/df_inf_psbi_trim3.rda")

df_inf_psbi$hb        <- round(pmin(pmax(df_inf_psbi$hb, 5), 18), 1)
df_inf_psbi_trim1$hb  <- round(pmin(pmax(df_inf_psbi_trim1$hb, 5), 18), 1)
df_inf_psbi_trim2$hb  <- round(pmin(pmax(df_inf_psbi_trim2$hb, 5), 18), 1)
df_inf_psbi_trim3$hb  <- round(pmin(pmax(df_inf_psbi_trim3$hb, 5), 18), 1)

# Establish the source of the functions
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/spline_binary_boot.R")
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/iso_binary_25_mod.R")
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/iso_binary_50_mod.R")
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/plot_boot.R")
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/outdata_mod.R")

# ───────────────────────────────────────
#ALL TRIMESTERS PSBI ----
# ───────────────────────────────────────

# ───────────────────────────────────────
## Fit Spline Models ----
# ───────────────────────────────────────
spline_psbi <- knot_fun_boot(df_inf_psbi, "hb", "inf_psbi")
saveRDS(spline_psbi, "iso_results/spline_psbi.rds")

# ───────────────────────────────────────
## Fit Isotonic Models ----
# ───────────────────────────────────────

iso_psbi_25 <- flexstepreg_glmer_25(df_inf_psbi$inf_psbi, df_inf_psbi$hb, df_inf_psbi$SITE, 
                                         covar2 = NULL, random_effect = NULL, alpha = 0.01)
saveRDS(iso_psbi_25, "iso_results/iso_psbi_25.rds")

iso_psbi_50 <- flexstepreg_glmer(df_inf_psbi$inf_psbi, df_inf_psbi$hb, df_inf_psbi$SITE, 
                                      covar2 = NULL, random_effect = NULL, alpha = 0.01)
saveRDS(iso_psbi_50, "iso_results/iso_psbi_50.rds")

# ───────────────────────────────────────
## Generate Plots ----
# ───────────────────────────────────────

# Save 0.50 plot
png("iso_results/vioplot_psbi_50.png", width = 800, height = 600)
# 0.50 Isotonic
vioplot_psbi_50 <- plot_boot_violin(
  x = df_inf_psbi$hb,
  y = df_inf_psbi$inf_psbi,
  xlab = "Hb",
  ylab = "PSBI",
  rcs_result = spline_psbi,
  iso_model = iso_psbi_50,
  outcome_var = df_inf_psbi$inf_psbi,
  title = "PSBI with Bootstrap Spline CI (0.50 Isotonic)"
)
dev.off()

# Save 0.25 plot
png("iso_results/vioplot_psbi_25.png", width = 800, height = 600)
# 0.25 Isotonic
vioplot_psbi_25 <- plot_boot_violin(
  x = df_inf_psbi$hb,
  y = df_inf_psbi$inf_psbi,
  xlab = "Hb",
  ylab = "PSBI",
  rcs_result = spline_psbi,
  iso_model = iso_psbi_25,
  outcome_var = df_inf_psbi$inf_psbi,
  title = "PSBI with Bootstrap Spline CI (0.25 Isotonic)"
)
dev.off()


# ───────────────────────────────────────
## Save Processed Isotonic Output Objects ----
# ───────────────────────────────────────

# Helper to forinf and save isotonic model outputs
save_iso_output <- function(model_obj, label, filename) {
  out <- outdata(model_obj, label)
  save(out, file = file.path("iso_results", filename))
  return(out)
}

# Save outputs
out_psbi_25  <- save_iso_output(iso_psbi_25,  "PSBI (.25)",  "out_psbi_25.rda")
out_psbi_50  <- save_iso_output(iso_psbi_50,  "PSBI (.50)",  "out_psbi_50.rda")


# ───────────────────────────────────────
#TRIMESTER 1 – PSBI Models ----
# ───────────────────────────────────────

# ───────────────────────────────────────
## Fit Spline Models ----
# ───────────────────────────────────────
spline_psbi_trim1 <- knot_fun_boot(df_inf_psbi_trim1, "hb", "inf_psbi")
saveRDS(spline_psbi_trim1, "iso_results/spline_psbi_trim1.rds")

# ───────────────────────────────────────
## Fit Isotonic Models ----
# ───────────────────────────────────────

iso_psbi_25_trim1 <- flexstepreg_glmer_25(df_inf_psbi_trim1$inf_psbi, df_inf_psbi_trim1$hb, df_inf_psbi_trim1$SITE, 
                                               covar2 = NULL, random_effect = NULL, alpha = 0.01)
saveRDS(iso_psbi_25_trim1, "iso_results/iso_psbi_25_trim1.rds")

iso_psbi_50_trim1 <- flexstepreg_glmer(df_inf_psbi_trim1$inf_psbi, df_inf_psbi_trim1$hb, df_inf_psbi_trim1$SITE, 
                                            covar2 = NULL, random_effect = NULL, alpha = 0.01)
saveRDS(iso_psbi_50_trim1, "iso_results/iso_psbi_50_trim1.rds")

# ───────────────────────────────────────
## Generate Plots ----
# ───────────────────────────────────────

# Save 0.50 plot
png("iso_results/vioplot_psbi_50_trim1.png", width = 800, height = 600)
# 0.50 Isotonic
vioplot_psbi_50_trim1 <- plot_boot_violin(
  x = df_inf_psbi_trim1$hb,
  y = df_inf_psbi_trim1$inf_psbi,
  xlab = "Hb",
  ylab = "PSBI",
  rcs_result = spline_psbi_trim1,
  iso_model = iso_psbi_50_trim1,
  outcome_var = df_inf_psbi_trim1$inf_psbi,
  title = "PSBI TRIM1 with Bootstrap Spline CI (0.50 Isotonic)"
)
dev.off()

# Save 0.25 plot
png("iso_results/vioplot_psbi_25_trim1.png", width = 800, height = 600)
# 0.25 Isotonic
vioplot_psbi_25_trim1 <- plot_boot_violin(
  x = df_inf_psbi_trim1$hb,
  y = df_inf_psbi_trim1$inf_psbi,
  xlab = "Hb",
  ylab = "PSBI",
  rcs_result = spline_psbi_trim1,
  iso_model = iso_psbi_25_trim1,
  outcome_var = df_inf_psbi_trim1$inf_psbi,
  title = "PSBI TRIM1 with Bootstrap Spline CI (0.25 Isotonic)"
)
dev.off()


# ───────────────────────────────────────
## Save Processed Isotonic Output Objects ----
# ───────────────────────────────────────
save_iso_output <- function(model_obj, label, filename) {
  out <- outdata(model_obj, label)
  save(out, file = file.path("iso_results", filename))
  return(out)
}

out_psbi_25_trim1  <- save_iso_output(iso_psbi_25_trim1,  "PSBI (.25, Trim1)",  "out_psbi_25_trim1.rda")
out_psbi_50_trim1  <- save_iso_output(iso_psbi_50_trim1,  "PSBI (.50, Trim1)",  "out_psbi_50_trim1.rda")




# ───────────────────────────────────────
#TRIMESTER 2 – PSBI Models ----
# ───────────────────────────────────────

# ───────────────────────────────────────
## Fit Spline Models ----
# ───────────────────────────────────────
spline_psbi_trim2 <- knot_fun_boot(df_inf_psbi_trim2, "hb", "inf_psbi")
saveRDS(spline_psbi_trim2, "iso_results/spline_psbi_trim2.rds")

# ───────────────────────────────────────
## Fit Isotonic Models ----
# ───────────────────────────────────────

iso_psbi_25_trim2 <- flexstepreg_glmer_25(df_inf_psbi_trim2$inf_psbi, df_inf_psbi_trim2$hb, df_inf_psbi_trim2$SITE, 
                                               covar2 = NULL, random_effect = NULL, alpha = 0.01)
saveRDS(iso_psbi_25_trim2, "iso_results/iso_psbi_25_trim2.rds")

iso_psbi_50_trim2 <- flexstepreg_glmer(df_inf_psbi_trim2$inf_psbi, df_inf_psbi_trim2$hb, df_inf_psbi_trim2$SITE, 
                                            covar2 = NULL, random_effect = NULL, alpha = 0.01)
saveRDS(iso_psbi_50_trim2, "iso_results/iso_psbi_50_trim2.rds")

# ───────────────────────────────────────
## Generate Plots ----
# ───────────────────────────────────────

# Save 0.50 plot
png("iso_results/vioplot_psbi_50_trim2.png", width = 800, height = 600)
# 0.50 Isotonic
vioplot_psbi_50_trim2 <- plot_boot_violin(
  x = df_inf_psbi_trim2$hb,
  y = df_inf_psbi_trim2$inf_psbi,
  xlab = "Hb",
  ylab = "PSBI",
  rcs_result = spline_psbi_trim2,
  iso_model = iso_psbi_50_trim2,
  outcome_var = df_inf_psbi_trim2$inf_psbi,
  title = "PSBI TRIM2 with Bootstrap Spline CI (0.50 Isotonic)"
)
dev.off()

# Save 0.25 plot
png("iso_results/vioplot_psbi_25_trim2.png", width = 800, height = 600)
# 0.25 Isotonic
vioplot_psbi_25_trim2 <- plot_boot_violin(
  x = df_inf_psbi_trim2$hb,
  y = df_inf_psbi_trim2$inf_psbi,
  xlab = "Hb",
  ylab = "PSBI",
  rcs_result = spline_psbi_trim2,
  iso_model = iso_psbi_25_trim2,
  outcome_var = df_inf_psbi_trim2$inf_psbi,
  title = "PSBI TRIM2 with Bootstrap Spline CI (0.25 Isotonic)"
)
dev.off()


# ───────────────────────────────────────
## Save Processed Isotonic Output Objects ----
# ───────────────────────────────────────
save_iso_output <- function(model_obj, label, filename) {
  out <- outdata(model_obj, label)
  save(out, file = file.path("iso_results", filename))
  return(out)
}

out_psbi_25_trim2  <- save_iso_output(iso_psbi_25_trim2,  "PSBI (.25, Trim2)",  "out_psbi_25_trim2.rda")
out_psbi_50_trim2  <- save_iso_output(iso_psbi_50_trim2,  "PSBI (.50, Trim2)",  "out_psbi_50_trim2.rda")




# ───────────────────────────────────────
#TRIMESTER 3 – PSBI Models ----
# ───────────────────────────────────────

# ───────────────────────────────────────
## Fit Spline Models ----
# ───────────────────────────────────────
spline_psbi_trim3 <- knot_fun_boot(df_inf_psbi_trim3, "hb", "inf_psbi")
saveRDS(spline_psbi_trim3, "iso_results/spline_psbi_trim3.rds")

# ───────────────────────────────────────
## Fit Isotonic Models ----
# ───────────────────────────────────────

iso_psbi_25_trim3 <- flexstepreg_glmer_25(df_inf_psbi_trim3$inf_psbi, df_inf_psbi_trim3$hb, df_inf_psbi_trim3$SITE, 
                                               covar2 = NULL, random_effect = NULL, alpha = 0.01)
saveRDS(iso_psbi_25_trim3, "iso_results/iso_psbi_25_trim3.rds")

iso_psbi_50_trim3 <- flexstepreg_glmer(df_inf_psbi_trim3$inf_psbi, df_inf_psbi_trim3$hb, df_inf_psbi_trim3$SITE, 
                                            covar2 = NULL, random_effect = NULL, alpha = 0.01)
saveRDS(iso_psbi_50_trim3, "iso_results/iso_psbi_50_trim3.rds")

# ───────────────────────────────────────
## Generate Plots ----
# ───────────────────────────────────────

# Save 0.50 plot
png("iso_results/vioplot_psbi_50_trim3.png", width = 800, height = 600)
# 0.50 Isotonic
vioplot_psbi_50_trim3 <- plot_boot_violin(
  x = df_inf_psbi_trim3$hb,
  y = df_inf_psbi_trim3$inf_psbi,
  xlab = "Hb",
  ylab = "PSBI",
  rcs_result = spline_psbi_trim3,
  iso_model = iso_psbi_50_trim3,
  outcome_var = df_inf_psbi_trim3$inf_psbi,
  title = "PSBI TRIM3 with Bootstrap Spline CI (0.50 Isotonic)"
)
dev.off()

# Save 0.25 plot
png("iso_results/vioplot_psbi_25_trim3.png", width = 800, height = 600)
# 0.25 Isotonic
vioplot_psbi_25_trim3 <- plot_boot_violin(
  x = df_inf_psbi_trim3$hb,
  y = df_inf_psbi_trim3$inf_psbi,
  xlab = "Hb",
  ylab = "PSBI",
  rcs_result = spline_psbi_trim3,
  iso_model = iso_psbi_25_trim3,
  outcome_var = df_inf_psbi_trim3$inf_psbi,
  title = "PSBI TRIM3 with Bootstrap Spline CI (0.25 Isotonic)"
)
dev.off()


# ───────────────────────────────────────
## Save Processed Isotonic Output Objects ----
# ───────────────────────────────────────
save_iso_output <- function(model_obj, label, filename) {
  out <- outdata(model_obj, label)
  save(out, file = file.path("iso_results", filename))
  return(out)
}

out_psbi_25_trim3  <- save_iso_output(iso_psbi_25_trim3,  "PSBI (.25, Trim3)",  "out_psbi_25_trim3.rda")
out_psbi_50_trim3  <- save_iso_output(iso_psbi_50_trim3,  "PSBI (.50, Trim3)",  "out_psbi_50_trim3.rda")
