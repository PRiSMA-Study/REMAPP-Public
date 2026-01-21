# ───────────────────────────────────────
# 📁 ReMAPP Aim 2 – PRECLAMP Analysis Pipeline
# 🔧 Modified by: Williams Precious
# 📅 Last updated: 2025-04-18
# 📄 Description: Fits spline and isotonic models for PRECLAMP, generates predictions and plots.
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
UploadDate <- "2025-10-31"

base_dir <- file.path("D:/Users/williams_pj/Documents/Analysis/ReMAPP/Aim2", UploadDate)
setwd(base_dir)

load("derived_data/df_mat_preclamp.rda")
load("derived_data/df_mat_preclamp_trim1.rda")
load("derived_data/df_mat_preclamp_trim2.rda")
load("derived_data/df_mat_preclamp_trim3.rda")


df_mat_preclamp$hb        <- round(pmin(pmax(df_mat_preclamp$hb, 5), 18), 1)
df_mat_preclamp_trim1$hb  <- round(pmin(pmax(df_mat_preclamp_trim1$hb, 5), 18), 1)
df_mat_preclamp_trim2$hb  <- round(pmin(pmax(df_mat_preclamp_trim2$hb, 5), 18), 1)
df_mat_preclamp_trim3$hb  <- round(pmin(pmax(df_mat_preclamp_trim3$hb, 5), 18), 1)

# Establish the source of the functions
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/spline_binary_boot.R")
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/iso_binary_25_mod.R")
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/iso_binary_50_mod.R")
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/plot_boot.R")
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/outdata_mod.R")

# ───────────────────────────────────────
#ALL TRIMESTERS PRECLAMP ----
# ───────────────────────────────────────

# ───────────────────────────────────────
## Fit Spline Models ----
# ───────────────────────────────────────
spline_preclamp <- knot_fun_boot(df_mat_preclamp, "hb", "preclamp")
saveRDS(spline_preclamp, "iso_results/spline_preclamp.rds")

# ───────────────────────────────────────
## Fit Isotonic Models ----
# ───────────────────────────────────────

iso_preclamp_25 <- flexstepreg_glmer_25(df_mat_preclamp$preclamp, df_mat_preclamp$hb, df_mat_preclamp$SITE, 
                                     covar2 = NULL, random_effect = NULL, alpha = 0.01)
saveRDS(iso_preclamp_25, "iso_results/iso_preclamp_25.rds")

iso_preclamp_50 <- flexstepreg_glmer(df_mat_preclamp$preclamp, df_mat_preclamp$hb, df_mat_preclamp$SITE, 
                                  covar2 = NULL, random_effect = NULL, alpha = 0.01)
saveRDS(iso_preclamp_50, "iso_results/iso_preclamp_50.rds")

# ───────────────────────────────────────
## Generate Plots ----
# ───────────────────────────────────────
# Save 0.50 plot
png("iso_results/vioplot_preclamp_50.png", width = 800, height = 600)
# 0.50 Isotonic
vioplot_preclamp_50 <- plot_boot_violin(
  x = df_mat_preclamp$hb,
  y = df_mat_preclamp$preclamp,
  xlab = "Hb",
  ylab = "PRECLAMP",
  rcs_result = spline_preclamp,
  iso_model = iso_preclamp_50,
  outcome_var = df_mat_preclamp$preclamp,
  title = "PRECLAMP with Bootstrap Spline CI (0.50 Isotonic)"
)
dev.off()

# Save 0.25 plot
png("iso_results/vioplot_preclamp_25.png", width = 800, height = 600)
# 0.25 Isotonic
vioplot_preclamp_25 <- plot_boot_violin(
  x = df_mat_preclamp$hb,
  y = df_mat_preclamp$preclamp,
  xlab = "Hb",
  ylab = "PRECLAMP",
  rcs_result = spline_preclamp,
  iso_model = iso_preclamp_25,
  outcome_var = df_mat_preclamp$preclamp,
  title = "PRECLAMP with Bootstrap Spline CI (0.25 Isotonic)"
)
dev.off()


# ───────────────────────────────────────
## Save Processed Isotonic Output Objects ----
# ───────────────────────────────────────

# Helper to format and save isotonic model outputs
save_iso_output <- function(model_obj, label, filename) {
  out <- outdata(model_obj, label)
  save(out, file = file.path("iso_results", filename))
  return(out)
}

# Save outputs
out_preclamp_25  <- save_iso_output(iso_preclamp_25,  "PRECLAMP (.25)",  "out_preclamp_25.rda")
out_preclamp_50  <- save_iso_output(iso_preclamp_50,  "PRECLAMP (.50)",  "out_preclamp_50.rda")


# ───────────────────────────────────────
#TRIMESTER 1 – PRECLAMP Models ----
# ───────────────────────────────────────

# ───────────────────────────────────────
## Fit Spline Models ----
# ───────────────────────────────────────
spline_preclamp_trim1 <- knot_fun_boot(df_mat_preclamp_trim1, "hb", "preclamp")
saveRDS(spline_preclamp_trim1, "iso_results/spline_preclamp_trim1.rds")

# ───────────────────────────────────────
## Fit Isotonic Models ----
# ───────────────────────────────────────

iso_preclamp_25_trim1 <- flexstepreg_glmer_25(df_mat_preclamp_trim1$preclamp, df_mat_preclamp_trim1$hb, df_mat_preclamp_trim1$SITE, 
                                           covar2 = NULL, random_effect = NULL, alpha = 0.01)
saveRDS(iso_preclamp_25_trim1, "iso_results/iso_preclamp_25_trim1.rds")

iso_preclamp_50_trim1 <- flexstepreg_glmer(df_mat_preclamp_trim1$preclamp, df_mat_preclamp_trim1$hb, df_mat_preclamp_trim1$SITE, 
                                        covar2 = NULL, random_effect = NULL, alpha = 0.01)
saveRDS(iso_preclamp_50_trim1, "iso_results/iso_preclamp_50_trim1.rds")

# ───────────────────────────────────────
## Generate Plots ----
# ───────────────────────────────────────

# Save 0.50 plot
png("iso_results/vioplot_preclamp_50_trim1.png", width = 800, height = 600)
# 0.50 Isotonic
vioplot_preclamp_50_trim1 <- plot_boot_violin(
  x = df_mat_preclamp_trim1$hb,
  y = df_mat_preclamp_trim1$preclamp,
  xlab = "Hb",
  ylab = "PRECLAMP",
  rcs_result = spline_preclamp,
  iso_model = iso_preclamp_50_trim1,
  outcome_var = df_mat_preclamp_trim1$preclamp,
  title = "PRECLAMP TRIM1 with Bootstrap Spline CI (0.50 Isotonic)"
)
dev.off()

# Save 0.25 plot
png("iso_results/vioplot_preclamp_25_trim1.png", width = 800, height = 600)
# 0.25 Isotonic
vioplot_preclamp_25_trim1 <- plot_boot_violin(
  x = df_mat_preclamp_trim1$hb,
  y = df_mat_preclamp_trim1$preclamp,
  xlab = "Hb",
  ylab = "PRECLAMP",
  rcs_result = spline_preclamp,
  iso_model = iso_preclamp_25_trim1,
  outcome_var = df_mat_preclamp_trim1$preclamp,
  title = "PRECLAMP TRIM1 with Bootstrap Spline CI (0.25 Isotonic)"
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

out_preclamp_25_trim1  <- save_iso_output(iso_preclamp_25_trim1,  "PRECLAMP (.25, Trim1)",  "out_preclamp_25_trim1.rda")
out_preclamp_50_trim1  <- save_iso_output(iso_preclamp_50_trim1,  "PRECLAMP (.50, Trim1)",  "out_preclamp_50_trim1.rda")



# ───────────────────────────────────────
#TRIMESTER 2 – PRECLAMP Models ----
# ───────────────────────────────────────

# ───────────────────────────────────────
## Fit Spline Models ----
# ───────────────────────────────────────
spline_preclamp_trim2 <- knot_fun_boot(df_mat_preclamp_trim2, "hb", "preclamp")
saveRDS(spline_preclamp_trim2, "iso_results/spline_preclamp_trim2.rds")

# ───────────────────────────────────────
## Fit Isotonic Models ----
# ───────────────────────────────────────

iso_preclamp_25_trim2 <- flexstepreg_glmer_25(df_mat_preclamp_trim2$preclamp, df_mat_preclamp_trim2$hb, df_mat_preclamp_trim2$SITE, 
                                           covar2 = NULL, random_effect = NULL, alpha = 0.01)
saveRDS(iso_preclamp_25_trim2, "iso_results/iso_preclamp_25_trim2.rds")

iso_preclamp_50_trim2 <- flexstepreg_glmer(df_mat_preclamp_trim2$preclamp, df_mat_preclamp_trim2$hb, df_mat_preclamp_trim2$SITE, 
                                        covar2 = NULL, random_effect = NULL, alpha = 0.01)
saveRDS(iso_preclamp_50_trim2, "iso_results/iso_preclamp_50_trim2.rds")

# ───────────────────────────────────────
## Generate Plots ----
# ───────────────────────────────────────
# Save 0.50 plot
png("iso_results/vioplot_preclamp_50_trim2.png", width = 800, height = 600)
# 0.50 Isotonic
vioplot_preclamp_50_trim2 <- plot_boot_violin(
  x = df_mat_preclamp_trim2$hb,
  y = df_mat_preclamp_trim2$preclamp,
  xlab = "Hb",
  ylab = "PRECLAMP",
  rcs_result = spline_preclamp,
  iso_model = iso_preclamp_50_trim2,
  outcome_var = df_mat_preclamp_trim2$preclamp,
  title = "PRECLAMP TRIM2 with Bootstrap Spline CI (0.50 Isotonic)"
)
dev.off()

# Save 0.25 plot
png("iso_results/vioplot_preclamp_25_trim2.png", width = 800, height = 600)
# 0.25 Isotonic
vioplot_preclamp_25_trim2 <- plot_boot_violin(
  x = df_mat_preclamp_trim2$hb,
  y = df_mat_preclamp_trim2$preclamp,
  xlab = "Hb",
  ylab = "PRECLAMP",
  rcs_result = spline_preclamp,
  iso_model = iso_preclamp_25_trim2,
  outcome_var = df_mat_preclamp_trim2$preclamp,
  title = "PRECLAMP TRIM2 with Bootstrap Spline CI (0.25 Isotonic)"
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

out_preclamp_25_trim2  <- save_iso_output(iso_preclamp_25_trim2,  "PRECLAMP (.25, Trim2)",  "out_preclamp_25_trim2.rda")
out_preclamp_50_trim2  <- save_iso_output(iso_preclamp_50_trim2,  "PRECLAMP (.50, Trim2)",  "out_preclamp_50_trim2.rda")




# ───────────────────────────────────────
#TRIMESTER 3 – PRECLAMP Models ----
# ───────────────────────────────────────

# ───────────────────────────────────────
## Fit Spline Models ----
# ───────────────────────────────────────
spline_preclamp_trim3 <- knot_fun_boot(df_mat_preclamp_trim3, "hb", "preclamp")
saveRDS(spline_preclamp_trim3, "iso_results/spline_preclamp_trim3.rds")

# ───────────────────────────────────────
## Fit Isotonic Models ----
# ───────────────────────────────────────

iso_preclamp_25_trim3 <- flexstepreg_glmer_25(df_mat_preclamp_trim3$preclamp, df_mat_preclamp_trim3$hb, df_mat_preclamp_trim3$SITE, 
                                           covar2 = NULL, random_effect = NULL, alpha = 0.01)
saveRDS(iso_preclamp_25_trim3, "iso_results/iso_preclamp_25_trim3.rds")

iso_preclamp_50_trim3 <- flexstepreg_glmer(df_mat_preclamp_trim3$preclamp, df_mat_preclamp_trim3$hb, df_mat_preclamp_trim3$SITE, 
                                        covar2 = NULL, random_effect = NULL, alpha = 0.01)
saveRDS(iso_preclamp_50_trim3, "iso_results/iso_preclamp_50_trim3.rds")

# ───────────────────────────────────────
## Generate Plots ----
# ───────────────────────────────────────

# Save 0.50 plot
png("iso_results/vioplot_preclamp_50_trim3.png", width = 800, height = 600)
# 0.50 Isotonic
vioplot_preclamp_50_trim3 <- plot_boot_violin(
  x = df_mat_preclamp_trim3$hb,
  y = df_mat_preclamp_trim3$preclamp,
  xlab = "Hb",
  ylab = "PRECLAMP",
  rcs_result = spline_preclamp,
  iso_model = iso_preclamp_50_trim3,
  outcome_var = df_mat_preclamp_trim3$preclamp,
  title = "PRECLAMP TRIM3 with Bootstrap Spline CI (0.50 Isotonic)"
)
dev.off()

# Save 0.25 plot
png("iso_results/vioplot_preclamp_25_trim3.png", width = 800, height = 600)
# 0.25 Isotonic
vioplot_preclamp_25_trim3 <- plot_boot_violin(
  x = df_mat_preclamp_trim3$hb,
  y = df_mat_preclamp_trim3$preclamp,
  xlab = "Hb",
  ylab = "PRECLAMP",
  rcs_result = spline_preclamp,
  iso_model = iso_preclamp_25_trim3,
  outcome_var = df_mat_preclamp_trim3$preclamp,
  title = "PRECLAMP TRIM3 with Bootstrap Spline CI (0.25 Isotonic)"
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

out_preclamp_25_trim3  <- save_iso_output(iso_preclamp_25_trim3,  "PRECLAMP (.25, Trim3)",  "out_preclamp_25_trim3.rda")
out_preclamp_50_trim3  <- save_iso_output(iso_preclamp_50_trim3,  "PRECLAMP (.50, Trim3)",  "out_preclamp_50_trim3.rda")
