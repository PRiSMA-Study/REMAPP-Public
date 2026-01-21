# ───────────────────────────────────────
# 📁 ReMAPP Aim 2 – PRETERM 37WKS Analysis Pipeline
# 🔧 Modified by: Williams Precious
# 📅 Last updated: 2025-04-18
# 📄 Description: Fits spline and isotonic models for PRETERM 37WKS, generates predictions and plots.
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
#Setup Paths & Load Data ----
# ───────────────────────────────────────
UploadDate <- "2025-10-31"

base_dir <- file.path("D:/Users/williams_pj/Documents/Analysis/ReMAPP/Aim2", UploadDate)
setwd(base_dir)

load("derived_data/df_inf_preterm37.rda")
load("derived_data/df_inf_preterm37_trim1.rda")
load("derived_data/df_inf_preterm37_trim2.rda")
load("derived_data/df_inf_preterm37_trim3.rda")

df_inf_preterm37$hb        <- round(pmin(pmax(df_inf_preterm37$hb, 5), 18), 1)
df_inf_preterm37_trim1$hb  <- round(pmin(pmax(df_inf_preterm37_trim1$hb, 5), 18), 1)
df_inf_preterm37_trim2$hb  <- round(pmin(pmax(df_inf_preterm37_trim2$hb, 5), 18), 1)
df_inf_preterm37_trim3$hb  <- round(pmin(pmax(df_inf_preterm37_trim3$hb, 5), 18), 1)

# Establish the source of the functions
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/spline_binary_boot.R")
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/iso_binary_25_mod.R")
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/iso_binary_50_mod.R")
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/plot_boot.R")
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/outdata_mod.R")

# ───────────────────────────────────────
#ALL TRIMESTERS PRETERM 37WKS ----
# ───────────────────────────────────────

# ───────────────────────────────────────
## Fit Spline Models ----
# ───────────────────────────────────────
spline_preterm37 <- knot_fun_boot(df_inf_preterm37, "hb", "preterm37")
saveRDS(spline_preterm37, "iso_results/spline_preterm37.rds")

# ───────────────────────────────────────
## Fit Isotonic Models ----
# ───────────────────────────────────────

iso_preterm37_25 <- flexstepreg_glmer_25(df_inf_preterm37$preterm37, df_inf_preterm37$hb, df_inf_preterm37$SITE, 
                                         covar2 = NULL, random_effect = NULL, alpha = 0.01)
saveRDS(iso_preterm37_25, "iso_results/iso_preterm37_25.rds")

iso_preterm37_50 <- flexstepreg_glmer(df_inf_preterm37$preterm37, df_inf_preterm37$hb, df_inf_preterm37$SITE, 
                                      covar2 = NULL, random_effect = NULL, alpha = 0.01)
saveRDS(iso_preterm37_50, "iso_results/iso_preterm37_50.rds")

# ───────────────────────────────────────
## Generate Plots ----
# ───────────────────────────────────────

# Save 0.50 plot
png("iso_results/vioplot_preterm37_50.png", width = 800, height = 600)
# 0.50 Isotonic
vioplot_preterm37_50 <- plot_boot_violin(
  x = df_inf_preterm37$hb,
  y = df_inf_preterm37$preterm37,
  xlab = "Hb",
  ylab = "PRETERM 37WKS",
  rcs_result = spline_preterm37,
  iso_model = iso_preterm37_50,
  outcome_var = df_inf_preterm37$preterm37,
  title = "PRETERM 37WKS with Bootstrap Spline CI (0.50 Isotonic)"
)
dev.off()

# Save 0.25 plot
png("iso_results/vioplot_preterm37_25.png", width = 800, height = 600)
# 0.25 Isotonic
vioplot_preterm37_25 <- plot_boot_violin(
  x = df_inf_preterm37$hb,
  y = df_inf_preterm37$preterm37,
  xlab = "Hb",
  ylab = "PRETERM 37WKS",
  rcs_result = spline_preterm37,
  iso_model = iso_preterm37_25,
  outcome_var = df_inf_preterm37$preterm37,
  title = "PRETERM 37WKS with Bootstrap Spline CI (0.25 Isotonic)"
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
out_preterm37_25  <- save_iso_output(iso_preterm37_25,  "PRETERM 37WKS (.25)",  "out_preterm37_25.rda")
out_preterm37_50  <- save_iso_output(iso_preterm37_50,  "PRETERM 37WKS (.50)",  "out_preterm37_50.rda")


# ───────────────────────────────────────
#TRIMESTER 1 – PRETERM 37WKS Models ----
# ───────────────────────────────────────

# ───────────────────────────────────────
## Fit Spline Models ----
# ───────────────────────────────────────
spline_preterm37_trim1 <- knot_fun_boot(df_inf_preterm37_trim1, "hb", "preterm37")
saveRDS(spline_preterm37_trim1, "iso_results/spline_preterm37_trim1.rds")

# ───────────────────────────────────────
## Fit Isotonic Models ----
# ───────────────────────────────────────

iso_preterm37_25_trim1 <- flexstepreg_glmer_25(df_inf_preterm37_trim1$preterm37, df_inf_preterm37_trim1$hb, df_inf_preterm37_trim1$SITE, 
                                               covar2 = NULL, random_effect = NULL, alpha = 0.01)
saveRDS(iso_preterm37_25_trim1, "iso_results/iso_preterm37_25_trim1.rds")

iso_preterm37_50_trim1 <- flexstepreg_glmer(df_inf_preterm37_trim1$preterm37, df_inf_preterm37_trim1$hb, df_inf_preterm37_trim1$SITE, 
                                            covar2 = NULL, random_effect = NULL, alpha = 0.01)
saveRDS(iso_preterm37_50_trim1, "iso_results/iso_preterm37_50_trim1.rds")

# ───────────────────────────────────────
## Generate Plots ----
# ───────────────────────────────────────

# Save 0.50 plot
png("iso_results/vioplot_preterm37_50_trim1.png", width = 800, height = 600)
# 0.50 Isotonic
vioplot_preterm37_50_trim1 <- plot_boot_violin(
  x = df_inf_preterm37_trim1$hb,
  y = df_inf_preterm37_trim1$preterm37,
  xlab = "Hb",
  ylab = "PRETERM 37WKS",
  rcs_result = spline_preterm37_trim1,
  iso_model = iso_preterm37_50_trim1,
  outcome_var = df_inf_preterm37_trim1$preterm37,
  title = "PRETERM 37WKS TRIM1 with Bootstrap Spline CI (0.50 Isotonic)"
)
dev.off()

# Save 0.25 plot
png("iso_results/vioplot_preterm37_25_trim1.png", width = 800, height = 600)
# 0.25 Isotonic
vioplot_preterm37_25_trim1 <- plot_boot_violin(
  x = df_inf_preterm37_trim1$hb,
  y = df_inf_preterm37_trim1$preterm37,
  xlab = "Hb",
  ylab = "PRETERM 37WKS",
  rcs_result = spline_preterm37_trim1,
  iso_model = iso_preterm37_25_trim1,
  outcome_var = df_inf_preterm37_trim1$preterm37,
  title = "PRETERM 37WKS TRIM1 with Bootstrap Spline CI (0.25 Isotonic)"
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

out_preterm37_25_trim1  <- save_iso_output(iso_preterm37_25_trim1,  "PRETERM 37WKS (.25, Trim1)",  "out_preterm37_25_trim1.rda")
out_preterm37_50_trim1  <- save_iso_output(iso_preterm37_50_trim1,  "PRETERM 37WKS (.50, Trim1)",  "out_preterm37_50_trim1.rda")




# ───────────────────────────────────────
#TRIMESTER 2 – PRETERM 37WKS Models ----
# ───────────────────────────────────────

# ───────────────────────────────────────
## Fit Spline Models ----
# ───────────────────────────────────────
spline_preterm37_trim2 <- knot_fun_boot(df_inf_preterm37_trim2, "hb", "preterm37")
saveRDS(spline_preterm37_trim2, "iso_results/spline_preterm37_trim2.rds")

# ───────────────────────────────────────
## Fit Isotonic Models ----
# ───────────────────────────────────────

iso_preterm37_25_trim2 <- flexstepreg_glmer_25(df_inf_preterm37_trim2$preterm37, df_inf_preterm37_trim2$hb, df_inf_preterm37_trim2$SITE, 
                                               covar2 = NULL, random_effect = NULL, alpha = 0.01)
saveRDS(iso_preterm37_25_trim2, "iso_results/iso_preterm37_25_trim2.rds")

iso_preterm37_50_trim2 <- flexstepreg_glmer(df_inf_preterm37_trim2$preterm37, df_inf_preterm37_trim2$hb, df_inf_preterm37_trim2$SITE, 
                                            covar2 = NULL, random_effect = NULL, alpha = 0.01)
saveRDS(iso_preterm37_50_trim2, "iso_results/iso_preterm37_50_trim2.rds")

# ───────────────────────────────────────
## Generate Plots ----
# ───────────────────────────────────────

# Save 0.50 plot
png("iso_results/vioplot_preterm37_50_trim2.png", width = 800, height = 600)
# 0.50 Isotonic
vioplot_preterm37_50_trim2 <- plot_boot_violin(
  x = df_inf_preterm37_trim2$hb,
  y = df_inf_preterm37_trim2$preterm37,
  xlab = "Hb",
  ylab = "PRETERM 37WKS",
  rcs_result = spline_preterm37_trim2,
  iso_model = iso_preterm37_50_trim2,
  outcome_var = df_inf_preterm37_trim2$preterm37,
  title = "PRETERM 37WKS TRIM2 with Bootstrap Spline CI (0.50 Isotonic)"
)
dev.off()

# Save 0.25 plot
png("iso_results/vioplot_preterm37_25_trim2.png", width = 800, height = 600)
# 0.25 Isotonic
vioplot_preterm37_25_trim2 <- plot_boot_violin(
  x = df_inf_preterm37_trim2$hb,
  y = df_inf_preterm37_trim2$preterm37,
  xlab = "Hb",
  ylab = "PRETERM 37WKS",
  rcs_result = spline_preterm37_trim2,
  iso_model = iso_preterm37_25_trim2,
  outcome_var = df_inf_preterm37_trim2$preterm37,
  title = "PRETERM 37WKS TRIM2 with Bootstrap Spline CI (0.25 Isotonic)"
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

out_preterm37_25_trim2  <- save_iso_output(iso_preterm37_25_trim2,  "PRETERM 37WKS (.25, Trim2)",  "out_preterm37_25_trim2.rda")
out_preterm37_50_trim2  <- save_iso_output(iso_preterm37_50_trim2,  "PRETERM 37WKS (.50, Trim2)",  "out_preterm37_50_trim2.rda")




# ───────────────────────────────────────
#TRIMESTER 3 – PRETERM 37WKS Models ----
# ───────────────────────────────────────

# ───────────────────────────────────────
## Fit Spline Models ----
# ───────────────────────────────────────
spline_preterm37_trim3 <- knot_fun_boot(df_inf_preterm37_trim3, "hb", "preterm37")
saveRDS(spline_preterm37_trim3, "iso_results/spline_preterm37_trim3.rds")

# ───────────────────────────────────────
## Fit Isotonic Models ----
# ───────────────────────────────────────

iso_preterm37_25_trim3 <- flexstepreg_glmer_25(df_inf_preterm37_trim3$preterm37, df_inf_preterm37_trim3$hb, df_inf_preterm37_trim3$SITE, 
                                               covar2 = NULL, random_effect = NULL, alpha = 0.01)
saveRDS(iso_preterm37_25_trim3, "iso_results/iso_preterm37_25_trim3.rds")

iso_preterm37_50_trim3 <- flexstepreg_glmer(df_inf_preterm37_trim3$preterm37, df_inf_preterm37_trim3$hb, df_inf_preterm37_trim3$SITE, 
                                            covar2 = NULL, random_effect = NULL, alpha = 0.01)
saveRDS(iso_preterm37_50_trim3, "iso_results/iso_preterm37_50_trim3.rds")

# ───────────────────────────────────────
## Generate Plots ----
# ───────────────────────────────────────

# Save 0.50 plot
png("iso_results/vioplot_preterm37_50_trim3.png", width = 800, height = 600)
# 0.50 Isotonic
vioplot_preterm37_50_trim3 <- plot_boot_violin(
  x = df_inf_preterm37_trim3$hb,
  y = df_inf_preterm37_trim3$preterm37,
  xlab = "Hb",
  ylab = "PRETERM 37WKS",
  rcs_result = spline_preterm37_trim3,
  iso_model = iso_preterm37_50_trim3,
  outcome_var = df_inf_preterm37_trim3$preterm37,
  title = "PRETERM 37WKS TRIM3 with Bootstrap Spline CI (0.50 Isotonic)"
)
dev.off()

# Save 0.25 plot
png("iso_results/vioplot_preterm37_25_trim3.png", width = 800, height = 600)
# 0.25 Isotonic
vioplot_preterm37_25_trim3 <- plot_boot_violin(
  x = df_inf_preterm37_trim3$hb,
  y = df_inf_preterm37_trim3$preterm37,
  xlab = "Hb",
  ylab = "PRETERM 37WKS",
  rcs_result = spline_preterm37,
  iso_model = iso_preterm37_25_trim3,
  outcome_var = df_inf_preterm37_trim3$preterm37,
  title = "PRETERM 37WKS TRIM3 with Bootstrap Spline CI (0.25 Isotonic)"
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

out_preterm37_25_trim3  <- save_iso_output(iso_preterm37_25_trim3,  "PRETERM 37WKS (.25, Trim3)",  "out_preterm37_25_trim3.rda")
out_preterm37_50_trim3  <- save_iso_output(iso_preterm37_50_trim3,  "PRETERM 37WKS (.50, Trim3)",  "out_preterm37_50_trim3.rda")
