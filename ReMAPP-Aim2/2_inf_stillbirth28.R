# ───────────────────────────────────────
# 📁 ReMAPP Aim 2 – STILLBIRTH 28WKS Analysis Pipeline
# 🔧 Modified by: Williams Precious
# 📅 Last updated: 2025-04-18
# 📄 Description: Fits spline and isotonic models for STILLBIRTH 28WKS, generates predictions and plots.
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

load("derived_data/df_inf_stillbirth28.rda")
load("derived_data/df_inf_stillbirth28_trim1.rda")
load("derived_data/df_inf_stillbirth28_trim2.rda")
#load("derived_data/df_inf_stillbirth28_trim3.rda")

df_inf_stillbirth28$hb        <- round(pmin(pmax(df_inf_stillbirth28$hb, 5), 18), 1)
df_inf_stillbirth28_trim1$hb  <- round(pmin(pmax(df_inf_stillbirth28_trim1$hb, 5), 18), 1)
df_inf_stillbirth28_trim2$hb  <- round(pmin(pmax(df_inf_stillbirth28_trim2$hb, 5), 18), 1)
#df_inf_stillbirth28_trim3$hb  <- round(pmin(pmax(df_inf_stillbirth28_trim3$hb, 5), 18), 1)

# Establish the source of the functions
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/spline_binary_boot.R")
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/iso_binary_25_mod.R")
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/iso_binary_50_mod.R")
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/plot_boot.R")
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/outdata_mod.R")

# ───────────────────────────────────────
#ALL TRIMESTERS STILLBIRTH 28WKS ----
# ───────────────────────────────────────

# ───────────────────────────────────────
## Fit Spline Models ----
# ───────────────────────────────────────
spline_stillbirth28 <- knot_fun_boot(df_inf_stillbirth28, "hb", "inf_stillbirth28")
saveRDS(spline_stillbirth28, "iso_results/spline_stillbirth28.rds")

# ───────────────────────────────────────
## Fit Isotonic Models ----
# ───────────────────────────────────────

iso_stillbirth28_25 <- flexstepreg_glmer_25(df_inf_stillbirth28$inf_stillbirth28, df_inf_stillbirth28$hb, df_inf_stillbirth28$SITE, 
                                            covar2 = NULL, random_effect = NULL, alpha = 0.01)
saveRDS(iso_stillbirth28_25, "iso_results/iso_stillbirth28_25.rds")

iso_stillbirth28_50 <- flexstepreg_glmer(df_inf_stillbirth28$inf_stillbirth28, df_inf_stillbirth28$hb, df_inf_stillbirth28$SITE, 
                                         covar2 = NULL, random_effect = NULL, alpha = 0.01)
saveRDS(iso_stillbirth28_50, "iso_results/iso_stillbirth28_50.rds")

# ───────────────────────────────────────
## Generate Plots ----
# ───────────────────────────────────────

# Save 0.50 plot
png("iso_results/vioplot_stillbirth28_50.png", width = 800, height = 600)
# 0.50 Isotonic
vioplot_stillbirth28_50 <- plot_boot_violin(
  x = df_inf_stillbirth28$hb,
  y = df_inf_stillbirth28$inf_stillbirth28,
  xlab = "Hb",
  ylab = "STILLBIRTH 28WKS",
  rcs_result = spline_stillbirth28,
  iso_model = iso_stillbirth28_50,
  outcome_var = df_inf_stillbirth28$inf_stillbirth28,
  title = "STILLBIRTH 28WKS with Bootstrap Spline CI (0.50 Isotonic)"
)
dev.off()

# Save 0.25 plot
png("iso_results/vioplot_stillbirth28_25.png", width = 800, height = 600)
# 0.25 Isotonic
vioplot_stillbirth28_25 <- plot_boot_violin(
  x = df_inf_stillbirth28$hb,
  y = df_inf_stillbirth28$inf_stillbirth28,
  xlab = "Hb",
  ylab = "STILLBIRTH 28WKS",
  rcs_result = spline_stillbirth28,
  iso_model = iso_stillbirth28_25,
  outcome_var = df_inf_stillbirth28$inf_stillbirth28,
  title = "STILLBIRTH 28WKS with Bootstrap Spline CI (0.25 Isotonic)"
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
out_stillbirth28_25  <- save_iso_output(iso_stillbirth28_25,  "STILLBIRTH 28WKS (.25)",  "out_stillbirth28_25.rda")
out_stillbirth28_50  <- save_iso_output(iso_stillbirth28_50,  "STILLBIRTH 28WKS (.50)",  "out_stillbirth28_50.rda")


# ───────────────────────────────────────
#TRIMESTER 1 – STILLBIRTH 28WKS Models ----
# ───────────────────────────────────────

# ───────────────────────────────────────
## Fit Spline Models ----
# ───────────────────────────────────────
spline_stillbirth28_trim1 <- knot_fun_boot(df_inf_stillbirth28_trim1, "hb", "inf_stillbirth28")
saveRDS(spline_stillbirth28_trim1, "iso_results/spline_stillbirth28_trim1.rds")

# ───────────────────────────────────────
## Fit Isotonic Models ----
# ───────────────────────────────────────

iso_stillbirth28_25_trim1 <- flexstepreg_glmer_25(df_inf_stillbirth28_trim1$inf_stillbirth28, df_inf_stillbirth28_trim1$hb, df_inf_stillbirth28_trim1$SITE, 
                                                  covar2 = NULL, random_effect = NULL, alpha = 0.01)
saveRDS(iso_stillbirth28_25_trim1, "iso_results/iso_stillbirth28_25_trim1.rds")

iso_stillbirth28_50_trim1 <- flexstepreg_glmer(df_inf_stillbirth28_trim1$inf_stillbirth28, df_inf_stillbirth28_trim1$hb, df_inf_stillbirth28_trim1$SITE, 
                                               covar2 = NULL, random_effect = NULL, alpha = 0.01)
saveRDS(iso_stillbirth28_50_trim1, "iso_results/iso_stillbirth28_50_trim1.rds")

# ───────────────────────────────────────
## Generate Plots ----
# ───────────────────────────────────────

# Save 0.50 plot
png("iso_results/vioplot_stillbirth28_50_trim1.png", width = 800, height = 600)
# 0.50 Isotonic
vioplot_stillbirth28_50_trim1 <- plot_boot_violin(
  x = df_inf_stillbirth28_trim1$hb,
  y = df_inf_stillbirth28_trim1$inf_stillbirth28,
  xlab = "Hb",
  ylab = "STILLBIRTH 28WKS",
  rcs_result = spline_stillbirth28_trim1,
  iso_model = iso_stillbirth28_50_trim1,
  outcome_var = df_inf_stillbirth28_trim1$inf_stillbirth28,
  title = "STILLBIRTH 28WKS TRIM1 with Bootstrap Spline CI (0.50 Isotonic)"
)
dev.off()

# Save 0.25 plot
png("iso_results/vioplot_stillbirth28_25_trim1.png", width = 800, height = 600)
# 0.25 Isotonic
vioplot_stillbirth28_25_trim1 <- plot_boot_violin(
  x = df_inf_stillbirth28_trim1$hb,
  y = df_inf_stillbirth28_trim1$inf_stillbirth28,
  xlab = "Hb",
  ylab = "STILLBIRTH 28WKS",
  rcs_result = spline_stillbirth28_trim1,
  iso_model = iso_stillbirth28_25_trim1,
  outcome_var = df_inf_stillbirth28_trim1$inf_stillbirth28,
  title = "STILLBIRTH 28WKS TRIM1 with Bootstrap Spline CI (0.25 Isotonic)"
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

out_stillbirth28_25_trim1  <- save_iso_output(iso_stillbirth28_25_trim1,  "STILLBIRTH 28WKS (.25, Trim1)",  "out_stillbirth28_25_trim1.rda")
out_stillbirth28_50_trim1  <- save_iso_output(iso_stillbirth28_50_trim1,  "STILLBIRTH 28WKS (.50, Trim1)",  "out_stillbirth28_50_trim1.rda")




# ───────────────────────────────────────
#TRIMESTER 2 – STILLBIRTH 28WKS Models ----
# ───────────────────────────────────────

# ───────────────────────────────────────
## Fit Spline Models ----
# ───────────────────────────────────────
spline_stillbirth28_trim2 <- knot_fun_boot(df_inf_stillbirth28_trim2, "hb", "inf_stillbirth28")
saveRDS(spline_stillbirth28_trim2, "iso_results/spline_stillbirth28_trim2.rds")

# ───────────────────────────────────────
## Fit Isotonic Models ----
# ───────────────────────────────────────

iso_stillbirth28_25_trim2 <- flexstepreg_glmer_25(df_inf_stillbirth28_trim2$inf_stillbirth28, df_inf_stillbirth28_trim2$hb, df_inf_stillbirth28_trim2$SITE, 
                                                  covar2 = NULL, random_effect = NULL, alpha = 0.01)
saveRDS(iso_stillbirth28_25_trim2, "iso_results/iso_stillbirth28_25_trim2.rds")

iso_stillbirth28_50_trim2 <- flexstepreg_glmer(df_inf_stillbirth28_trim2$inf_stillbirth28, df_inf_stillbirth28_trim2$hb, df_inf_stillbirth28_trim2$SITE, 
                                               covar2 = NULL, random_effect = NULL, alpha = 0.01)
saveRDS(iso_stillbirth28_50_trim2, "iso_results/iso_stillbirth28_50_trim2.rds")

# ───────────────────────────────────────
## Generate Plots ----
# ───────────────────────────────────────

# Save 0.50 plot
png("iso_results/vioplot_stillbirth28_50_trim2.png", width = 800, height = 600)
# 0.50 Isotonic
vioplot_stillbirth28_50_trim2 <- plot_boot_violin(
  x = df_inf_stillbirth28_trim2$hb,
  y = df_inf_stillbirth28_trim2$inf_stillbirth28,
  xlab = "Hb",
  ylab = "STILLBIRTH 28WKS",
  rcs_result = spline_stillbirth28_trim2,
  iso_model = iso_stillbirth28_50_trim2,
  outcome_var = df_inf_stillbirth28_trim2$inf_stillbirth28,
  title = "STILLBIRTH 28WKS TRIM2 with Bootstrap Spline CI (0.50 Isotonic)"
)
dev.off()

# Save 0.25 plot
png("iso_results/vioplot_stillbirth28_25_trim2.png", width = 800, height = 600)
# 0.25 Isotonic
vioplot_stillbirth28_25_trim2 <- plot_boot_violin(
  x = df_inf_stillbirth28_trim2$hb,
  y = df_inf_stillbirth28_trim2$inf_stillbirth28,
  xlab = "Hb",
  ylab = "STILLBIRTH 28WKS",
  rcs_result = spline_stillbirth28_trim2,
  iso_model = iso_stillbirth28_25_trim2,
  outcome_var = df_inf_stillbirth28_trim2$inf_stillbirth28,
  title = "STILLBIRTH 28WKS TRIM2 with Bootstrap Spline CI (0.25 Isotonic)"
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

out_stillbirth28_25_trim2  <- save_iso_output(iso_stillbirth28_25_trim2,  "STILLBIRTH 28WKS (.25, Trim2)",  "out_stillbirth28_25_trim2.rda")
out_stillbirth28_50_trim2  <- save_iso_output(iso_stillbirth28_50_trim2,  "STILLBIRTH 28WKS (.50, Trim2)",  "out_stillbirth28_50_trim2.rda")


