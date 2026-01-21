# ───────────────────────────────────────
# 📁 ReMAPP Aim 2 – PPA PNC6 Analysis Pipeline
# 🔧 Modified by: Williams Precious
# 📄 Description: Fits spline and isotonic models for PPA PNC6, generates predictions and plots.
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

load("derived_data/df_mat_ppa_pnc6.rda")
load("derived_data/df_mat_ppa_pnc6_trim1.rda")
load("derived_data/df_mat_ppa_pnc6_trim2.rda")
load("derived_data/df_mat_ppa_pnc6_trim3.rda")

df_mat_ppa_pnc6$hb        <- round(pmin(pmax(df_mat_ppa_pnc6$hb, 5), 18), 1)
df_mat_ppa_pnc6_trim1$hb  <- round(pmin(pmax(df_mat_ppa_pnc6_trim1$hb, 5), 18), 1)
df_mat_ppa_pnc6_trim2$hb  <- round(pmin(pmax(df_mat_ppa_pnc6_trim2$hb, 5), 18), 1)
df_mat_ppa_pnc6_trim3$hb  <- round(pmin(pmax(df_mat_ppa_pnc6_trim3$hb, 5), 18), 1)

# Establish the source of the functions
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/spline_binary_boot.R")
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/iso_binary_25_mod.R")
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/iso_binary_50_mod.R")
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/plot_boot.R")
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/outdata_mod.R")

# ───────────────────────────────────────
#ALL TRIMESTERS PPA PNC6 ----
# ───────────────────────────────────────

# ───────────────────────────────────────
## Fit Spline Models ----
# ───────────────────────────────────────
spline_ppa_pnc6 <- knot_fun_boot(df_mat_ppa_pnc6, "hb", "ppa_pnc6")
saveRDS(spline_ppa_pnc6, "iso_results/spline_ppa_pnc6.rds")

# ───────────────────────────────────────
## Fit Isotonic Models ----
# ───────────────────────────────────────

iso_ppa_pnc6_25 <- flexstepreg_glmer_25(df_mat_ppa_pnc6$ppa_pnc6, df_mat_ppa_pnc6$hb, df_mat_ppa_pnc6$SITE, 
                                 covar2 = NULL, random_effect = NULL, alpha = 0.01)
saveRDS(iso_ppa_pnc6_25, "iso_results/iso_ppa_pnc6_25.rds")

iso_ppa_pnc6_50 <- flexstepreg_glmer(df_mat_ppa_pnc6$ppa_pnc6, df_mat_ppa_pnc6$hb, df_mat_ppa_pnc6$SITE, 
                                 covar2 = NULL, random_effect = NULL, alpha = 0.01)
saveRDS(iso_ppa_pnc6_50, "iso_results/iso_ppa_pnc6_50.rds")

# ───────────────────────────────────────
## Generate Plots ----
# ───────────────────────────────────────

# Save 0.50 plot
png("iso_results/vioplot_ppa_pnc6_50.png", width = 800, height = 600)
# 0.50 Isotonic
vioplot_ppa_pnc6_50 <- plot_boot_violin(
  x = df_mat_ppa_pnc6$hb,
  y = df_mat_ppa_pnc6$ppa_pnc6,
  xlab = "Hb",
  ylab = "PPA PNC6",
  rcs_result = spline_ppa_pnc6,
  iso_model = iso_ppa_pnc6_50,
  outcome_var = df_mat_ppa_pnc6$ppa_pnc6,
  title = "PPA PNC6 with Bootstrap Spline CI (0.50 Isotonic)"
)
dev.off()

# Save 0.25 plot
png("iso_results/vioplot_ppa_pnc6_25.png", width = 800, height = 600)
# 0.25 Isotonic
vioplot_ppa_pnc6_25 <- plot_boot_violin(
  x = df_mat_ppa_pnc6$hb,
  y = df_mat_ppa_pnc6$ppa_pnc6,
  xlab = "Hb",
  ylab = "PPA PNC6",
  rcs_result = spline_ppa_pnc6,
  iso_model = iso_ppa_pnc6_25,
  outcome_var = df_mat_ppa_pnc6$ppa_pnc6,
  title = "PPA PNC6 with Bootstrap Spline CI (0.25 Isotonic)"
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
out_ppa_pnc6_25  <- save_iso_output(iso_ppa_pnc6_25,  "PPA PNC6 (.25)",  "out_ppa_pnc6_25.rda")
out_ppa_pnc6_50  <- save_iso_output(iso_ppa_pnc6_50,  "PPA PNC6 (.50)",  "out_ppa_pnc6_50.rda")


# ───────────────────────────────────────
#TRIMESTER 1 – PPA PNC6 Models ----
# ───────────────────────────────────────

# ───────────────────────────────────────
## Fit Spline Models ----
# ───────────────────────────────────────
spline_ppa_pnc6_trim1 <- knot_fun_boot(df_mat_ppa_pnc6_trim1, "hb", "ppa_pnc6")
saveRDS(spline_ppa_pnc6_trim1, "iso_results/spline_ppa_pnc6_trim1.rds")

# ───────────────────────────────────────
## Fit Isotonic Models ----
# ───────────────────────────────────────

iso_ppa_pnc6_25_trim1 <- flexstepreg_glmer_25(df_mat_ppa_pnc6_trim1$ppa_pnc6, df_mat_ppa_pnc6_trim1$hb, df_mat_ppa_pnc6_trim1$SITE, 
                                        covar2 = NULL, random_effect = NULL, alpha = 0.01)
saveRDS(iso_ppa_pnc6_25_trim1, "iso_results/iso_ppa_pnc6_25_trim1.rds")

iso_ppa_pnc6_50_trim1 <- flexstepreg_glmer(df_mat_ppa_pnc6_trim1$ppa_pnc6, df_mat_ppa_pnc6_trim1$hb, df_mat_ppa_pnc6_trim1$SITE, 
                                     covar2 = NULL, random_effect = NULL, alpha = 0.01)
saveRDS(iso_ppa_pnc6_50_trim1, "iso_results/iso_ppa_pnc6_50_trim1.rds")

# ───────────────────────────────────────
## Generate Plots ----
# ───────────────────────────────────────

# Save 0.50 plot
png("iso_results/vioplot_ppa_pnc6_50_trim1.png", width = 800, height = 600)
# 0.50 Isotonic
vioplot_ppa_pnc6_50_trim1 <- plot_boot_violin(
  x = df_mat_ppa_pnc6_trim1$hb,
  y = df_mat_ppa_pnc6_trim1$ppa_pnc6,
  xlab = "Hb",
  ylab = "PPA PNC6",
  rcs_result = spline_ppa_pnc6_trim1,
  iso_model = iso_ppa_pnc6_50_trim1,
  outcome_var = df_mat_ppa_pnc6_trim1$ppa_pnc6,
  title = "PPA PNC6 TRIM1 with Bootstrap Spline CI (0.50 Isotonic)"
)
dev.off()

# Save 0.25 plot
png("iso_results/vioplot_ppa_pnc6_25_trim1.png", width = 800, height = 600)
# 0.25 Isotonic
vioplot_ppa_pnc6_25_trim1 <- plot_boot_violin(
  x = df_mat_ppa_pnc6_trim1$hb,
  y = df_mat_ppa_pnc6_trim1$ppa_pnc6,
  xlab = "Hb",
  ylab = "PPA PNC6",
  rcs_result = spline_ppa_pnc6_trim1,
  iso_model = iso_ppa_pnc6_25_trim1,
  outcome_var = df_mat_ppa_pnc6_trim1$ppa_pnc6,
  title = "PPA PNC6 TRIM1 with Bootstrap Spline CI (0.25 Isotonic)"
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

out_ppa_pnc6_25_trim1  <- save_iso_output(iso_ppa_pnc6_25_trim1,  "PPA PNC6 (.25, Trim1)",  "out_ppa_pnc6_25_trim1.rda")
out_ppa_pnc6_50_trim1  <- save_iso_output(iso_ppa_pnc6_50_trim1,  "PPA PNC6 (.50, Trim1)",  "out_ppa_pnc6_50_trim1.rda")




# ───────────────────────────────────────
#TRIMESTER 2 – PPA PNC6 Models ----
# ───────────────────────────────────────

# ───────────────────────────────────────
## Fit Spline Models ----
# ───────────────────────────────────────
spline_ppa_pnc6_trim2 <- knot_fun_boot(df_mat_ppa_pnc6_trim2, "hb", "ppa_pnc6")
saveRDS(spline_ppa_pnc6_trim2, "iso_results/spline_ppa_pnc6_trim2.rds")

# ───────────────────────────────────────
## Fit Isotonic Models ----
# ───────────────────────────────────────

iso_ppa_pnc6_25_trim2 <- flexstepreg_glmer_25(df_mat_ppa_pnc6_trim2$ppa_pnc6, df_mat_ppa_pnc6_trim2$hb, df_mat_ppa_pnc6_trim2$SITE, 
                                              covar2 = NULL, random_effect = NULL, alpha = 0.01)
saveRDS(iso_ppa_pnc6_25_trim2, "iso_results/iso_ppa_pnc6_25_trim2.rds")

iso_ppa_pnc6_50_trim2 <- flexstepreg_glmer(df_mat_ppa_pnc6_trim2$ppa_pnc6, df_mat_ppa_pnc6_trim2$hb, df_mat_ppa_pnc6_trim2$SITE, 
                                           covar2 = NULL, random_effect = NULL, alpha = 0.01)
saveRDS(iso_ppa_pnc6_50_trim2, "iso_results/iso_ppa_pnc6_50_trim2.rds")

# ───────────────────────────────────────
## Generate Plots ----
# ───────────────────────────────────────

# Save 0.50 plot
png("iso_results/vioplot_ppa_pnc6_50_trim2.png", width = 800, height = 600)
# 0.50 Isotonic
vioplot_ppa_pnc6_50_trim2 <- plot_boot_violin(
  x = df_mat_ppa_pnc6_trim2$hb,
  y = df_mat_ppa_pnc6_trim2$ppa_pnc6,
  xlab = "Hb",
  ylab = "PPA PNC6",
  rcs_result = spline_ppa_pnc6_trim2,
  iso_model = iso_ppa_pnc6_50_trim2,
  outcome_var = df_mat_ppa_pnc6_trim2$ppa_pnc6,
  title = "PPA PNC6 TRIM2 with Bootstrap Spline CI (0.50 Isotonic)"
)
dev.off()

# Save 0.25 plot
png("iso_results/vioplot_ppa_pnc6_25_trim2.png", width = 800, height = 600)
# 0.25 Isotonic
vioplot_ppa_pnc6_25_trim2 <- plot_boot_violin(
  x = df_mat_ppa_pnc6_trim2$hb,
  y = df_mat_ppa_pnc6_trim2$ppa_pnc6,
  xlab = "Hb",
  ylab = "PPA PNC6",
  rcs_result = spline_ppa_pnc6_trim2,
  iso_model = iso_ppa_pnc6_25_trim2,
  outcome_var = df_mat_ppa_pnc6_trim2$ppa_pnc6,
  title = "PPA PNC6 TRIM2 with Bootstrap Spline CI (0.25 Isotonic)"
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

out_ppa_pnc6_25_trim2  <- save_iso_output(iso_ppa_pnc6_25_trim2,  "PPA PNC6 (.25, Trim2)",  "out_ppa_pnc6_25_trim2.rda")
out_ppa_pnc6_50_trim2  <- save_iso_output(iso_ppa_pnc6_50_trim2,  "PPA PNC6 (.50, Trim2)",  "out_ppa_pnc6_50_trim2.rda")




# ───────────────────────────────────────
#TRIMESTER 3 – PPA PNC6 Models ----
# ───────────────────────────────────────

# ───────────────────────────────────────
## Fit Spline Models ----
# ───────────────────────────────────────
spline_ppa_pnc6_trim3 <- knot_fun_boot(df_mat_ppa_pnc6_trim3, "hb", "ppa_pnc6")
saveRDS(spline_ppa_pnc6_trim3, "iso_results/spline_ppa_pnc6_trim3.rds")

# ───────────────────────────────────────
## Fit Isotonic Models ----
# ───────────────────────────────────────

iso_ppa_pnc6_25_trim3 <- flexstepreg_glmer_25(df_mat_ppa_pnc6_trim3$ppa_pnc6, df_mat_ppa_pnc6_trim3$hb, df_mat_ppa_pnc6_trim3$SITE, 
                                              covar2 = NULL, random_effect = NULL, alpha = 0.01)
saveRDS(iso_ppa_pnc6_25_trim3, "iso_results/iso_ppa_pnc6_25_trim3.rds")

iso_ppa_pnc6_50_trim3 <- flexstepreg_glmer(df_mat_ppa_pnc6_trim3$ppa_pnc6, df_mat_ppa_pnc6_trim3$hb, df_mat_ppa_pnc6_trim3$SITE, 
                                           covar2 = NULL, random_effect = NULL, alpha = 0.01)
saveRDS(iso_ppa_pnc6_50_trim3, "iso_results/iso_ppa_pnc6_50_trim3.rds")

# ───────────────────────────────────────
## Generate Plots ----
# ───────────────────────────────────────

# Save 0.50 plot
png("iso_results/vioplot_ppa_pnc6_50_trim3.png", width = 800, height = 600)
# 0.50 Isotonic
vioplot_ppa_pnc6_50_trim3 <- plot_boot_violin(
  x = df_mat_ppa_pnc6_trim3$hb,
  y = df_mat_ppa_pnc6_trim3$ppa_pnc6,
  xlab = "Hb",
  ylab = "PPA PNC6",
  rcs_result = spline_ppa_pnc6_trim3,
  iso_model = iso_ppa_pnc6_50_trim3,
  outcome_var = df_mat_ppa_pnc6_trim3$ppa_pnc6,
  title = "PPA PNC6 TRIM3 with Bootstrap Spline CI (0.50 Isotonic)"
)
dev.off()

# Save 0.25 plot
png("iso_results/vioplot_ppa_pnc6_25_trim3.png", width = 800, height = 600)
# 0.25 Isotonic
vioplot_ppa_pnc6_25_trim3 <- plot_boot_violin(
  x = df_mat_ppa_pnc6_trim3$hb,
  y = df_mat_ppa_pnc6_trim3$ppa_pnc6,
  xlab = "Hb",
  ylab = "PPA PNC6",
  rcs_result = spline_ppa_pnc6_trim3,
  iso_model = iso_ppa_pnc6_25_trim3,
  outcome_var = df_mat_ppa_pnc6_trim3$ppa_pnc6,
  title = "PPA PNC6 TRIM3 with Bootstrap Spline CI (0.25 Isotonic)"
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

out_ppa_pnc6_25_trim3  <- save_iso_output(iso_ppa_pnc6_25_trim3,  "PPA PNC6 (.25, Trim3)",  "out_ppa_pnc6_25_trim3.rda")
out_ppa_pnc6_50_trim3  <- save_iso_output(iso_ppa_pnc6_50_trim3,  "PPA PNC6 (.50, Trim3)",  "out_ppa_pnc6_50_trim3.rda")
