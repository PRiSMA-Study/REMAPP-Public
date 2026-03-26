# ───────────────────────────────────────
# 📁 ReMAPP Aim 2 – Fatigue score Analysis Pipeline
# 🔧 Modified by: Williams Precious
# 📅 Last updated: 2025-04-18
# 📄 Description: Fits spline and isotonic models for Fatigue at ANC and PNC, generates predictions and plots.
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

#prepare data
# Load the dataset
# Load data
load("derived_data/df_mat_anc_fat.rda")
load("derived_data/df_mat_anc_fat_trim1.rda")
load("derived_data/df_mat_anc_fat_trim2.rda")
load("derived_data/df_mat_anc_fat_trim3.rda")
load("derived_data/df_mat_pnc_fat_pnc6.rda")

df_mat_anc_fat$hb        <- round(pmin(pmax(df_mat_anc_fat$hb, 5), 18), 1)
df_mat_anc_fat_trim1$hb  <- round(pmin(pmax(df_mat_anc_fat_trim1$hb, 5), 18), 1)
df_mat_anc_fat_trim2$hb  <- round(pmin(pmax(df_mat_anc_fat_trim2$hb, 5), 18), 1)
df_mat_anc_fat_trim3$hb  <- round(pmin(pmax(df_mat_anc_fat_trim3$hb, 5), 18), 1)
df_mat_pnc_fat_pnc6$hb  <- round(pmin(pmax(df_mat_pnc_fat_pnc6$hb, 5), 18), 1)

# Establish the source of the functions
# source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/spline_binary_boot.R")
# source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/iso_binary_25_mod.R")
# source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/iso_binary_50_mod.R")
# source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/plot_boot.R")
# source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/outdata_mod.R")


source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/spline_continuous_boot.R")
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/iso_continuous_50.R")
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/iso_continuous_25.R")
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/plot_boot_continuous.R")
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/outdata_mod.R")
# ───────────────────────────────────────
#ALL TRIMESTERS FATIGUE ANC ----
# ───────────────────────────────────────

# ───────────────────────────────────────
## Fit Spline Models ----
# ───────────────────────────────────────
spline_anc_fat <- knot_fun_boot(df_mat_anc_fat, "hb", "fatigue_score")
saveRDS(spline_anc_fat, "iso_results/spline_anc_fat.rds")

# ───────────────────────────────────────
## Fit Isotonic Models ----
# ───────────────────────────────────────

iso_anc_fat_25 <- flexstepreg_lmer25(df_mat_anc_fat$fatigue_score, df_mat_anc_fat$hb, df_mat_anc_fat$SITE,
                                       covar2 = NULL, random_effect = df_mat_anc_fat$MOMID, alpha = 0.01)
saveRDS(iso_anc_fat_25, "iso_results/iso_anc_fat_25.rds")

iso_anc_fat_50 <- flexstepreg_lmer(df_mat_anc_fat$fatigue_score, df_mat_anc_fat$hb, df_mat_anc_fat$SITE, 
                                    covar2 = NULL, random_effect = df_mat_anc_fat$MOMID, alpha = 0.01)
saveRDS(iso_anc_fat_50, "iso_results/iso_anc_fat_50.rds")

# ───────────────────────────────────────
## Generate Plots ----
# ───────────────────────────────────────

# Save 0.50 plot
png("iso_results/vioplot_anc_fat_50.png", width = 800, height = 600)
# 0.50 Isotonic
vioplot_anc_fat_50 <- plot_boot_continuous(
  x = df_mat_anc_fat$hb,
  y = df_mat_anc_fat$fatigue_score,
  xlab = "Hb",
  ylab = "FATIGUE ANC",
  rcs_result = spline_anc_fat,
  iso_model = iso_anc_fat_50,
  outcome_var = df_mat_anc_fat$fatigue_score,
  title = "FATIGUE ANC with Bootstrap Spline CI (0.50 Isotonic)"
)

dev.off()

# Save 0.25 plot
png("iso_results/vioplot_anc_fat_25.png", width = 800, height = 600)
# 0.25 Isotonic
vioplot_anc_fat_25 <- plot_boot_continuous(
  x = df_mat_anc_fat$hb,
  y = df_mat_anc_fat$fatigue_score,
  xlab = "Hb",
  ylab = "FATIGUE ANC",
  rcs_result = spline_anc_fat,
  iso_model = iso_anc_fat_25,
  outcome_var = df_mat_anc_fat$fatigue_score,
  title = "FATIGUE ANC with Bootstrap Spline CI (0.25 Isotonic)"
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
out_anc_fat_25  <- save_iso_output(iso_anc_fat_25,  "FATIGUE ANC (.25)",  "out_anc_fat_25.rda")
out_anc_fat_50  <- save_iso_output(iso_anc_fat_50,  "FATIGUE ANC (.50)",  "out_anc_fat_50.rda")


#PJW Comment - There is only one 1st trimester Fatigue measurement
# # ───────────────────────────────────────
# #TRIMESTER 1 – FATIGUE ANC Models ----
# # ───────────────────────────────────────
# 
# # ───────────────────────────────────────
# ## Fit Spline Models ----
# # ───────────────────────────────────────
# spline_anc_fat_trim1 <- knot_fun_boot(df_mat_anc_fat_trim1, "hb", "fatigue_score")
# saveRDS(spline_anc_fat_trim1, "iso_results/spline_anc_fat_trim1.rds")
# 
# # ───────────────────────────────────────
# ## Fit Isotonic Models ----
# # ───────────────────────────────────────
# 
# iso_anc_fat_25_trim1 <- flexstepreg_lmer25(df_mat_anc_fat_trim1$fatigue_score, df_mat_anc_fat_trim1$hb, df_mat_anc_fat_trim1$SITE, 
#                                              covar2 = NULL, random_effect = df_mat_anc_fat_trim1$MOMID, alpha = 0.01)
# saveRDS(iso_anc_fat_25_trim1, "iso_results/iso_anc_fat_25_trim1.rds")
# 
# iso_anc_fat_50_trim1 <- flexstepreg_lmer(df_mat_anc_fat_trim1$fatigue_score, df_mat_anc_fat_trim1$hb, df_mat_anc_fat_trim1$SITE, 
#                                           covar2 = NULL, random_effect = df_mat_anc_fat_trim1$MOMID, alpha = 0.01)
# saveRDS(iso_anc_fat_50_trim1, "iso_results/iso_anc_fat_50_trim1.rds")
# 
# # ───────────────────────────────────────
# ## Generate Plots ----
# # ───────────────────────────────────────
# 
# # Save 0.50 plot
# png("iso_results/vioplot_anc_fat_50_trim1.png", width = 800, height = 600)
# # 0.50 Isotonic
# vioplot_anc_fat_50_trim1 <- plot_boot_continuous(
#   x = df_mat_anc_fat_trim1$hb,
#   y = df_mat_anc_fat_trim1$fatigue_score,
#   xlab = "Hb",
#   ylab = "FATIGUE ANC",
#   rcs_result = spline_anc_fat_trim1,
#   iso_model = iso_anc_fat_50_trim1,
#   outcome_var = df_mat_anc_fat_trim1$fatigue_score,
#   title = "FATIGUE ANC TRIM1 with Bootstrap Spline CI (0.50 Isotonic)"
# )
# dev.off()
# 
# # Save 0.25 plot
# png("iso_results/vioplot_anc_fat_25_trim1.png", width = 800, height = 600)
# # 0.25 Isotonic
# vioplot_anc_fat_25_trim1 <- plot_boot_continuous(
#   x = df_mat_anc_fat_trim1$hb,
#   y = df_mat_anc_fat_trim1$fatigue_score,
#   xlab = "Hb",
#   ylab = "FATIGUE ANC",
#   rcs_result = spline_anc_fat_trim1,
#   iso_model = iso_anc_fat_25_trim1,
#   outcome_var = df_mat_anc_fat_trim1$fatigue_score,
#   title = "FATIGUE ANC TRIM1 with Bootstrap Spline CI (0.25 Isotonic)"
# )
# 
# dev.off()
# 
# 
# # ───────────────────────────────────────
# ## Save Processed Isotonic Output Objects ----
# # ───────────────────────────────────────
# save_iso_output <- function(model_obj, label, filename) {
#   out <- outdata(model_obj, label)
#   save(out, file = file.path("iso_results", filename))
#   return(out)
# }
# 
# out_anc_fat_25_trim1  <- save_iso_output(iso_anc_fat_25_trim1,  "FATIGUE ANC(.25, Trim1)",  "out_anc_fat_25_trim1.rda")
# out_anc_fat_50_trim1  <- save_iso_output(iso_anc_fat_50_trim1,  "FATIGUE ANC(.50, Trim1)",  "out_anc_fat_50_trim1.rda")
# 
# 


# ───────────────────────────────────────
#TRIMESTER 2 – FATIGUE ANC Models ----
# ───────────────────────────────────────

# ───────────────────────────────────────
## Fit Spline Models ----
# ───────────────────────────────────────
spline_anc_fat_trim2 <- knot_fun_boot(df_mat_anc_fat_trim2, "hb", "fatigue_score")
saveRDS(spline_anc_fat_trim2, "iso_results/spline_anc_fat_trim2.rds")

# ───────────────────────────────────────
## Fit Isotonic Models ----
# ───────────────────────────────────────

iso_anc_fat_25_trim2 <- flexstepreg_lmer25(df_mat_anc_fat_trim2$fatigue_score, df_mat_anc_fat_trim2$hb, df_mat_anc_fat_trim2$SITE, 
                                             covar2 = NULL, random_effect = df_mat_anc_fat_trim2$MOMID, alpha = 0.01)
saveRDS(iso_anc_fat_25_trim2, "iso_results/iso_anc_fat_25_trim2.rds")

iso_anc_fat_50_trim2 <- flexstepreg_lmer(df_mat_anc_fat_trim2$fatigue_score, df_mat_anc_fat_trim2$hb, df_mat_anc_fat_trim2$SITE, 
                                          covar2 = NULL, random_effect = df_mat_anc_fat_trim2$MOMID, alpha = 0.01)
saveRDS(iso_anc_fat_50_trim2, "iso_results/iso_anc_fat_50_trim2.rds")

# ───────────────────────────────────────
## Generate Plots ----
# ───────────────────────────────────────

# Save 0.50 plot
png("iso_results/vioplot_anc_fat_50_trim2.png", width = 800, height = 600)
# 0.50 Isotonic
vioplot_anc_fat_50_trim2 <- plot_boot_continuous(
  x = df_mat_anc_fat_trim2$hb,
  y = df_mat_anc_fat_trim2$fatigue_score,
  xlab = "Hb",
  ylab = "FATIGUE ANC",
  rcs_result = spline_anc_fat_trim2,
  iso_model = iso_anc_fat_50_trim2,
  outcome_var = df_mat_anc_fat_trim2$fatigue_score,
  title = "FATIGUE ANC TRIM2 with Bootstrap Spline CI (0.50 Isotonic)"
)
dev.off()

# Save 0.25 plot
png("iso_results/vioplot_anc_fat_25_trim2.png", width = 800, height = 600)
# 0.25 Isotonic
vioplot_anc_fat_25_trim2 <- plot_boot_continuous(
  x = df_mat_anc_fat_trim2$hb,
  y = df_mat_anc_fat_trim2$fatigue_score,
  xlab = "Hb",
  ylab = "FATIGUE ANC",
  rcs_result = spline_anc_fat_trim2,
  iso_model = iso_anc_fat_25_trim2,
  outcome_var = df_mat_anc_fat_trim2$fatigue_score,
  title = "FATIGUE ANC TRIM2 with Bootstrap Spline CI (0.25 Isotonic)"
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

out_anc_fat_25_trim2  <- save_iso_output(iso_anc_fat_25_trim2,  "FATIGUE ANC(.25, Trim2)",  "out_anc_fat_25_trim2.rda")
out_anc_fat_50_trim2  <- save_iso_output(iso_anc_fat_50_trim2,  "FATIGUE ANC(.50, Trim2)",  "out_anc_fat_50_trim2.rda")




# ───────────────────────────────────────
#TRIMESTER 3 – FATIGUE ANC Models ----
# ───────────────────────────────────────

# ───────────────────────────────────────
## Fit Spline Models ----
# ───────────────────────────────────────
spline_anc_fat_trim3 <- knot_fun_boot(df_mat_anc_fat_trim3, "hb", "fatigue_score")
saveRDS(spline_anc_fat_trim3, "iso_results/spline_anc_fat_trim3.rds")

# ───────────────────────────────────────
## Fit Isotonic Models ----
# ───────────────────────────────────────

iso_anc_fat_25_trim3 <- flexstepreg_lmer25(df_mat_anc_fat_trim3$fatigue_score, df_mat_anc_fat_trim3$hb, df_mat_anc_fat_trim3$SITE, 
                                             covar2 = NULL, random_effect = df_mat_anc_fat_trim3$MOMID, alpha = 0.01)
saveRDS(iso_anc_fat_25_trim3, "iso_results/iso_anc_fat_25_trim3.rds")

iso_anc_fat_50_trim3 <- flexstepreg_lmer(df_mat_anc_fat_trim3$fatigue_score, df_mat_anc_fat_trim3$hb, df_mat_anc_fat_trim3$SITE, 
                                          covar2 = NULL, random_effect = df_mat_anc_fat_trim3$MOMID, alpha = 0.01)
saveRDS(iso_anc_fat_50_trim3, "iso_results/iso_anc_fat_50_trim3.rds")

# ───────────────────────────────────────
## Generate Plots ----
# ───────────────────────────────────────

# Save 0.50 plot
png("iso_results/vioplot_anc_fat_50_trim3.png", width = 800, height = 600)
# 0.50 Isotonic
vioplot_anc_fat_50_trim3 <- plot_boot_continuous(
  x = df_mat_anc_fat_trim3$hb,
  y = df_mat_anc_fat_trim3$fatigue_score,
  xlab = "Hb",
  ylab = "FATIGUE",
  rcs_result = spline_anc_fat_trim3,
  iso_model = iso_anc_fat_50_trim3,
  outcome_var = df_mat_anc_fat_trim3$fatigue_score,
  title = "FATIGUE ANC TRIM3 with Bootstrap Spline CI (0.50 Isotonic)"
)
dev.off()

# Save 0.25 plot
png("iso_results/vioplot_anc_fat_25_trim3.png", width = 800, height = 600)
# 0.25 Isotonic
vioplot_anc_fat_25_trim3 <- plot_boot_continuous(
  x = df_mat_anc_fat_trim3$hb,
  y = df_mat_anc_fat_trim3$fatigue_score,
  xlab = "Hb",
  ylab = "FATIGUE",
  rcs_result = spline_anc_fat_trim3,
  iso_model = iso_anc_fat_25_trim3,
  outcome_var = df_mat_anc_fat_trim3$fatigue_score,
  title = "FATIGUE ANC TRIM3 with Bootstrap Spline CI (0.25 Isotonic)"
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

out_anc_fat_25_trim3  <- save_iso_output(iso_anc_fat_25_trim3,  "FATIGUE ANC(.25, Trim3)",  "out_anc_fat_25_trim3.rda")
out_anc_fat_50_trim3  <- save_iso_output(iso_anc_fat_50_trim3,  "FATIGUE ANC(.50, Trim3)",  "out_anc_fat_50_trim3.rda")


# ───────────────────────────────────────
#PNC6 – FATIGUE PNC Models ----
# ───────────────────────────────────────

# ───────────────────────────────────────
## Fit Spline Models ----
# ───────────────────────────────────────
spline_pnc_fat_pnc6 <- knot_fun_boot(df_mat_pnc_fat_pnc6, "hb", "fatigue_score")
saveRDS(spline_pnc_fat_pnc6, "iso_results/spline_pnc_fat_pnc6.rds")

# ───────────────────────────────────────
## Fit Isotonic Models ----
# ───────────────────────────────────────

iso_pnc_fat_25_pnc6 <- flexstepreg_lmer25(df_mat_pnc_fat_pnc6$fatigue_score, df_mat_pnc_fat_pnc6$hb, df_mat_pnc_fat_pnc6$SITE, 
                                           covar2 = NULL, random_effect = df_mat_pnc_fat_pnc6$MOMID, alpha = 0.01)
saveRDS(iso_pnc_fat_25_pnc6, "iso_results/iso_pnc_fat_25_pnc6.rds")

iso_pnc_fat_50_pnc6 <- flexstepreg_lmer(df_mat_pnc_fat_pnc6$fatigue_score, df_mat_pnc_fat_pnc6$hb, df_mat_pnc_fat_pnc6$SITE, 
                                         covar2 = NULL, random_effect = df_mat_pnc_fat_pnc6$MOMID, alpha = 0.01)
saveRDS(iso_pnc_fat_50_pnc6, "iso_results/iso_pnc_fat_50_pnc6.rds")

# ───────────────────────────────────────
## Generate Plots ----
# ───────────────────────────────────────

# Save 0.50 plot
png("iso_results/vioplot_pnc_fat_50_pnc6.png", width = 800, height = 600)
# 0.50 Isotonic
vioplot_pnc_fat_50_pnc6 <- plot_boot_continuous(
  x = df_mat_pnc_fat_pnc6$hb,
  y = df_mat_pnc_fat_pnc6$fatigue_score,
  xlab = "Hb",
  ylab = "FATIGUE",
  rcs_result = spline_pnc_fat_pnc6,
  iso_model = iso_pnc_fat_50_pnc6,
  outcome_var = df_mat_pnc_fat_pnc6$fatigue_score,
  title = "PP 6wks Fatigue vs PNC6 Hb with Bootstrap Spline CI (0.50 Isotonic)"
)
dev.off()

# Save 0.25 plot
png("iso_results/vioplot_pnc_fat_25_pnc6.png", width = 800, height = 600)
# 0.25 Isotonic
vioplot_pnc_fat_25_pnc6 <- plot_boot_continuous(
  x = df_mat_pnc_fat_pnc6$hb,
  y = df_mat_pnc_fat_pnc6$fatigue_score,
  xlab = "Hb",
  ylab = "FATIGUE",
  rcs_result = spline_pnc_fat_pnc6,
  iso_model = iso_pnc_fat_25_pnc6,
  outcome_var = df_mat_pnc_fat_pnc6$fatigue_score,
  title = "PP 6wks Fatigue vs PNC6 Hb with Bootstrap Spline CI (0.25 Isotonic)"
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

out_pnc_fat_25_pnc6  <- save_iso_output(iso_pnc_fat_25_pnc6,  "FATIGUE AT 6WKS(.25, PNC6)",  "out_pnc_fat_25_pnc6.rda")
out_pnc_fat_50_pnc6  <- save_iso_output(iso_pnc_fat_50_pnc6,  "FATIGUE AT 6WKS(.50, PNC6)",  "out_pnc_fat_50_pnc6.rda")


