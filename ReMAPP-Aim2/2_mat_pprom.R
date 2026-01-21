# ───────────────────────────────────────
# 📁 ReMAPP Aim 2 – PPROM Analysis Pipeline
# 🔧 Modified by: Williams Precious
# 📅 Last updated: 2025-04-18
# 📄 Description: Fits spline and isotonic models for PPROM, generates predictions and plots.
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

load("derived_data/df_mat_pprom.rda")
load("derived_data/df_mat_pprom_trim1.rda")
load("derived_data/df_mat_pprom_trim2.rda")
load("derived_data/df_mat_pprom_trim3.rda")


# for (df_name in c("df_mat_pprom_trim1", "df_mat_pprom_trim2", "df_mat_pprom_trim3")) {
#   df <- get(df_name)
#   
#   # If both hb.x and hb.y exist, keep hb.x
#   if ("hb.x" %in% names(df)) {
#     df$hb <- df$hb.x
#     df$hb.x <- NULL
#   }
#   if ("hb.y" %in% names(df)) {
#     df$hb.y <- NULL
#   }
#   
#   # Apply the transformation
#   if ("hb" %in% names(df)) {
#     df$hb <- round(pmin(pmax(df$hb, 5), 18), 1)
#   }
#   
#   assign(df_name, df)
# }
# 
# saveRDS(df_mat_pprom_trim1, "df_mat_pprom_trim1.rds")
# saveRDS(df_mat_pprom_trim2, "df_mat_pprom_trim2.rds")
# saveRDS(df_mat_pprom_trim3, "df_mat_pprom_trim3.rds")

df_mat_pprom$hb        <- round(pmin(pmax(df_mat_pprom$hb, 5), 18), 1)
df_mat_pprom_trim1$hb  <- round(pmin(pmax(df_mat_pprom_trim1$hb, 5), 18), 1)
df_mat_pprom_trim2$hb  <- round(pmin(pmax(df_mat_pprom_trim2$hb, 5), 18), 1)
df_mat_pprom_trim3$hb  <- round(pmin(pmax(df_mat_pprom_trim3$hb, 5), 18), 1)

# Establish the source of the functions
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/spline_binary_boot.R")
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/iso_binary_25_mod.R")
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/iso_binary_50_mod.R")
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/plot_boot.R")
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/outdata_mod.R")

# ───────────────────────────────────────
#ALL TRIMESTERS PPROM ----
# ───────────────────────────────────────

# ───────────────────────────────────────
## Fit Spline Models ----
# ───────────────────────────────────────
spline_pprom <- knot_fun_boot(df_mat_pprom, "hb", "pprom")
saveRDS(spline_pprom, "iso_results/spline_pprom.rds")

# ───────────────────────────────────────
## Fit Isotonic Models ----
# ───────────────────────────────────────

iso_pprom_25 <- flexstepreg_glmer_25(df_mat_pprom$pprom, df_mat_pprom$hb, df_mat_pprom$SITE, 
                                   covar2 = NULL, random_effect = NULL, alpha = 0.01)
saveRDS(iso_pprom_25, "iso_results/iso_pprom_25.rds")

iso_pprom_50 <- flexstepreg_glmer(df_mat_pprom$pprom, df_mat_pprom$hb, df_mat_pprom$SITE, 
                                covar2 = NULL, random_effect = NULL, alpha = 0.01)
saveRDS(iso_pprom_50, "iso_results/iso_pprom_50.rds")

# ───────────────────────────────────────
## Generate Plots ----
# ───────────────────────────────────────
# Save 0.50 plot
png("iso_results/vioplot_pprom_50.png", width = 800, height = 600)
# 0.50 Isotonic
vioplot_pprom_50 <- plot_boot_violin(
  x = df_mat_pprom$hb,
  y = df_mat_pprom$pprom,
  xlab = "Hb",
  ylab = "PPROM",
  rcs_result = spline_pprom,
  iso_model = iso_pprom_50,
  outcome_var = df_mat_pprom$pprom,
  title = "PPROM with Bootstrap Spline CI (0.50 Isotonic)"
)
dev.off()

# Save 0.25 plot
png("iso_results/vioplot_pprom_25.png", width = 800, height = 600)
# 0.25 Isotonic
vioplot_pprom_25 <- plot_boot_violin(
  x = df_mat_pprom$hb,
  y = df_mat_pprom$pprom,
  xlab = "Hb",
  ylab = "PPROM",
  rcs_result = spline_pprom,
  iso_model = iso_pprom_25,
  outcome_var = df_mat_pprom$pprom,
  title = "PPROM with Bootstrap Spline CI (0.25 Isotonic)"
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
out_pprom_25  <- save_iso_output(iso_pprom_25,  "PPROM (.25)",  "out_pprom_25.rda")
out_pprom_50  <- save_iso_output(iso_pprom_50,  "PPROM (.50)",  "out_pprom_50.rda")


# ───────────────────────────────────────
#TRIMESTER 1 – PPROM Models ----
# ───────────────────────────────────────

# ───────────────────────────────────────
## Fit Spline Models ----
# ───────────────────────────────────────
spline_pprom_trim1 <- knot_fun_boot(df_mat_pprom_trim1, "hb", "pprom")
saveRDS(spline_pprom_trim1, "iso_results/spline_pprom_trim1.rds")

# ───────────────────────────────────────
## Fit Isotonic Models ----
# ───────────────────────────────────────

iso_pprom_25_trim1 <- flexstepreg_glmer_25(df_mat_pprom_trim1$pprom, df_mat_pprom_trim1$hb, df_mat_pprom_trim1$SITE, 
                                         covar2 = NULL, random_effect = NULL, alpha = 0.01)
saveRDS(iso_pprom_25_trim1, "iso_results/iso_pprom_25_trim1.rds")

iso_pprom_50_trim1 <- flexstepreg_glmer(df_mat_pprom_trim1$pprom, df_mat_pprom_trim1$hb, df_mat_pprom_trim1$SITE, 
                                      covar2 = NULL, random_effect = NULL, alpha = 0.01)
saveRDS(iso_pprom_50_trim1, "iso_results/iso_pprom_50_trim1.rds")

# ───────────────────────────────────────
## Generate Plots ----
# ───────────────────────────────────────

# Save 0.50 plot
png("iso_results/vioplot_pprom_50_trim1.png", width = 800, height = 600)
# 0.50 Isotonic
vioplot_pprom_50_trim1 <- plot_boot_violin(
  x = df_mat_pprom_trim1$hb,
  y = df_mat_pprom_trim1$pprom,
  xlab = "Hb",
  ylab = "PPROM",
  rcs_result = spline_pprom_trim1,
  iso_model = iso_pprom_50_trim1,
  outcome_var = df_mat_pprom_trim1$pprom,
  title = "PPROM TRIM1 with Bootstrap Spline CI (0.50 Isotonic)"
)
dev.off()

# Save 0.25 plot
png("iso_results/vioplot_pprom_25_trim1.png", width = 800, height = 600)
# 0.25 Isotonic
vioplot_pprom_25_trim1 <- plot_boot_violin(
  x = df_mat_pprom_trim1$hb,
  y = df_mat_pprom_trim1$pprom,
  xlab = "Hb",
  ylab = "PPROM",
  rcs_result = spline_pprom_trim1,
  iso_model = iso_pprom_25_trim1,
  outcome_var = df_mat_pprom_trim1$pprom,
  title = "PPROM TRIM1 with Bootstrap Spline CI (0.25 Isotonic)"
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

out_pprom_25_trim1  <- save_iso_output(iso_pprom_25_trim1,  "PPROM (.25, Trim1)",  "out_pprom_25_trim1.rda")
out_pprom_50_trim1  <- save_iso_output(iso_pprom_50_trim1,  "PPROM (.50, Trim1)",  "out_pprom_50_trim1.rda")




# ───────────────────────────────────────
#TRIMESTER 2 – PPROM Models ----
# ───────────────────────────────────────

# ───────────────────────────────────────
## Fit Spline Models ----
# ───────────────────────────────────────
spline_pprom_trim2 <- knot_fun_boot(df_mat_pprom_trim2, "hb", "pprom")
saveRDS(spline_pprom_trim2, "iso_results/spline_pprom_trim2.rds")

# ───────────────────────────────────────
## Fit Isotonic Models ----
# ───────────────────────────────────────

iso_pprom_25_trim2 <- flexstepreg_glmer_25(df_mat_pprom_trim2$pprom, df_mat_pprom_trim2$hb, df_mat_pprom_trim2$SITE, 
                                         covar2 = NULL, random_effect = NULL, alpha = 0.01)
saveRDS(iso_pprom_25_trim2, "iso_results/iso_pprom_25_trim2.rds")

iso_pprom_50_trim2 <- flexstepreg_glmer(df_mat_pprom_trim2$pprom, df_mat_pprom_trim2$hb, df_mat_pprom_trim2$SITE, 
                                      covar2 = NULL, random_effect = NULL, alpha = 0.01)
saveRDS(iso_pprom_50_trim2, "iso_results/iso_pprom_50_trim2.rds")

# ───────────────────────────────────────
## Generate Plots ----
# ───────────────────────────────────────
# Save 0.50 plot
png("iso_results/vioplot_pprom_50_trim2.png", width = 800, height = 600)
# 0.50 Isotonic
vioplot_pprom_50_trim2 <- plot_boot_violin(
  x = df_mat_pprom_trim2$hb,
  y = df_mat_pprom_trim2$pprom,
  xlab = "Hb",
  ylab = "PPROM",
  rcs_result = spline_pprom_trim2,
  iso_model = iso_pprom_50_trim2,
  outcome_var = df_mat_pprom_trim2$pprom,
  title = "PPROM TRIM2 with Bootstrap Spline CI (0.50 Isotonic)"
)
dev.off()

# Save 0.25 plot
png("iso_results/vioplot_pprom_25_trim2.png", width = 800, height = 600)
# 0.25 Isotonic
vioplot_pprom_25_trim2 <- plot_boot_violin(
  x = df_mat_pprom_trim2$hb,
  y = df_mat_pprom_trim2$pprom,
  xlab = "Hb",
  ylab = "PPROM",
  rcs_result = spline_pprom_trim2,
  iso_model = iso_pprom_25_trim2,
  outcome_var = df_mat_pprom_trim2$pprom,
  title = "PPROM TRIM2 with Bootstrap Spline CI (0.25 Isotonic)"
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

out_pprom_25_trim2  <- save_iso_output(iso_pprom_25_trim2,  "PPROM (.25, Trim2)",  "out_pprom_25_trim2.rda")
out_pprom_50_trim2  <- save_iso_output(iso_pprom_50_trim2,  "PPROM (.50, Trim2)",  "out_pprom_50_trim2.rda")




# ───────────────────────────────────────
#TRIMESTER 3 – PPROM Models ----
# ───────────────────────────────────────

# ───────────────────────────────────────
## Fit Spline Models ----
# ───────────────────────────────────────
spline_pprom_trim3 <- knot_fun_boot(df_mat_pprom_trim3, "hb", "pprom")
saveRDS(spline_pprom_trim3, "iso_results/spline_pprom_trim3.rds")

# ───────────────────────────────────────
## Fit Isotonic Models ----
# ───────────────────────────────────────

iso_pprom_25_trim3 <- flexstepreg_glmer_25(df_mat_pprom_trim3$pprom, df_mat_pprom_trim3$hb, df_mat_pprom_trim3$SITE, 
                                         covar2 = NULL, random_effect = NULL, alpha = 0.01)
saveRDS(iso_pprom_25_trim3, "iso_results/iso_pprom_25_trim3.rds")

iso_pprom_50_trim3 <- flexstepreg_glmer(df_mat_pprom_trim3$pprom, df_mat_pprom_trim3$hb, df_mat_pprom_trim3$SITE, 
                                      covar2 = NULL, random_effect = NULL, alpha = 0.01)
saveRDS(iso_pprom_50_trim3, "iso_results/iso_pprom_50_trim3.rds")

# ───────────────────────────────────────
## Generate Plots ----
# ───────────────────────────────────────

# Save 0.50 plot
png("iso_results/vioplot_pprom_50_trim3.png", width = 800, height = 600)
# 0.50 Isotonic
vioplot_pprom_50_trim3 <- plot_boot_violin(
  x = df_mat_pprom_trim3$hb,
  y = df_mat_pprom_trim3$pprom,
  xlab = "Hb",
  ylab = "PPROM",
  rcs_result = spline_pprom,
  iso_model = iso_pprom_50_trim3,
  outcome_var = df_mat_pprom_trim3$pprom,
  title = "PPROM TRIM3 with Bootstrap Spline CI (0.50 Isotonic)"
)
dev.off()

# Save 0.25 plot
png("iso_results/vioplot_pprom_25_trim3.png", width = 800, height = 600)
# 0.25 Isotonic
vioplot_pprom_25_trim3 <- plot_boot_violin(
  x = df_mat_pprom_trim3$hb,
  y = df_mat_pprom_trim3$pprom,
  xlab = "Hb",
  ylab = "PPROM",
  rcs_result = spline_pprom,
  iso_model = iso_pprom_25_trim3,
  outcome_var = df_mat_pprom_trim3$pprom,
  title = "PPROM TRIM3 with Bootstrap Spline CI (0.25 Isotonic)"
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

out_pprom_25_trim3  <- save_iso_output(iso_pprom_25_trim3,  "PPROM (.25, Trim3)",  "out_pprom_25_trim3.rda")
out_pprom_50_trim3  <- save_iso_output(iso_pprom_50_trim3,  "PPROM (.50, Trim3)",  "out_pprom_50_trim3.rda")
