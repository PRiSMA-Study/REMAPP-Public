# ───────────────────────────────────────
# 📁 ReMAPP Aim 2 – SGA Analysis Pipeline
# 🔧 Modified by: Williams Precious
# 📅 Last updated: 2025-04-18
# 📄 Description: Fits spline and isotonic models for SGA, generates predictions and plots.
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
# 📁 Setup Paths & Load Data ----
# ───────────────────────────────────────
UploadDate <- "2025-04-18"
base_dir <- file.path("D:/Users/williams_pj/Documents/Analysis/ReMAPP/Aim2", UploadDate)
setwd(base_dir)

load("derived_data/df_inf_sga10.rda")
load("derived_data/df_inf_sga10_trim1.rda")
load("derived_data/df_inf_sga10_trim2.rda")
load("derived_data/df_inf_sga10_trim3.rda")
# ───────────────────────────────────────
# 📈 ALL TRIMESTERS SGA ----
# ───────────────────────────────────────

# ───────────────────────────────────────
# 📈 Fit Spline Models ----
# ───────────────────────────────────────
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/spline_binary_5_18.R")
spline_sga10_5_18 <- knot_fun(data = df_inf_sga10,
                             hb_var = "hb",
                             outcome_var = "sga10")
saveRDS(spline_sga10_5_18, "iso_results/spline_sga10.rds")

source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/spline_binary_boot.R")
spline_sga10_boot <- knot_fun(df_inf_sga10, "hb", "sga10")
spline_sga10_bmodel <- spline_sga10_boot$model
spline_sga10_boot_ci <- spline_sga10_boot$boot_ci

saveRDS(spline_sga10_bmodel, "iso_results/spline_sga10_bmodel.rds")
saveRDS(spline_sga10_boot_ci, "iso_results/spline_sga10_boot_ci.rds")
save(spline_sga10_boot, file = "iso_results/spline_sga10_boot.rda")
# ───────────────────────────────────────
# 📉 Fit Isotonic Models ----
# ───────────────────────────────────────
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/iso_binary.R")
iso_sga10 <- flexstepreg_glmer(df_inf_sga10$sga10, df_inf_sga10$hb, df_inf_sga10$SITE, 
                              covar2 = NULL, random_effect = NULL, alpha = 0.01)
saveRDS(iso_sga10, "iso_results/iso_sga10.rds")

source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/iso_binary_25.R")
iso_sga10_25 <- flexstepreg_glmer(df_inf_sga10$sga10, df_inf_sga10$hb, df_inf_sga10$SITE, 
                                 covar2 = NULL, random_effect = NULL, alpha = 0.01)
saveRDS(iso_sga10_25, "iso_results/iso_sga10_25.rds")

source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/iso_binary_50.R")
iso_sga10_50 <- flexstepreg_glmer(df_inf_sga10$sga10, df_inf_sga10$hb, df_inf_sga10$SITE, 
                                 covar2 = NULL, random_effect = NULL, alpha = 0.01)
saveRDS(iso_sga10_50, "iso_results/iso_sga10_50.rds")

# ───────────────────────────────────────
# 📊 Generate Plots ----
# ───────────────────────────────────────
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/plot_binary.R")
isoplot_sga10 <- iso_fun_glmer(
  x = df_inf_sga10$hb,
  y = df_inf_sga10$sga10,
  ylab = "SGA <10rd %ile",
  spline_model = spline_sga10_5_18,
  iso_model = iso_sga10_50
)

source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/plot_violin.R")
vioplot_sga10 <- spline_and_violin_plot(
  x = df_inf_sga10$hb,
  y = df_inf_sga10$sga10,
  xlab = "Hb",
  ylab = "SGA 10rd %ile",
  spline_model = spline_sga10_5_18,
  iso_model = iso_sga10,
  outcome_var = df_inf_sga10$sga10,
  title = "SGA 10rd %ile"
)

source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/plot_boot.R")

plot_boot_violin(
  x = df_inf_sga10$hb,
  y = df_inf_sga10$sga10,
  xlab = "Hb",
  ylab = "SGA 10rd %ile",
  rcs_result = spline_sga10_boot,
  iso_model = iso_sga10_50,
  outcome_var = df_inf_sga10$sga10,
  title = "SGA 10rd %ile with Bootstrap Spline CI"
)


# ───────────────────────────────────────
# 💾 Save Plots ----
# ───────────────────────────────────────
png("iso_results/isoplot_sga10.png")
print(isoplot_sga10)
dev.off()

png("iso_results/vioplot_sga10.png")
print(vioplot_sga10)
dev.off()


# ───────────────────────────────────────
# 💾 Save Processed Isotonic Output Objects ----
# ───────────────────────────────────────
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/outdata.R")

# Helper to format and save isotonic model outputs
save_iso_output <- function(model_obj, label, filename) {
  out <- outdata(model_obj, label)
  save(out, file = file.path("iso_results", filename))
  return(out)
}

# Save outputs
out_sga10     <- save_iso_output(iso_sga10,     "SGA 10rd %ile",        "out_sga10.rda")
out_sga10_25  <- save_iso_output(iso_sga10_25,  "SGA 10rd %ile (.25)",  "out_sga10_25.rda")
out_sga10_50  <- save_iso_output(iso_sga10_50,  "SGA 10rd %ile (.50)",  "out_sga10_50.rda")


# ───────────────────────────────────────
# 📈 TRIMESTER 1 – SGA Models
# ───────────────────────────────────────

# ───────────────────────────────────────
# 📈 Fit Spline Models ----
# ───────────────────────────────────────
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/spline_binary_5_18.R")
spline_sga10_5_18_trim1 <- knot_fun(data = df_inf_sga10_trim1,
                                   hb_var = "hb",
                                   outcome_var = "sga10")
saveRDS(spline_sga10_5_18_trim1, "iso_results/spline_sga10_trim1.rds")

source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/spline_binary_boot.R")
spline_sga10_boot_trim1 <- knot_fun(df_inf_sga10_trim1, "hb", "sga10")
spline_sga10_bmodel_trim1 <- spline_sga10_boot_trim1$model
spline_sga10_boot_ci_trim1 <- spline_sga10_boot_trim1$boot_ci

saveRDS(spline_sga10_bmodel_trim1, "iso_results/spline_sga10_bmodel_trim1.rds")
saveRDS(spline_sga10_boot_ci_trim1, "iso_results/spline_sga10_boot_ci_trim1.rds")
save(spline_sga10_boot_trim1, file = "iso_results/spline_sga10_boot_trim1.rda")

# ───────────────────────────────────────
# 📉 Fit Isotonic Models ----
# ───────────────────────────────────────
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/iso_binary.R")
iso_sga10_trim1 <- flexstepreg_glmer(df_inf_sga10_trim1$sga10, df_inf_sga10_trim1$hb, df_inf_sga10_trim1$SITE, 
                                    covar2 = NULL, random_effect = NULL, alpha = 0.01)
saveRDS(iso_sga10_trim1, "iso_results/iso_sga10_trim1.rds")

source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/iso_binary_25.R")
iso_sga10_25_trim1 <- flexstepreg_glmer(df_inf_sga10_trim1$sga10, df_inf_sga10_trim1$hb, df_inf_sga10_trim1$SITE, 
                                       covar2 = NULL, random_effect = NULL, alpha = 0.01)
saveRDS(iso_sga10_25_trim1, "iso_results/iso_sga10_25_trim1.rds")

source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/iso_binary_50.R")
iso_sga10_50_trim1 <- flexstepreg_glmer(df_inf_sga10_trim1$sga10, df_inf_sga10_trim1$hb, df_inf_sga10_trim1$SITE, 
                                       covar2 = NULL, random_effect = NULL, alpha = 0.01)
saveRDS(iso_sga10_50_trim1, "iso_results/iso_sga10_50_trim1.rds")

# ───────────────────────────────────────
# 📊 Generate Plots ----
# ───────────────────────────────────────
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/plot_binary.R")
isoplot_sga10_trim1 <- iso_fun_glmer(
  x = df_inf_sga10_trim1$hb,
  y = df_inf_sga10_trim1$sga10,
  ylab = "SGA <10rd %ile",
  spline_model = spline_sga10_5_18_trim1,
  iso_model = iso_sga10_50_trim1
)

source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/plot_violin.R")
vioplot_sga10_trim1 <- spline_and_violin_plot(
  x = df_inf_sga10_trim1$hb,
  y = df_inf_sga10_trim1$sga10,
  xlab = "Hb",
  ylab = "SGA 10rd %ile",
  spline_model = spline_sga10_5_18_trim1,
  iso_model = iso_sga10_trim1,
  outcome_var = df_inf_sga10_trim1$sga10,
  title = "SGA 10rd %ile (Trimester 1)"
)

# ───────────────────────────────────────
# 💾 Save Plots ----
# ───────────────────────────────────────
png("iso_results/isoplot_sga10_trim1.png")
print(isoplot_sga10_trim1)
dev.off()

png("iso_results/vioplot_sga10_trim1.png")
print(vioplot_sga10_trim1)
dev.off()

# ───────────────────────────────────────
# 💾 Save Processed Isotonic Output Objects ----
# ───────────────────────────────────────
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/outdata.R")

save_iso_output <- function(model_obj, label, filename) {
  out <- outdata(model_obj, label)
  save(out, file = file.path("iso_results", filename))
  return(out)
}

out_sga10_trim1     <- save_iso_output(iso_sga10_trim1,     "SGA 10rd %ile (Trim1)",       "out_sga10_trim1.rda")
out_sga10_25_trim1  <- save_iso_output(iso_sga10_25_trim1,  "SGA 10rd %ile (.25, Trim1)",  "out_sga10_25_trim1.rda")
out_sga10_50_trim1  <- save_iso_output(iso_sga10_50_trim1,  "SGA 10rd %ile (.50, Trim1)",  "out_sga10_50_trim1.rda")




# ───────────────────────────────────────
# 📈 TRIMESTER 2 – SGA Models
# ───────────────────────────────────────

# ───────────────────────────────────────
# 📈 Fit Spline Models ----
# ───────────────────────────────────────
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/spline_binary_5_18.R")
spline_sga10_5_18_trim2 <- knot_fun(data = df_inf_sga10_trim2,
                                   hb_var = "hb",
                                   outcome_var = "sga10")
saveRDS(spline_sga10_5_18_trim2, "iso_results/spline_sga10_trim2.rds")

source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/spline_binary_boot.R")
spline_sga10_boot_trim2 <- knot_fun(df_inf_sga10_trim2, "hb", "sga10")
spline_sga10_bmodel_trim2 <- spline_sga10_boot_trim2$model
spline_sga10_boot_ci_trim2 <- spline_sga10_boot_trim2$boot_ci

saveRDS(spline_sga10_bmodel_trim2, "iso_results/spline_sga10_bmodel_trim2.rds")
saveRDS(spline_sga10_boot_ci_trim2, "iso_results/spline_sga10_boot_ci_trim2.rds")
save(spline_sga10_boot_trim2, file = "iso_results/spline_sga10_boot_trim2.rda")

# ───────────────────────────────────────
# 📉 Fit Isotonic Models ----
# ───────────────────────────────────────
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/iso_binary.R")
iso_sga10_trim2 <- flexstepreg_glmer(df_inf_sga10_trim2$sga10, df_inf_sga10_trim2$hb, df_inf_sga10_trim2$SITE, 
                                    covar2 = NULL, random_effect = NULL, alpha = 0.01)
saveRDS(iso_sga10_trim2, "iso_results/iso_sga10_trim2.rds")

source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/iso_binary_25.R")
iso_sga10_25_trim2 <- flexstepreg_glmer(df_inf_sga10_trim2$sga10, df_inf_sga10_trim2$hb, df_inf_sga10_trim2$SITE, 
                                       covar2 = NULL, random_effect = NULL, alpha = 0.01)
saveRDS(iso_sga10_25_trim2, "iso_results/iso_sga10_25_trim2.rds")

source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/iso_binary_50.R")
iso_sga10_50_trim2 <- flexstepreg_glmer(df_inf_sga10_trim2$sga10, df_inf_sga10_trim2$hb, df_inf_sga10_trim2$SITE, 
                                       covar2 = NULL, random_effect = NULL, alpha = 0.01)
saveRDS(iso_sga10_50_trim2, "iso_results/iso_sga10_50_trim2.rds")

# ───────────────────────────────────────
# 📊 Generate Plots ----
# ───────────────────────────────────────
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/plot_binary.R")
isoplot_sga10_trim2 <- iso_fun_glmer(
  x = df_inf_sga10_trim2$hb,
  y = df_inf_sga10_trim2$sga10,
  ylab = "SGA <10rd %ile",
  spline_model = spline_sga10_5_18_trim2,
  iso_model = iso_sga10_50_trim2
)

source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/plot_violin.R")
vioplot_sga10_trim2 <- spline_and_violin_plot(
  x = df_inf_sga10_trim2$hb,
  y = df_inf_sga10_trim2$sga10,
  xlab = "Hb",
  ylab = "SGA 10rd %ile",
  spline_model = spline_sga10_5_18_trim2,
  iso_model = iso_sga10_trim2,
  outcome_var = df_inf_sga10_trim2$sga10,
  title = "SGA 10rd %ile (Trimester 1)"
)

# ───────────────────────────────────────
# 💾 Save Plots ----
# ───────────────────────────────────────
png("iso_results/isoplot_sga10_trim2.png")
print(isoplot_sga10_trim2)
dev.off()

png("iso_results/vioplot_sga10_trim2.png")
print(vioplot_sga10_trim2)
dev.off()

# ───────────────────────────────────────
# 💾 Save Processed Isotonic Output Objects ----
# ───────────────────────────────────────
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/outdata.R")

save_iso_output <- function(model_obj, label, filename) {
  out <- outdata(model_obj, label)
  save(out, file = file.path("iso_results", filename))
  return(out)
}

out_sga10_trim2     <- save_iso_output(iso_sga10_trim2,     "SGA 10rd %ile (Trim2)",       "out_sga10_trim2.rda")
out_sga10_25_trim2  <- save_iso_output(iso_sga10_25_trim2,  "SGA 10rd %ile (.25, Trim2)",  "out_sga10_25_trim2.rda")
out_sga10_50_trim2  <- save_iso_output(iso_sga10_50_trim2,  "SGA 10rd %ile (.50, Trim2)",  "out_sga10_50_trim2.rda")



# ───────────────────────────────────────
# 📈 TRIMESTER 10 – SGA Models
# ───────────────────────────────────────

# ───────────────────────────────────────
# 📈 Fit Spline Models ----
# ───────────────────────────────────────
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/spline_binary_5_18.R")
spline_sga10_5_18_trim3 <- knot_fun(data = df_inf_sga10_trim3,
                                   hb_var = "hb",
                                   outcome_var = "sga10")
saveRDS(spline_sga10_5_18_trim3, "iso_results/spline_sga10_trim3.rds")

source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/spline_binary_boot.R")
spline_sga10_boot_trim3 <- knot_fun(df_inf_sga10_trim3, "hb", "sga10")
spline_sga10_bmodel_trim3 <- spline_sga10_boot_trim3$model
spline_sga10_boot_ci_trim3 <- spline_sga10_boot_trim3$boot_ci

saveRDS(spline_sga10_bmodel_trim3, "iso_results/spline_sga10_bmodel_trim3.rds")
saveRDS(spline_sga10_boot_ci_trim3, "iso_results/spline_sga10_boot_ci_trim3.rds")
save(spline_sga10_boot_trim3, file = "iso_results/spline_sga10_boot_trim3.rda")

# ───────────────────────────────────────
# 📉 Fit Isotonic Models ----
# ───────────────────────────────────────
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/iso_binary.R")
iso_sga10_trim3 <- flexstepreg_glmer(df_inf_sga10_trim3$sga10, df_inf_sga10_trim3$hb, df_inf_sga10_trim3$SITE, 
                                    covar2 = NULL, random_effect = NULL, alpha = 0.01)
saveRDS(iso_sga10_trim3, "iso_results/iso_sga10_trim3.rds")

source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/iso_binary_25.R")
iso_sga10_25_trim3 <- flexstepreg_glmer(df_inf_sga10_trim3$sga10, df_inf_sga10_trim3$hb, df_inf_sga10_trim3$SITE, 
                                       covar2 = NULL, random_effect = NULL, alpha = 0.01)
saveRDS(iso_sga10_25_trim3, "iso_results/iso_sga10_25_trim3.rds")

source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/iso_binary_50.R")
iso_sga10_50_trim3 <- flexstepreg_glmer(df_inf_sga10_trim3$sga10, df_inf_sga10_trim3$hb, df_inf_sga10_trim3$SITE, 
                                       covar2 = NULL, random_effect = NULL, alpha = 0.01)
saveRDS(iso_sga10_50_trim3, "iso_results/iso_sga10_50_trim3.rds")

# ───────────────────────────────────────
# 📊 Generate Plots ----
# ───────────────────────────────────────
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/plot_binary.R")
isoplot_sga10_trim3 <- iso_fun_glmer(
  x = df_inf_sga10_trim3$hb,
  y = df_inf_sga10_trim3$sga10,
  ylab = "SGA <10rd %ile",
  spline_model = spline_sga10_5_18_trim3,
  iso_model = iso_sga10_50_trim3
)

source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/plot_violin.R")
vioplot_sga10_trim3 <- spline_and_violin_plot(
  x = df_inf_sga10_trim3$hb,
  y = df_inf_sga10_trim3$sga10,
  xlab = "Hb",
  ylab = "SGA 10rd %ile",
  spline_model = spline_sga10_5_18_trim3,
  iso_model = iso_sga10_trim3,
  outcome_var = df_inf_sga10_trim3$sga10,
  title = "SGA 10rd %ile (Trimester 1)"
)

# ───────────────────────────────────────
# 💾 Save Plots ----
# ───────────────────────────────────────
png("iso_results/isoplot_sga10_trim3.png")
print(isoplot_sga10_trim3)
dev.off()

png("iso_results/vioplot_sga10_trim3.png")
print(vioplot_sga10_trim3)
dev.off()

# ───────────────────────────────────────
# 💾 Save Processed Isotonic Output Objects ----
# ───────────────────────────────────────
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/outdata.R")

save_iso_output <- function(model_obj, label, filename) {
  out <- outdata(model_obj, label)
  save(out, file = file.path("iso_results", filename))
  return(out)
}

out_sga10_trim3     <- save_iso_output(iso_sga10_trim3,     "SGA 10rd %ile (Trim3)",       "out_sga10_trim3.rda")
out_sga10_25_trim3  <- save_iso_output(iso_sga10_25_trim3,  "SGA 10rd %ile (.25, Trim3)",  "out_sga10_25_trim3.rda")
out_sga10_50_trim3  <- save_iso_output(iso_sga10_50_trim3,  "SGA 10rd %ile (.50, Trim3)",  "out_sga10_50_trim3.rda")


