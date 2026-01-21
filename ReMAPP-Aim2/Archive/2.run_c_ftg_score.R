#****************************************************************************
#Fatigue score
#****************************************************************************
library(tidyverse)

source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/spline_continuous.R")
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/iso_continuous.R")


# Load data
load("derived_data/df_mat_ftg.rda")

# Prepare ANC data
df_ftg_anc <- df_mat_ftg %>%
  filter(visit %in% c("ANC20", "ANC32")) %>%
  filter(!is.na(hb) & !is.na(ftg_score)) %>%
  arrange(MOMID, visit)


# Prepare PNC data
df_ftg_pnc <- df_mat_ftg %>%
  filter(visit == "PNC6") %>%
  filter(!is.na(hb) & !is.na(ftg_score))

save(df_ftg_anc, file = "derived_data/df_ftg_anc.rda")
save(df_ftg_pnc, file = "derived_data/df_ftg_pnc.rda")

# ****************************************************************************
# ANC models (hemoglobin predicts depression during pregnancy)

# Spline model
spline_ftg_score <- knot_fun_continue(df_ftg_anc, "hb", "ftg_score")
saveRDS(spline_ftg_score, "iso_results/spline_ftg_score.rds")

# Isotonic model
iso_ftg_score <- flexstepreg_lmer(df_ftg_anc$ftg_score, df_ftg_anc$hb, df_ftg_anc$SITE,
                                  covar2 = df_ftg_anc$visit,
                                  random_effect = df_ftg_anc$MOMID,
                                  alpha = 0.01)
saveRDS(iso_ftg_score, "iso_results/iso_ftg_score.rds")

# ****************************************************************************
# PNC models (postpartum depression)

# Spline model
spline_ftg_score_pnc <- knot_fun_continue(df_ftg_pnc, "hb", "ftg_score")
saveRDS(spline_ftg_score_pnc, "iso_results/spline_ftg_score_pnc.rds")

# Isotonic model
iso_ftg_score_pnc <- flexstepreg_lmer(df_ftg_pnc$ftg_score, df_ftg_pnc$hb, df_ftg_pnc$SITE,
                                      covar2 =NULL,
                                      random_effect = df_ftg_pnc$MOMID,
                                      alpha = 0.01)
saveRDS(iso_ftg_score_pnc, "iso_results/iso_ftg_score_pnc.rds")

