#****************************************************************************
#EPDS score
#****************************************************************************
library(tidyverse)

source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/spline_continuous.R")
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/iso_continuous.R")
#****************************************************************************
#*******all
#prepare data
# Load the dataset
# Load data
load("derived_data/df_mat_dpr.rda")

# Prepare ANC data
df_dpr_anc <- df_mat_dpr %>%
  filter(visit %in% c("ANC20", "ANC32")) %>%
  filter(!is.na(hb) & !is.na(dpr_score)) %>%
  arrange(MOMID, visit)

# Prepare PNC data
df_dpr_pnc <- df_mat_dpr %>%
  filter(visit == "PNC6") %>%
  filter(!is.na(hb) & !is.na(dpr_score))

save(df_dpr_anc, file = "derived_data/df_dpr_anc.rda")
save(df_dpr_pnc, file = "derived_data/df_dpr_pnc.rda")


# ****************************************************************************
# ANC models (hemoglobin predicts depression during pregnancy)

# Spline model
spline_dpr_score <- knot_fun_continue(df_dpr_anc, "hb", "dpr_score")
saveRDS(spline_dpr_score, "iso_results/spline_dpr_score.rds")

# Isotonic model
iso_dpr_score <- flexstepreg_lmer(df_dpr_anc$dpr_score, df_dpr_anc$hb, df_dpr_anc$SITE,
                                  covar2 = df_dpr_anc$visit,
                                  random_effect = df_dpr_anc$MOMID,
                                  alpha = 0.01)
saveRDS(iso_dpr_score, "iso_results/iso_dpr_score.rds")

# ****************************************************************************
# PNC models (postpartum depression)

# Spline model
spline_dpr_score_pnc <- knot_fun_continue(df_dpr_pnc, "hb", "dpr_score")
saveRDS(spline_dpr_score_pnc, "iso_results/spline_dpr_score_pnc.rds")

# Isotonic model
iso_dpr_score_pnc <- flexstepreg_lmer(df_dpr_pnc$dpr_score, df_dpr_pnc$hb, df_dpr_pnc$SITE,
                                  covar2 =NULL,
                                  random_effect = df_dpr_pnc$MOMID,
                                  alpha = 0.01)
saveRDS(iso_dpr_score_pnc, "iso_results/iso_dpr_score_pnc.rds")

# # Overall counts
# summary_table <- df_mat_dpr %>%
#   group_by(visit) %>%
#   summarise(
#     total_records = n(),
#     hb_present = sum(!is.na(hb)),
#     dpr_score_present = sum(!is.na(dpr_score)),
#     complete_pair = sum(!is.na(hb) & !is.na(dpr_score))
#   )
# 
# print(summary_table)