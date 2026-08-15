#*****************************************************************************
# ReMAPP Aim 1: Fractional polynomial regression ----
#*****************************************************************************

# This script reproduces the primary ReMAPP Aim 1 analysis using the
# publicly available MAT_REMAPP_AIM1.csv dataset.

#*****************************************************************************
# 1. Packages ----
#*****************************************************************************

library(dplyr)
library(tidyr)
library(stringr)
library(lubridate)
library(gamlss)
library(tidyverse)
library(naniar)
library(flextable)
library(officer)
library(scales) 
library(wesanderson)



#*****************************************************************************
# 2. Set folders and load analysis dataset ----
#*****************************************************************************

## Data upload date
UploadDate <- "2026-01-30"

## TNT outcome data
path_to_tnt <- file.path(
  "Z:/Outcome Data",
  UploadDate
)

## Main ReMAPP analysis directory
# FOLDER SETUP (READ THIS BEFORE RUNNING!!!)
#Create your version!
# Before running this script, choose a location on your computer where the
# ReMAPP analysis files and results should be stored.
#
# The script will use the following folder structure:
#
# ReMAPP/
# └── Aim1/
#     └── results/
# Update base_dir below to the location where you would like the ReMAPP
# analysis folder to be stored.
#
# Examples:
# Windows: "C:/Users/your_name/Documents/ReMAPP/"
# Mac:     "~/Documents/ReMAPP/"

base_dir <- "D:/Users/williams_pj/Documents/Analysis/ReMAPP/"

## Working directory
working_dir <- file.path(base_dir, "Aim1")

## Subfolders
results_path <- file.path(working_dir, "results")

## Create testing folders if they do not already exist
dir.create(
  working_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

dir.create(
  results_path,
  showWarnings = FALSE,
  recursive = TRUE
)

# Set working directory
setwd(working_dir)

#Load remapp df
remapp_df <-  read.csv(
  file.path(path_to_tnt, "MAT_REMAPP_AIM1.csv"),
  stringsAsFactors = FALSE
)


#*****************************************************************************
# 3. Convert hemoglobin data to long format ----
#*****************************************************************************

# HB_LIST contains all hemoglobin measurements for each pregnancy in the
# following format:
#
# YYYY-MM-DD (gestational_age_weeks w): hemoglobin
#
# Multiple measurements are separated by semicolons.

reconstruct_hb_long <- function(data) {
  
  data %>%
    filter(
      !is.na(HB_LIST),
      HB_LIST != "",
      HB_LIST != "No Hb",
      HEALTHY_ELIGIBLE_G6PD == 1 #filter only analysis dataset 
    ) %>%
    
    # One row per hemoglobin measurement
    separate_rows(
      HB_LIST,
      sep = ";\\s*"
    ) %>%
    
    mutate(
      HB_LIST = str_trim(HB_LIST),
      
      # Hemoglobin measurement date
      HB_DATE = ymd(
        str_match(
          HB_LIST,
          "^([0-9]{4}-[0-9]{2}-[0-9]{2})"
        )[, 2]
      ),
      
      # Gestational age at hemoglobin measurement
      ga_wks = as.numeric(
        str_match(
          HB_LIST,
          "\\(([^w]+)w\\):"
        )[, 2]
      ),
      
      # Hemoglobin concentration in g/dL
      hb = as.numeric(
        str_match(
          HB_LIST,
          ":\\s*([0-9.]+)$"
        )[, 2]
      ),
      
      # Gestational age in days
      ga_days_hb = round(ga_wks * 7),
      
      # Trimester at measurement
      trimester = case_when(
        ga_wks > 0  & ga_wks < 14  ~ 1L,
        ga_wks >= 14 & ga_wks < 28 ~ 2L,
        ga_wks >= 28               ~ 3L,
        TRUE                        ~ NA_integer_
      )
    ) %>%
    
    select(
      -HB_LIST,
      -HB_N
    )
}


remapp_df_long <- reconstruct_hb_long(remapp_df)



#*****************************************************************************
# 4. Prepare FPR analytical dataset ----
#*****************************************************************************

df_fpr <- remapp_df_long %>%
  select(
    SITE,
    MOMID,
    PREGID,
    ga_wks,
    hb
  ) %>%
  filter(
    !is.na(hb),
    !is.na(ga_wks),
    ga_wks >= 6,
    ga_wks <= 40
  )


#*****************************************************************************
# 5. Fit fractional polynomial regression ----
#*****************************************************************************

# At gestational age t:
#
# Hb(t) ~ N(mu(t), sigma(t)^2)
#
# Second-order fractional polynomial models are fitted for both
# the mean and standard deviation.

fpr_control <- gamlss.control(
  trace = FALSE
)

fpr_model <- gamlss(
  hb ~ fp(ga_wks, 2),
  sigma.fo = ~ fp(ga_wks, 2),
  data = df_fpr,
  control = fpr_control
)


#*****************************************************************************
# 6. Extract model parameters ----
#*****************************************************************************

fpr_parameters <- data.frame(
  Parameter = c("Mean (mu)", "Standard deviation (sigma)"),
  
  Intercept = c(
    fpr_model$mu.coefficients[[1]],
    fpr_model$sigma.coefficients[[1]]
  ),
  
  FP_Intercept = c(
    getSmo(fpr_model, what = "mu")$coef[[1]],
    getSmo(fpr_model, what = "sigma")$coef[[1]]
  ),
  
  Coefficient_1 = c(
    getSmo(fpr_model, what = "mu")$coef[[2]],
    getSmo(fpr_model, what = "sigma")$coef[[2]]
  ),
  
  Coefficient_2 = c(
    getSmo(fpr_model, what = "mu")$coef[[3]],
    getSmo(fpr_model, what = "sigma")$coef[[3]]
  ),
  
  Power_1 = c(
    getSmo(fpr_model, what = "mu")$power[[1]],
    getSmo(fpr_model, what = "sigma")$power[[1]]
  ),
  
  Power_2 = c(
    getSmo(fpr_model, what = "mu")$power[[2]],
    getSmo(fpr_model, what = "sigma")$power[[2]]
  ),
  
  Link = c(
    fpr_model$mu.link,
    fpr_model$sigma.link
  ),
  
  check.names = FALSE
)


#*****************************************************************************
# 7. Estimate weekly FPR centiles ----
#*****************************************************************************

# Generate predicted 5th, 50th, and 95th centiles for exact
# gestational weeks 6 through 40.

fpr_centiles <- centiles.pred(
  fpr_model,
  data = df_fpr,
  xname = "ga_wks",
  xvalues = 6:40,
  type = "centiles",
  cent = c(5, 50, 95)
) %>%
  as.data.frame() %>%
  rename(
    GA_WEEKS = x,
    CENTILE_5 = `5`,
    CENTILE_50 = `50`,
    CENTILE_95 = `95`
  )


#*****************************************************************************
# 8. Add gestational-week sample size ----
#*****************************************************************************

ga_sample_size <- df_fpr %>%
  mutate(
    GA_DAYS = round(ga_wks * 7),
    GA_WEEKS = floor(GA_DAYS / 7)
  ) %>%
  count(
    GA_WEEKS,
    name = "N"
  )


fpr_centiles <- fpr_centiles %>%
  left_join(
    ga_sample_size,
    by = "GA_WEEKS"
  ) %>%
  mutate(
    `Gestational age (exact week)` =
      paste0(GA_WEEKS, " weeks + 0 days"),
    
    `5th centile` = round(CENTILE_5, 1),
    `50th centile` = round(CENTILE_50, 1),
    `95th centile` = round(CENTILE_95, 1)
  ) %>%
  select(
    `Gestational age (exact week)`,
    `Sample size` = N,
    `5th centile`,
    `50th centile`,
    `95th centile`
  )


#*****************************************************************************
# 9. Summarise predicted centiles ----
#*****************************************************************************

# For each gestational period, report the median of the gestation-specific
# predicted 5th, 50th, and 95th centiles.

get_summary_hb_values_median <- function(centile_table) {
  
  weeks <- as.numeric(
    gsub(
      " weeks.*",
      "",
      centile_table$`Gestational age (exact week)`
    )
  )
  
  overall_summary <- data.frame(
    `Time period` = "Overall",
    `GA weeks` = "6-40",
    `5th centile (g/dL)` =
      round(median(centile_table$`5th centile`, na.rm = TRUE), 1),
    `50th centile (g/dL)` =
      round(median(centile_table$`50th centile`, na.rm = TRUE), 1),
    `95th centile (g/dL)` =
      round(median(centile_table$`95th centile`, na.rm = TRUE), 1),
    check.names = FALSE
  )
  
  first_trimester <- data.frame(
    `Time period` = "First trimester",
    `GA weeks` = "6-13",
    `5th centile (g/dL)` =
      round(median(
        centile_table$`5th centile`[weeks >= 6 & weeks <= 13],
        na.rm = TRUE
      ), 1),
    `50th centile (g/dL)` =
      round(median(
        centile_table$`50th centile`[weeks >= 6 & weeks <= 13],
        na.rm = TRUE
      ), 1),
    `95th centile (g/dL)` =
      round(median(
        centile_table$`95th centile`[weeks >= 6 & weeks <= 13],
        na.rm = TRUE
      ), 1),
    check.names = FALSE
  )
  
  second_trimester <- data.frame(
    `Time period` = "Second trimester",
    `GA weeks` = "14-27",
    `5th centile (g/dL)` =
      round(median(
        centile_table$`5th centile`[weeks >= 14 & weeks <= 27],
        na.rm = TRUE
      ), 1),
    `50th centile (g/dL)` =
      round(median(
        centile_table$`50th centile`[weeks >= 14 & weeks <= 27],
        na.rm = TRUE
      ), 1),
    `95th centile (g/dL)` =
      round(median(
        centile_table$`95th centile`[weeks >= 14 & weeks <= 27],
        na.rm = TRUE
      ), 1),
    check.names = FALSE
  )
  
  third_trimester <- data.frame(
    `Time period` = "Third trimester",
    `GA weeks` = "28-40",
    `5th centile (g/dL)` =
      round(median(
        centile_table$`5th centile`[weeks >= 28 & weeks <= 40],
        na.rm = TRUE
      ), 1),
    `50th centile (g/dL)` =
      round(median(
        centile_table$`50th centile`[weeks >= 28 & weeks <= 40],
        na.rm = TRUE
      ), 1),
    `95th centile (g/dL)` =
      round(median(
        centile_table$`95th centile`[weeks >= 28 & weeks <= 40],
        na.rm = TRUE
      ), 1),
    check.names = FALSE
  )
  
  bind_rows(
    overall_summary,
    first_trimester,
    second_trimester,
    third_trimester
  )
}


hb_summary_median <- get_summary_hb_values_median(
  fpr_centiles
)


#*****************************************************************************
# 10. WHO-style anemia thresholds ----
#*****************************************************************************

# WHO-style anemia severity thresholds are derived from the model-estimated
# 5th centile.
#
# Mild anemia:     80% to <100% of the 5th centile
# Moderate anemia: 60% to <80% of the 5th centile
# Severe anemia:   <60% of the 5th centile

create_who_anemia_table <- function(centile_summary) {
  
  centile_summary %>%
    transmute(
      `Time period`,
      `GA weeks`,
      
      `5th centile Hb (g/dL)` =
        `5th centile (g/dL)`,
      
      `95th centile Hb (g/dL)` =
        `95th centile (g/dL)`,
      
      `Mild anemia` = sprintf(
        "%.1f-<%.1f g/dL",
        `5th centile (g/dL)` * 0.80,
        `5th centile (g/dL)`
      ),
      
      `Moderate anemia` = sprintf(
        "%.1f-<%.1f g/dL",
        `5th centile (g/dL)` * 0.60,
        `5th centile (g/dL)` * 0.80
      ),
      
      `Severe anemia` = paste0(
        "<",
        sprintf(
          "%.1f",
          `5th centile (g/dL)` * 0.60
        ),
        " g/dL"
      )
    )
}


hc_anemia_who <- create_who_anemia_table(
  hb_summary_median
)


#*****************************************************************************
# 11. Plot FPR centile curves ----
#*****************************************************************************

# Retain the plotting colors and appearance used in the primary
# manuscript analysis:
#
# Observed Hb = grey86
# 5th centile = red
# 50th centile = green
# 95th centile = blue
# WHO threshold = black dashed line

png(
  filename = "results/FPR_centile_curves.png",
  width = 1800,
  height = 1400,
  res = 200
)

centiles(
  fpr_model,
  
  xvar = df_fpr$ga_wks,
  
  cent = c(
    5,
    50,
    95
  ),
  
  legend = FALSE,
  
  xlab = "Gestational Age (weeks)",
  ylab = "CBC HB (g/dL)",
  
  main = paste0(
    "Centile Curves Using FPR Model\n",
    "(5, 50, 95)th Centiles vs WHO Threshold\n",
    "Pooled"
  ),
  
  bg = "transparent",
  
  xlim = c(
    6,
    40
  ),
  
  ylim = c(
    5,
    18
  ),
  
  pch = 15,
  
  cex = 0.2,
  
  col = "grey86",
  
  col.cent = c(
    "red",
    "green",
    "blue"
  )
)


# WHO reference thresholds

segments(
  x0 = 6,
  y0 = 11,
  x1 = 14,
  y1 = 11,
  col = "black",
  lty = 2
)

segments(
  x0 = 14,
  y0 = 10.5,
  x1 = 27,
  y1 = 10.5,
  col = "black",
  lty = 2
)

segments(
  x0 = 27,
  y0 = 11,
  x1 = 36,
  y1 = 11,
  col = "black",
  lty = 2
)


# WHO threshold labels

text(
  x = 6.5,
  y = 11.5,
  "11",
  col = "black",
  cex = 0.66
)

text(
  x = 20.5,
  y = 11,
  "10.5",
  col = "black",
  cex = 0.66
)

text(
  x = 34.5,
  y = 11.5,
  "11",
  col = "black",
  cex = 0.66
)


# Figure legend

legend(
  "topright",
  
  inset = c(
    0,
    0
  ),
  
  legend = c(
    "5th centile",
    "50th centile",
    "95th centile",
    "WHO threshold"
  ),
  
  cex = 0.6,
  
  col = c(
    "red",
    "green",
    "blue",
    "black"
  ),
  
  lty = c(
    1,
    1,
    1,
    2
  )
)

dev.off()


#*****************************************************************************
# 12. Create results workbook ----
#*****************************************************************************

# Save all numerical results in one workbook.
# The FPR figure is saved separately as a PNG.

results_workbook <- createWorkbook()


# Model parameters ----

addWorksheet(
  results_workbook,
  "Model Parameters"
)

writeData(
  results_workbook,
  sheet = "Model Parameters",
  x = fpr_parameters,
  withFilter = TRUE
)


# Weekly centiles ----

addWorksheet(
  results_workbook,
  "Weekly Centiles"
)

writeData(
  results_workbook,
  sheet = "Weekly Centiles",
  x = fpr_centiles,
  withFilter = TRUE
)


# Centile summary ----

addWorksheet(
  results_workbook,
  "Centile Summary"
)

writeData(
  results_workbook,
  sheet = "Centile Summary",
  x = hb_summary_median,
  withFilter = TRUE
)


# WHO thresholds ----

addWorksheet(
  results_workbook,
  "WHO Thresholds"
)

writeData(
  results_workbook,
  sheet = "WHO Thresholds",
  x = hc_anemia_who,
  withFilter = TRUE
)


#*****************************************************************************
# 13. Format results workbook ----
#*****************************************************************************

header_style <- createStyle(
  textDecoration = "bold",
  halign = "center",
  valign = "center",
  border = "Bottom"
)

for (
  sheet_name in names(results_workbook)
) {
  
  addStyle(
    results_workbook,
    sheet = sheet_name,
    style = header_style,
    rows = 1,
    cols = 1:ncol(
      switch(
        sheet_name,
        "Model Parameters" = fpr_parameters,
        "Weekly Centiles" = fpr_centiles,
        "Centile Summary" = hb_summary_median,
        "WHO Thresholds" = hc_anemia_who
      )
    ),
    gridExpand = TRUE
  )
  
  freezePane(
    results_workbook,
    sheet = sheet_name,
    firstRow = TRUE
  )
  
  setColWidths(
    results_workbook,
    sheet = sheet_name,
    cols = 1:ncol(
      switch(
        sheet_name,
        "Model Parameters" = fpr_parameters,
        "Weekly Centiles" = fpr_centiles,
        "Centile Summary" = hb_summary_median,
        "WHO Thresholds" = hc_anemia_who
      )
    ),
    widths = "auto"
  )
}


#*****************************************************************************
# 14. Save results ----
#*****************************************************************************

saveWorkbook(
  results_workbook,
  file = "results/ReMAPP_Aim1_FPR_Results.xlsx",
  overwrite = TRUE
)

saveRDS(
  fpr_model,
  file = "results/FPR_model.rds"
)


#*****************************************************************************
# 15. Analysis summary ----
#*****************************************************************************

cat(
  "\nReMAPP Aim 1 FPR analysis complete\n",
  "----------------------------------\n",
  "Participants: ",
  n_distinct(df_fpr$MOMID),
  "\n",
  "Pregnancies: ",
  n_distinct(df_fpr$PREGID),
  "\n",
  "Hemoglobin measurements: ",
  nrow(df_fpr),
  "\n",
  "Gestational age range: ",
  round(min(df_fpr$ga_wks), 1),
  "-",
  round(max(df_fpr$ga_wks), 1),
  " weeks\n\n",
  "Outputs:\n",
  "results/ReMAPP_Aim1_FPR_Results.xlsx\n",
  "results/FPR_centile_curves.png\n",
  "results/FPR_model.rds\n"
)