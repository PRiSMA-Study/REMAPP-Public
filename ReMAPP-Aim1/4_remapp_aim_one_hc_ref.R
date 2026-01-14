REMAPP <- all_screen_df %>%
  filter(screen_pass == "Yes") %>%
  select(MOMID, PREGID, SITE, HEALTHY_ELIGIBLE, HB_DATA, HB_VALID_DATA)

# Dichotomous vars
dichotomous_vars <- c("MARRIED", "PAID_WORK", "NULLIPAROUS", 
                      "MISCARRIAGE", "LIVEBIRTH", "SEX",
                      "PRETERMBIRTH_LT37", "LBW2500_ANY", "INF_DTH")

# Merge
demo_df <- REMAPP %>%
  left_join(inf_perinatal_df, by = c("MOMID", "PREGID", "SITE")) %>%
  left_join(MAT_DEMOGRAPHIC, by = c("MOMID", "PREGID", "SITE")) %>%
  mutate(
    # new grouping variable with two levels
    COHORT_CAT = case_when(
      HEALTHY_ELIGIBLE == 1 & HB_DATA == "Yes" & HB_VALID_DATA == "Yes" ~ "Healthy cohort",
      TRUE ~ "Full Cohort"
    ),
    # second flag to force Full Cohort group to include everyone
    ALL_FLAG = "Full Cohort",
    BWEIGHT_KG = round(BWEIGHT_ANY/1000, 1)
  )

# Stack rows for both groups
demo_df_long <- bind_rows(
  demo_df %>% mutate(GROUP = ALL_FLAG),
  demo_df %>% filter(COHORT_CAT == "Healthy cohort") %>% mutate(GROUP = "Healthy cohort")
)

# Now create table
demo_table <- demo_df_long %>%
  ungroup() %>%
  select(GROUP, MAT_AGE, GA_WKS_ENROLL, BMI_ENROLL, MUAC_ENROLL, SCHOOL_YRS,
         MARRIED, PAID_WORK, NULLIPAROUS, MISCARRIAGE) %>%
  tbl_summary(
    by = GROUP,
    type = list(
      all_of(c("MARRIED", "PAID_WORK", "NULLIPAROUS", "MISCARRIAGE")) ~ "dichotomous",
      all_continuous() ~ "continuous"
    ),
    value = list(
      all_of(c("MARRIED", "PAID_WORK", "NULLIPAROUS", "MISCARRIAGE")) ~ "1"
    ),
    label = list(
      MAT_AGE ~ "Age (years); mean (SD)",
      GA_WKS_ENROLL ~ "Gestational age at first visit (weeks); mean (SD)",
      BMI_ENROLL ~ "Body mass index (kg/m²); mean (SD)",
      MUAC_ENROLL ~ "MUAC at first visit (cm); mean (SD)",
      SCHOOL_YRS ~ "Years of formal education; mean (SD)",
      MARRIED ~ "Married or cohabiting; n (%)",
      PAID_WORK ~ "Engaged in paid work; n (%)",
      NULLIPAROUS ~ "Never gave birth before; n (%)",
      MISCARRIAGE ~ "Previous miscarriage; n (%)"
    ),
    missing = "no",
    statistic = list(
      all_continuous() ~ "{mean} ± {sd}",
      all_categorical() ~ "{p}%"
    ),
    digits = list(
      all_categorical() ~ 1
    )
  ) %>%
  modify_header(
    label ~ "**Demographic Characteristics**"
  ) %>%
  modify_caption(
    "**<div align='left'><b>Maternal Baseline Demographic Characteristics**"
  ) %>%
  modify_footnote_body(
    footnote = "MUAC = Mid-upper arm circumference",
    columns = "label",
    rows = variable == "MUAC_ENROLL"
  )

# Save & trim
demo_table %>%
  as_gt() %>%
  gtsave("baseline_demo.png")

img_base <- image_read("baseline_demo.png") %>% image_trim()


image_write(img_base, path = "baseline_demo_trimmed.png")
keep_vars <- c(
  # continuous
  "MAT_AGE","GA_WKS_ENROLL","SCHOOL_YRS",
  # candidate categorical/dichotomous
  "MARRIED","GPARITY_CAT",
  "SUBSTANCE_USE","MALARIA","HEP_B","HEP_C","IRON_DEF",
  "RBC_HEMO","G6PD","PREV_LBW","PREV_STILLBIRTH","PREV_UNPLANNED_CS",
  # multi-category
  "AGE_CAT","BMI_CAT","MUAC_CAT","HEIGHT_CAT"
)