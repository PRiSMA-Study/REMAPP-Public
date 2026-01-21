

# Packages and cleaning codes ----
suppressPackageStartupMessages({
  library(dplyr); library(stringr); library(tidyr); library(ggplot2); library(ggrepel); library(forcats); library(purrr)
})

## outcome df cleaning ----
clean_outcome <- function(x) {
  x %>%
    str_replace("\\s*\\(\\s*\\.50\\s*(,\\s*Trim\\s*[123])?\\s*\\)\\s*$", "") %>% 
    str_squish()
}

## map raw grid names ----
relabel_outcome <- function(x) {
  x <- toupper(x)
  dplyr::recode(x,
                "COMPO"                = "Infant Composite",
                "PRETERM 37WKS"        = "Preterm<37weeks",
                "PRETERM 34WKS"        = "Preterm<34weeks",
                "SGA10"                = "SGA<10th",
                "SGA3"                 = "SGA<3th",
                "LBW<2500G"            = "LBW<2500g",
                "LBW<1500G"            = "LBW<1500g",
                "STILLBIRTH 20WKS"     = "Stillbirth>=20weeks",
                "STILLBIRTH 28WKS"     = "Stillbirth>=28weeks",
                "HYPERBILIRUBINEMIA"   = "Hyperbilirubinemia",
                "NEONATAL MORTALITY"   = "Neontal Mortality",
                "PSBI"                 = "PSBI",
                "ASPH"                 = "ASPH",
                "MAT COMPO"            = "Maternal Composite",
                "PRECLAMP"             = "Preeclampsia",
                "PREECLAMP"            = "Preeclampsia",
                "PREECLAMPSIA"         = "Preeclampsia",
                "PPROM"                = "PPROM",
                "PPH"                  = "PPH",
                "PPA PNC6"             = "PPA-PNC6",
                "PPA PNC26"            = "PPA-PNC26",
                .default = x
  )
}

## axis breaks using df_hb_long2 if available; else from hb_grid
5 <- tryCatch(floor(min(df_hb_long2$hb, na.rm = TRUE)), error = function(e) floor(min(hb_grid$hb, na.rm = TRUE)))
18 <- tryCatch(ceiling(max(df_hb_long2$hb, na.rm = TRUE)), error = function(e) ceiling(max(hb_grid$hb, na.rm = TRUE)))

## Build df_new from hb_grid to mimic your old structure ----
df_new <- hb_grid %>%
  transmute(
    Outcome_raw = clean_outcome(Outcome),
    Outcome_lbl = relabel_outcome(Outcome_raw),
    Hb_sequence = hb,
    risk = risk
  ) %>%
  # drop NA risks and EPDS/Fatigue score outcomes
  filter(!is.na(risk)) %>%
  filter(!str_detect(Outcome_lbl, regex("EPDS\\s*score|FATIGUE\\s*score|FAT|EPDS|DPR", ignore_case = TRUE))) %>%
  rename(Outcome = Outcome_lbl)

# All outcomes plot ----

## palette & ordering -----
colors_20_custom <- c(
  "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728",
  "#9467bd", "#8c564b", "#e377c2", "#7f7f7f",
  "#bcbd22", "#17becf", "#9E9D24", "#f39c12",
  "#2ecc71", "#e74c3c", "#2196F3", "#1abc9c",
  "#8E44AD", "#34495e", "#16a085", "#f1c40f",
  "#2d3987","#C52C1A"
)

if (length(unique(df_new$Outcome)) > length(colors_20_custom)) {
  stop("Not enough colors for the number of unique outcomes!")
}

## order outcomes by risk at max Hb so labels stack nicely at the right
end_hb <- max(df_new$Hb_sequence, na.rm = TRUE)
ord <- df_new %>%
  filter(Hb_sequence == end_hb) %>%
  group_by(Outcome) %>%
  summarise(risk_end = mean(risk, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(risk_end)) %>% pull(Outcome)
df_new <- df_new %>% mutate(Outcome = factor(Outcome, levels = ord))

## create the plot
ggplot(df_new, aes(x = Hb_sequence, y = risk, color = Outcome)) +
  geom_line() +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) +
  scale_x_continuous(   breaks = seq(5, 18, by = 1),   
                        minor_breaks = seq(5, 18, by = 0.5),   
                        limits = c(5, 18) )+
  labs(
    title = "Overview of Risk of Events vs Hemoglobin",
    x = "Hemoglobin",
    y = "Risk of Event",
    color = "Outcome"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 1, hjust = 0.5),
    axis.text = element_text(size = 8),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(size = 8),
    plot.title = element_text(size = 12)
  ) +
  scale_color_manual(values = colors_20_custom) +
  geom_label_repel(
    data = df_new %>%
      group_by(Outcome) %>%
      filter(Hb_sequence == max(Hb_sequence, na.rm = TRUE)),
    aes(label = Outcome, color = Outcome),
    box.padding = 0.25,
    point.padding = 0.25,
    segment.color = "grey50",
    size = 2,
    max.overlaps = 13
  )

plot_group_1 <- function(df_cat, title_text) {
  ggplot(df_cat, aes(x = Hb_sequence, y = risk, color = Outcome)) +
    geom_line() +
    scale_y_continuous(breaks = seq(0, 1, 0.1)) +
    scale_x_continuous(   breaks = seq(5, 18, by = 1),   
                          minor_breaks = seq(5, 18, by = 0.5),   
                          limits = c(5, 18) ) +
    labs(
      title = title_text,
      x = "Hemoglobin",
      y = "Risk of Event",
      color = "Outcome"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 1, hjust = 0.5),
      axis.text = element_text(size = 8),
      legend.position = "none",
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "white"),
      strip.text = element_text(size = 8),
      plot.title = element_text(size = 12)
    ) +
    scale_color_manual(values = colors_20_custom) +
    geom_label_repel(
      data = df_cat %>%
        group_by(Outcome) %>%
        filter(Hb_sequence == max(Hb_sequence, na.rm = TRUE)),
      aes(label = Outcome, color = Outcome),
      box.padding = 0.25,
      point.padding = 0.25,
      segment.color = "grey50",
      size = 2,
      max.overlaps = 0
    )
}

plot_group_2 <- function(df_cat, title_text) {
  ggplot(df_cat, aes(x = Hb_sequence, y = risk, color = Outcome)) +
    geom_line() +
    scale_y_continuous(breaks = seq(0, 1, 0.05)) +
    scale_x_continuous(   breaks = seq(5, 18, by = 1),   
                          minor_breaks = seq(5, 18, by = 0.5),   
                          limits = c(5, 18) ) +
    labs(
      title = title_text,
      x = "Hemoglobin",
      y = "Risk of Event",
      color = "Outcome"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 1, hjust = 0.5),
      axis.text = element_text(size = 8),
      legend.position = "none",
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "white"),
      strip.text = element_text(size = 8),
      plot.title = element_text(size = 12)
    ) +
    scale_color_manual(values = colors_20_custom) +
    geom_label_repel(
      data = df_cat %>%
        group_by(Outcome) %>%
        filter(Hb_sequence == max(Hb_sequence, na.rm = TRUE)),
      aes(label = Outcome, color = Outcome),
      box.padding = 0.25,
      point.padding = 0.25,
      segment.color = "grey50",
      size = 2,
      max.overlaps = 0
    )
}
# Subset 1: Infant (Composite & Size @ Birth + Preterm) ----
df_cat1 <- df_new %>%
  filter(Outcome %in% c("Infant Composite","LBW<2500g","LBW<1500g",
                        "SGA<10th","SGA<3th","Preterm<34weeks","Preterm<37weeks"))

plot_group_1(df_cat1, "Risk vs Hemoglobin - Infant Birth Outcomes")

# Subset 2: Birth outcomes (Stillbirths, Hyperbili, ASPH, PSBI) ----
df_cat2 <- df_new %>%
  filter(Outcome %in% c("Stillbirth>=20weeks","Stillbirth>=28weeks",
                        "Hyperbilirubinemia","ASPH","PSBI"))

plot_group_2(df_cat2, "Risk vs Hemoglobin - Infant Birth Timing and Health COnditions")


# Subset 3: Maternal complications (PPROM, Preeclampsia, PPH) ----
df_cat3 <- df_new %>%
  filter(Outcome %in% c("Maternal Composite","PPROM","Preeclampsia","PPH"))

plot_group_2(df_cat3, "Risk vs Hemoglobin - Maternal Complications")

# Subset 4: Postpartum anemia (6w & 26w) ----
df_cat4 <- df_new %>%
  filter(Outcome %in% c("PPA-PNC26","PPA-PNC6")) %>%
  mutate(across(starts_with("Risk"), ~ round(., 1), .names = "{.col}"))

plot_group_1(df_cat4, "Risk vs Hemoglobin - Postpartum Anemia")

#Subset 5: Depression ----
df_dpr <- hb_grid %>%
  transmute(
    Outcome_raw = clean_outcome(Outcome),
    Outcome_lbl = relabel_outcome(Outcome_raw),
    Hb_sequence = hb,
    risk = risk
  ) %>%
  # drop NA risks and EPDS/Fatigue score outcomes
  filter(!is.na(risk)) %>%
  filter(str_detect(Outcome_lbl, regex("DPR\\s*score|EPDS|DPR", ignore_case = TRUE))) %>%
  rename(Outcome = Outcome_lbl)

ggplot(df_dpr, aes(x = Hb_sequence, y = risk, color = Outcome)) +
  geom_line() +
  #scale_y_continuous(breaks = seq(0, 1, 0.1)) +
  scale_x_continuous(   breaks = seq(5, 18, by = 1),   
                        minor_breaks = seq(5, 18, by = 0.5),   
                        limits = c(5, 18) )+
  labs(
    title = "Overview of Depression Score (EPDS) vs Hemoglobin",
    x = "Hemoglobin",
    y = "Score",
    color = "Outcome"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 1, hjust = 0.5),
    axis.text = element_text(size = 8),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(size = 8),
    plot.title = element_text(size = 12)
  ) +
  scale_color_manual(values = colors_20_custom) +
  geom_label_repel(
    data = df_dpr %>%
      group_by(Outcome) %>%
      slice_max(Hb_sequence, n = 1, with_ties = FALSE) %>%
      ungroup(),
    aes(label = Outcome, color = Outcome),
    box.padding = 0.25,
    point.padding = 0.25,
    segment.color = "grey50",
    size = 2
  )

#Subset 6: Fatigue ----
df_fat <- hb_grid %>%
  transmute(
    Outcome_raw = clean_outcome(Outcome),
    Outcome_lbl = relabel_outcome(Outcome_raw),
    Hb_sequence = hb,
    risk = risk
  ) %>%
  # drop NA risks and EPDS/Fatigue score outcomes
  filter(!is.na(risk)) %>%
  filter(str_detect(Outcome_lbl, regex("FAT\\s*score|FATIGUE\\s*score|FAT", ignore_case = TRUE))) %>%
  rename(Outcome = Outcome_lbl)

ggplot(df_fat, aes(x = Hb_sequence, y = risk, color = Outcome)) +
  geom_line() +
  #scale_y_continuous(breaks = seq(0, 1, 0.1)) +
  scale_x_continuous(   breaks = seq(5, 18, by = 1),   
                        minor_breaks = seq(5, 18, by = 0.5),   
                        limits = c(5, 18) )+
  labs(
    title = "Overview of Fatigue Score vs Hemoglobin",
    x = "Hemoglobin",
    y = "Score",
    color = "Outcome"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 1, hjust = 0.5),
    axis.text = element_text(size = 8),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(size = 8),
    plot.title = element_text(size = 12)
  ) +
  scale_color_manual(values = colors_20_custom) +
  geom_label_repel(
    data = df_fat %>%
      group_by(Outcome) %>%
      slice_max(Hb_sequence, n = 1, with_ties = FALSE) %>%
      ungroup(),
    aes(label = Outcome, color = Outcome),
    box.padding = 0.25,
    point.padding = 0.25,
    segment.color = "grey50",
    size = 2
  )
