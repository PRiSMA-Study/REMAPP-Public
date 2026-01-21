\newpage 

## Summary Tables of Hemoglobin Thresholds by Risk Rank
```{r, echo=FALSE, include=FALSE}

# code -> pretty label
outcome_label_map <- c(
  # Infant
  "PRETERM 37WKS"        = "Preterm <37wks",
  "PRETERM 34WKS"        = "Preterm <34wks",
  "LBW<2500g"            = "LBW <2500g",
  "LBW<1500g"            = "LBW <1500g",
  "SGA10"                = "SGA <10th percentile",
  "SGA3"                 = "SGA <3rd percentile",
  "PSBI"                 = "Neonatal PSBI",
  "ASPH"                 = "Neonatal Asphyxia",
  "NEONATAL MORTALITY"   = "Neonatal Mortality",
  "STILLBIRTH 20WKS"     = "Stillbirth 20wks",
  "STILLBIRTH 28WKS"     = "Stillbirth 28wks",
  "HYPERBILIRUBINEMIA"   = "Hyperbilirubinemia",
  "COMPO"                = "Infant composite",
  
  # Maternal
  "PPH"                  = "Postpartum hemorrhage",
  "PPA PNC6"             = "PPA by PNC6",
  "PPA PNC26"            = "PPA by PNC26",
  "PPROM"                = "Premature Rupture of Membranes",
  "PRECLAMP"             = "Preeclampsia",
  "MAT COMPO"            = "Maternal Composite",
  "FATIGUE ANC"          = "Fatigued (ANC)",
  "FATIGUE PNC"          = "Fatigued (PNC)",
  "DPR ANC"              = "Depression (ANC)",
  "DPR PNC"              = "Depression (PNC)"
)


heatmap_logic_to_table <- function(hb_grid,
                                   trimester = c("All","Trim1","Trim2","Trim3","Pp6wks"),
                                   low_tail_cut = 8,
                                   high_tail_cut = 13,
                                   near_max_frac = 0.90,
                                   risk_digits = 2,
                                   drop_mono_risk = TRUE,
                                   drop_single_threshold = TRUE,
                                   exclude_outcomes = c("WHO Guideline","Healthy Cohort")) {
  stopifnot(all(c("Outcome","hb","risk") %in% names(hb_grid)))
  suppressPackageStartupMessages({
    library(dplyr); library(stringr); library(tidyr); library(tibble)
  })
  trimester <- match.arg(trimester)
  tol <- 1e-9
  
  # 1) Parse trimester & clean outcome names
  df0 <- hb_grid %>%
    mutate(
      TrimTag = case_when(
        str_detect(Outcome, "Trim\\s*1") ~ "Trim1",
        str_detect(Outcome, "Trim\\s*2") ~ "Trim2",
        str_detect(Outcome, "Trim\\s*3") ~ "Trim3",
        str_detect(Outcome, "Pp\\s*6") ~ "Pp6wks",
        TRUE ~ "All"
      ),
      Outcome_clean = sub("\\s*\\(\\s*\\.50\\s*(,\\s*Trim\\s*[123])?\\s*\\)\\s*$", "", Outcome)
    ) %>%
    filter(!(Outcome_clean %in% exclude_outcomes))
  
  # If any Trim tags exist, filter by the requested one; else keep all and tag as that trimester
  if (any(df0$TrimTag != "All")) df0 <- df0 %>% filter(TrimTag == trimester) else df0 <- df0 %>% mutate(TrimTag = trimester)
  if (nrow(df0) == 0) {
    return(tibble(Outcome = character(),
                  `Normal (lowest risk)` = character(),
                  `Highest Risk - Low Threshold` = character(),
                  `Highest Risk - High Threshold` = character()))
  }
  
  # 2) Round for ranking, compute per-outcome min/max and unique count
  df <- df0 %>%
    mutate(risk_key = ifelse(is.na(risk), NA_real_, round(risk, risk_digits))) %>%
    group_by(Outcome_clean) %>%
    mutate(
      nuniq    = n_distinct(risk_key, na.rm = TRUE),
      rmin_key = suppressWarnings(min(risk_key, na.rm = TRUE)),
      rmax_key = suppressWarnings(max(risk_key, na.rm = TRUE))
    ) %>%
    ungroup()
  
  # 3) Optional drops
  if (isTRUE(drop_mono_risk)) {
    df <- df %>% group_by(Outcome_clean) %>% filter(nuniq > 1) %>% ungroup()
  }
  if (isTRUE(drop_single_threshold) && "group" %in% names(df0)) {
    df <- df %>% group_by(Outcome_clean) %>% filter(n_distinct(group[!is.na(group)]) > 1) %>% ungroup()
  }
  if (nrow(df) == 0) {
    return(tibble(Outcome = character(),
                  `Normal (lowest risk)` = character(),
                  `Highest Risk - Low Threshold` = character(),
                  `Highest Risk - High Threshold` = character()))
  }
  
  # 4) Identify the EARLIEST contiguous min segment safely (no scalar if() in mutate)
  df_min <- df %>%
    arrange(Outcome_clean, hb) %>%
    group_by(Outcome_clean) %>%
    mutate(
      is_min = !is.na(risk_key) & (risk_key == rmin_key),
      # start of a "min" run: when entering is_min or hb jump != 0.5
      min_run_start = is_min &
        (is.na(lag(is_min)) | !lag(is_min) | abs(hb - lag(hb) - 0.5) > tol),
      # assign run IDs only to min rows; NA for non-min rows
      min_seg_id = ifelse(is_min, cumsum(replace_na(min_run_start, FALSE)), NA_integer_)
    ) %>%
    ungroup()
  
  # Earliest min segment per outcome (lowest Hb start)
  first_min_seg <- df_min %>%
    filter(is_min) %>%
    group_by(Outcome_clean, min_seg_id) %>%
    summarise(seg_start_hb = min(hb, na.rm = TRUE), .groups = "drop") %>%
    group_by(Outcome_clean) %>%
    arrange(seg_start_hb, min_seg_id, .by_group = TRUE) %>%
    slice(1) %>%
    ungroup() %>%
    select(Outcome_clean, first_min_seg_id = min_seg_id)
  
  df <- df_min %>%
    left_join(first_min_seg, by = "Outcome_clean") %>%
    mutate(is_normal = is_min & !is.na(first_min_seg_id) & (min_seg_id == first_min_seg_id))
  
  # 5) Tail logic (near-max in low/high Hb), never override Normal
  df <- df %>%
    mutate(
      in_low_tail   = hb <  low_tail_cut,
      in_high_tail  = hb >= high_tail_cut,
      is_max        = !is.na(risk_key) & (risk_key == rmax_key),
      
      # For >=3 risks use near-max rule; for exactly 2 risks, the max itself qualifies
      near_max_ok = dplyr::case_when(
        nuniq >= 3 ~ (!is.na(risk_key) & (risk_key >= near_max_frac * rmax_key)),
        nuniq == 2 ~ is_max,
        TRUE       ~ FALSE
      ),
      
      is_high_low   = !is_normal & in_low_tail  & near_max_ok,
      is_high_high  = !is_normal & in_high_tail & near_max_ok,
      
      class = dplyr::case_when(
        is_normal    ~ "Normal (lowest risk)",
        is_high_low  ~ "Highest Risk - Low Threshold",
        is_high_high ~ "Highest Risk - High Threshold",
        TRUE         ~ NA_character_
      )
    )
  
  # 6) Build ranges per class and pivot wide
  out <- df %>%
    filter(!is.na(class)) %>%
    group_by(Outcome_clean, class) %>%
    summarise(hb_range = .collapse_ranges(hb), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = class, values_from = hb_range) %>%
    arrange(Outcome_clean) %>%
    rename(Outcome = Outcome_clean)
  
  # Ensure all columns exist
  need_cols <- c("Normal (lowest risk)",
                 "Highest Risk - Low Threshold",
                 "Highest Risk - High Threshold")
  for (nm in need_cols) if (!nm %in% names(out)) out[[nm]] <- NA_character_
  
  out %>% select(Outcome, all_of(need_cols))
}


table_all  <- heatmap_logic_to_table(hb_grid,    trimester = "All")
table_t1   <- heatmap_logic_to_table(hb_grid_t1, trimester = "Trim1")
table_t2   <- heatmap_logic_to_table(hb_grid_t2, trimester = "Trim2")
table_t3   <- heatmap_logic_to_table(hb_grid_t3, trimester = "Trim3")
table_p6   <- heatmap_logic_to_table(hb_grid_p6, trimester = "Pp6wks")

suppressPackageStartupMessages({
  library(dplyr)
  library(gt)
})

# --- order buckets exactly as listed ---
maternal_outcomes <- c("MAT COMPO",
                       "PRECLAMP","PPROM","PPH","PPA PNC6","PPA PNC26",
                       "FATIGUE PNC","FATIGUE ANC","DPR PNC","DPR ANC"
)
infant_outcomes <- c(
  "STILLBIRTH 28WKS","STILLBIRTH 20WKS", "NEONATAL MORTALITY","SGA3","SGA10","PSBI",
  "PRETERM 37WKS","PRETERM 34WKS","LBW<2500g","LBW<1500g",
  "HYPERBILIRUBINEMIA","COMPO"
)

# display order for labels (matches your code order)
display_order <- c(outcome_label_map[maternal_outcomes],
                   outcome_label_map[infant_outcomes])

# ---- helper: make a complete, ordered, gray-missing summary table ----
make_summary_gt <- function(tbl, title, group_rows = TRUE) {
  # canonical column names from heatmap_logic_to_table()
  need_cols <- c(
    "Outcome",
    "Normal (lowest risk)",
    "Highest Risk-Low Threshold",
    "Highest Risk-High Threshold"
  )
  
  # accept either spacing and normalize
  tbl <- tbl %>%
    dplyr::rename(
      `Highest Risk-Low Threshold`  = dplyr::any_of("Highest Risk - Low Threshold"),
      `Highest Risk-High Threshold` = dplyr::any_of("Highest Risk - High Threshold")
    )
  
  # ensure columns exist
  for (nm in setdiff(need_cols, names(tbl))) tbl[[nm]] <- NA_character_
  
  listed   <- c(maternal_outcomes, infant_outcomes)
  extras   <- setdiff(tbl$Outcome, listed)
  full_ord <- c(listed, sort(extras))
  
  complete <- tibble::tibble(Code = full_ord) %>%
    dplyr::left_join(tbl, by = c("Code" = "Outcome")) %>%
    dplyr::mutate(
      Group = dplyr::case_when(
        Code %in% maternal_outcomes ~ "Maternal Outcomes",
        Code %in% infant_outcomes   ~ "Infant Outcomes",
        TRUE                        ~ "Other"
      ),
      # coerce to character for safe lookup
      Display = dplyr::coalesce(outcome_label_map[as.character(Code)], Code),
      # vectorized missing-row flag across the three range cols
      missing_row = dplyr::if_all(
        dplyr::all_of(need_cols[-1]),
        ~ dplyr::coalesce(as.character(.x), "") == ""
      ),
      Code    = factor(Code, levels = full_ord),
      Display = factor(
        Display,
        levels = c(display_order, as.character(setdiff(Display, display_order)))
      )
    ) %>%
    dplyr::arrange(Code)
  
  # ---- inline theme (stable across gt versions) ----
  local_theme <- function(gt_tbl) {
    gt_tbl %>%
      gt::gt() %>%
      gt::tab_options(
        table.font.names = c("Inter","Arial","Helvetica","sans-serif"),
        table.font.size = gt::px(13),
        data_row.padding = gt::px(6),
        row.striping.include_table_body = TRUE,
        row.striping.background_color = "#F9FAFB",
        column_labels.background.color = "#9edee0",
        table.border.top.style = "none",
        table.border.bottom.style = "solid",
        table.border.bottom.width = gt::px(1),
        table.border.bottom.color = "#E5E7EB",
        table.width = gt::px(650)   # width belongs in tab_options, not tab_style
      ) %>%
      # Title styling (both title & subtitle areas)
      gt::tab_style(
        style = list(
          gt::cell_fill(color = "#016795"),
          gt::cell_text(weight = "bold", color = "white", v_align = "middle")
        ),
        locations = gt::cells_title(groups = c("title","subtitle"))
      ) %>%
      # Stub head, spanners, and column labels
      gt::tab_style(
        style = list(
          gt::cell_fill(color = "#9edee0"),
          gt::cell_text(weight = "bold", color = "black", v_align = "middle")
        ),
        locations = list(
          gt::cells_stubhead(),
          gt::cells_column_spanners(),
          gt::cells_column_labels(columns = gt::everything())
        )
      )
  }
  
  g <- complete %>%
    dplyr::select(
      Code, Display,
      `Normal (lowest risk)`,
      `Highest Risk-Low Threshold`,
      `Highest Risk-High Threshold`,
      Group, missing_row
    ) %>%
    local_theme() %>%
    # shorter, intuitive headers
    gt::cols_label(
      Display = "Outcome",
      `Normal (lowest risk)`        = "Normal",
      `Highest Risk-Low Threshold`  = "High risk (low Hb)",
      `Highest Risk-High Threshold` = "High risk (high Hb)"
    ) %>%
    gt::tab_header(title = gt::md(title)) %>%
    # format markdown in all cells (lets you pass md() in values)
    gt::fmt_markdown(columns = gt::everything()) %>%
    # hide NA text everywhere
    gt::sub_missing(columns = gt::everything(), rows = gt::everything(), missing_text = "") %>%
    # highlight the first visible column (Outcome/Display)
    gt::tab_style(
      style = list(
        gt::cell_fill(color = "#f2fafa"),
        gt::cell_text(weight = "bold")
      ),
      locations = gt::cells_body(columns = "Display")
    ) %>%
    # subtle tints on risk columns (only where non-empty)
    gt::tab_style(
      style = list(gt::cell_text(color = "#065F46")),
      locations = gt::cells_body(
        columns = "Normal (lowest risk)",
        rows = !is.na(`Normal (lowest risk)`) & `Normal (lowest risk)` != ""
      )
    ) %>%
    gt::tab_style(
      style = list(gt::cell_text(color = "#991B1B")),
      locations = gt::cells_body(
        columns = "Highest Risk-Low Threshold",
        rows = !is.na(`Highest Risk-Low Threshold`) & `Highest Risk-Low Threshold` != ""
      )
    ) %>%
    gt::tab_style(
      style = list(gt::cell_text(color = "#9A3412")),
      locations = gt::cells_body(
        columns = "Highest Risk-High Threshold",
        rows = !is.na(`Highest Risk-High Threshold`) & `Highest Risk-High Threshold` != ""
      )
    )
  
  # row groups + bold group labels (apply after groups exist)
  if (isTRUE(group_rows)) {
    g <- g %>%
      gt::tab_row_group(group = "Maternal Outcomes", rows = Code %in% maternal_outcomes) %>%
      gt::tab_row_group(group = "Infant Outcomes",   rows = Code %in% infant_outcomes) %>%
      { if (any(complete$Group == "Other"))
        gt::tab_row_group(., group = "Other", rows = !(Code %in% c(maternal_outcomes, infant_outcomes)))
        else . } %>%
      gt::tab_style(
        style = gt::cell_text(weight = "bold"),
        locations = gt::cells_row_groups()   # all row-group labels
      ) %>% opt_align_table_header(align = "left")
  }
  
  # gray-out rows with no ranges (overrides any tints)
  g %>%
    gt::tab_style(
      style = list(
        gt::cell_fill(color = "#E6E6E6"),
        gt::cell_text(color = "#7D7D7D")
      ),
      locations = gt::cells_body(rows = missing_row, columns = gt::everything())
    ) %>%
    gt::cols_hide(columns = c("missing_row","Group","Code")) %>%
    gt::cols_align(align = "left", columns = "Display")
}

# ---- Build the four summary tables (All, Trim1, Trim2, Trim3) ----
gt_sum_all <- make_summary_gt(table_all, "**Hb Threshold Ranges - All Trimesters**")
gt_sum_t1  <- make_summary_gt(table_t1,  "**Hb Threshold Ranges — Trimester 1**")
gt_sum_t2  <- make_summary_gt(table_t2,  "**Hb Threshold Ranges — Trimester 2**")
gt_sum_t3  <- make_summary_gt(table_t3,  "**Hb Threshold Ranges — Trimester 3**")

gtsave(gt_sum_all, "gt_sum_all.png")
gtsave(gt_sum_t1,  "gt_sum_t1.png")
gtsave(gt_sum_t2,  "gt_sum_t2.png")
gtsave(gt_sum_t3,  "gt_sum_t3.png")
```


### All Trimesters Summary Table
```{r, out.width='90%', fig.align = "left", echo=FALSE, include=TRUE}

knitr::include_graphics("gt_sum_all.png")

```

\newpage
### Trimester 1 Summary Table 
```{r, out.width='90%',fig.align = "left", echo=FALSE, include=TRUE}
knitr::include_graphics("gt_sum_t1.png")
```

\newpage
### Trimester 2 Summary Table 
```{r, out.width='90%', fig.align = "left", echo=FALSE, include=TRUE}
knitr::include_graphics("gt_sum_t2.png")
```

\newpage
### Trimester 3 Summary Table
```{r, out.width='90%', fig.align = "left", echo=FALSE, include=TRUE}
knitr::include_graphics("gt_sum_t3.png")
```