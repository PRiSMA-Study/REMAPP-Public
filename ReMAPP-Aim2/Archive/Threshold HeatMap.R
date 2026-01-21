
suppressPackageStartupMessages({
  library(dplyr); library(stringr); library(purrr); library(tidyr)
})

#All General Codes ----
# 1. Parse the exact label formats produced by outdata() ----
parse_iso_label <- function(label, hb_min = 5, hb_max = 18, eps = 1e-6) {
  if (length(label) == 0 || is.na(label) || !nzchar(label)) {
    return(tibble(lower = NA_real_, upper = NA_real_))
  }
  
  # Convert to string and normalize
  lab <- as.character(label) |>
    stringi::stri_trim_both() |>
    stringi::stri_replace_all_regex("[\\p{Z}\\s]+", " ") |>  # Normalize spaces
    stringi::stri_replace_all_regex("[–—]", "-")              # Normalize dashes
  
  # Define patterns (using explicit Unicode in regex)
  patterns <- list(
    # Pattern 1: "≤ b" or "<= b"
    list(
      regex = "^[\\u2264<=]\\s*(-?\\d+(?:\\.\\d+)?)$",
      handler = function(m) tibble(lower = hb_min, upper = as.numeric(m[,2]))
    ),
    # Pattern 2: "> a"
    list(
      regex = "^>\\s*(-?\\d+(?:\\.\\d+)?)$",
      handler = function(m) tibble(lower = as.numeric(m[,2]) + eps, upper = hb_max)
    ),
    # Pattern 3: "≥ a" or ">= a"
    list(
      regex = "^[\\u2265>=]\\s*(-?\\d+(?:\\.\\d+)?)$",
      handler = function(m) tibble(lower = as.numeric(m[,2]), upper = hb_max)
    ),
    # Pattern 4: "> a - ≤ b"
    list(
      regex = "^>\\s*(-?\\d+(?:\\.\\d+)?)\\s*-\\s*[\\u2264<=]\\s*(-?\\d+(?:\\.\\d+)?)$",
      handler = function(m) tibble(lower = as.numeric(m[,2]) + eps, upper = as.numeric(m[,3]))
    ),
    # Pattern 5: "≥ a - ≤ b"
    list(
      regex = "^[\\u2265>=]\\s*(-?\\d+(?:\\.\\d+)?)\\s*-\\s*[\\u2264<=]\\s*(-?\\d+(?:\\.\\d+)?)$",
      handler = function(m) tibble(lower = as.numeric(m[,2]), upper = as.numeric(m[,3]))
    ),
    # Pattern 6: "All values"
    list(
      regex = "^all\\s*values?$",
      handler = function(m) tibble(lower = hb_min, upper = hb_max),
      ignore_case = TRUE
    )
  )
  
  # Try each pattern
  for (p in patterns) {
    opts <- if (isTRUE(p$ignore_case)) stringi::stri_opts_regex(case_insensitive = TRUE) else NULL
    if (stringi::stri_detect_regex(lab, p$regex, opts_regex = opts)) {
      m <- stringi::stri_match_first_regex(lab, p$regex, opts_regex = opts)
      return(p$handler(m))
    }
  }
  
  warning("Unrecognized label format: '", lab, "'")
  tibble(lower = NA_real_, upper = NA_real_)
}

# 2. Tidy one outdata() row -> interval rows (no guessing, no summarizing) ----
tidy_outdata_iso <- function(out_df, hb_min = 5, hb_max = 18, eps = 1e-6) {
  # First find threshold columns
  thres_cols <- grep("^Thres\\d+$", names(out_df), value = TRUE)
  if (!length(thres_cols)) stop("No Thres<i> columns found.")
  
  # Order them numerically
  thres_cols <- thres_cols[order(as.integer(str_extract(thres_cols, "\\d+")))]
  
  outcome <- out_df$Outcome[1]
  
  map_dfr(thres_cols, function(col) {
    i <- as.integer(str_extract(col, "\\d+"))
    parsed <- parse_iso_label(out_df[[col]], hb_min, hb_max, eps)
    parsed |>
      mutate(
        Outcome = outcome,
        group   = i,
        risk    = as.numeric(out_df[[paste0("Risk/Mscore", i)]]),
        n       = as.numeric(out_df[[paste0("N", i)]])
      ) |>
      relocate(Outcome, group, lower, upper, risk, n)
  })
}


# 3.Interval creation from environment ----
collect_intervals_from_env <- function(
    which = c("All","Trim1","Trim2","Trim3"),
    where = .GlobalEnv,
    hb_min = 5, hb_max = 18, eps = 1e-6,
    pattern = NULL,                 # optional: override regex manually
    assign_intermediates = TRUE
) {
  which <- match.arg(which)
  
  # Build an anchored, case-insensitive pattern if not provided
  if (is.null(pattern)) {
    base <- "(?i)^out.*_50"
    suffix <- switch(which,
                     "All"   = "$",
                     "Trim1" = "_?trim1$",
                     "Trim2" = "_?trim2$",
                     "Trim3" = "_?trim3$"
    )
    pattern <- paste0(base, suffix)
  }
  
  # Find matching objects
  objs  <- ls(envir = where, all.names = TRUE)
  picks <- grep(pattern, objs, value = TRUE, perl = TRUE)
  if (!length(picks)) stop("No objects matched pattern: ", pattern)
  
  outs <- mget(picks, envir = where, inherits = FALSE)
  
  # Keep only data.frames that look like outdata() results
  dfs  <- outs[vapply(outs, is.data.frame, logical(1))]
  dfs  <- dfs[vapply(dfs, function(x) any(grepl("^Thres\\d+$", names(x))), logical(1))]
  if (!length(dfs)) stop("Matched objects, but none look like outdata() frames for: ", which)
  
  # Tidy each and optionally expose intermediates
  intermediate_dfs <- lapply(dfs, tidy_outdata_iso, hb_min = hb_min, hb_max = hb_max, eps = eps)
  if (isTRUE(assign_intermediates)) {
    list2env(setNames(intermediate_dfs, paste0("intermediate_", names(intermediate_dfs))), envir = where)
  }
  
  dplyr::bind_rows(intermediate_dfs)
}

#4. Expand to grid ----
expand_to_grid <- function(intervals, hb_seq = seq(5, 18, by = 0.5)) {
  by_outcome <- split(intervals, intervals$Outcome)
  bind_rows(lapply(names(by_outcome), function(outc) {
    intv <- arrange(by_outcome[[outc]], lower, upper, group)
    idx <- vapply(hb_seq, function(hb) {
      hit <- which(!is.na(intv$lower) & !is.na(intv$upper) &
                     hb >= intv$lower & hb <= intv$upper)
      if (length(hit)) hit[1] else NA_integer_
    }, integer(1))
    tibble(Outcome = outc, hb = hb_seq, risk = ifelse(is.na(idx), NA_real_, intv$risk[idx]))
  }))
}

# 5. Risk Plot Function ----
# --- groups (COMPO in Infant) ---
maternal_outcomes <- c(
  "PRECLAMP","PPROM","PPH","PPA PNC6","PPA PNC26",
  "FATIGUE PNC","FATIGUE ANC","DPR PNC","DPR ANC"
)
infant_outcomes <- c(
  "STILLBIRTH 28WKS","STILLBIRTH 20WKS","SGA3","SGA10","PSBI",
  "PRETERM 37WKS","PRETERM 34WKS","LBW<2500g","LBW<1500g",
  "HYPERBILIRUBINEMIA","COMPO"
)

plot_risk_heatmap <- function(grid_df,
                              trimester = NULL,
                              title = NULL,
                              show_values = FALSE,
                              risk_digits = 2,
                              group_order = c("Reference","Maternal","Infant"),
                              drop_mono_risk = TRUE,
                              drop_single_threshold = TRUE,
                              non_mono_last = TRUE,
                              add_reference_rows = TRUE,
                              reference_group = "Reference") {
  
  # grid_df = hb_grid
  # trimester = "All"
  # title = NULL
  # show_values = FALSE
  # risk_digits = 2
  # group_order = c("Reference","Maternal","Infant")
  # drop_mono_risk = TRUE
  # drop_single_threshold = TRUE
  # non_mono_last = TRUE
  # add_reference_rows = TRUE
  # reference_group = "Reference"
  
  
  stopifnot(all(c("Outcome","hb","risk") %in% names(grid_df)))
  suppressPackageStartupMessages({
    library(dplyr); library(ggplot2); library(forcats); library(stringr); library(tidyr); library(tibble)
  })
  
  
  mono_grey <- "#bfbfbf"
  
  # --- parse trimester + clean base outcome name ---
  df0 <- grid_df %>%
    mutate(
      Trimester = str_match(Outcome, "Trim\\s*(\\d)")[,2],
      Trimester = dplyr::recode(Trimester, `1`="Trim1", `2`="Trim2", `3`="Trim3", .default = "All"),
      Outcome_base = Outcome %>%
        str_replace("\\s*\\(\\s*\\.50\\s*(,\\s*Trim\\s*[123])?\\s*\\)\\s*$", "") %>%
        str_squish()
    )
  
  # --- build the dataset for the requested trimester ---
  df_trim <- df0 
  df <- df_trim %>%
    mutate(Outcome_clean = Outcome_base,
           Group = case_when(
             Outcome_clean %in% maternal_outcomes ~ "Maternal",
             Outcome_clean %in% infant_outcomes   ~ "Infant",
             TRUE ~ NA
           )) %>%
    select(Outcome_clean, Group, Trimester, hb, risk, dplyr::any_of("group"))
  # threshold-count for this trimester only
  thresh_tbl <- df_trim %>%
    group_by(Outcome_base) %>%
    summarise(n_thresh = if ("group" %in% names(df_trim)) n_distinct(group[!is.na(group)]) else NA_integer_,
              .groups = "drop")
  df <- df %>% left_join(thresh_tbl, by = c("Outcome_clean" = "Outcome_base"))
  
  # if nothing in this trimester (e.g., no Trim2 data), stop gracefully unless reference rows requested
  if (nrow(df) == 0 && !add_reference_rows) stop(paste("No data for", trimester))
  
  # --- ranking + display formatting ---
  # if (nrow(df) > 0) {
  #   df <- df %>%
  #     mutate(
  #       risk_key = ifelse(is.na(risk), NA_real_, round(risk, risk_digits)),
  #       risk_fmt = ifelse(is.na(risk), "", sprintf(paste0("%.", risk_digits, "f"), risk))
  #     )
  # 
  #   # per-outcome rank classes
  #   stats <- df %>%
  #     group_by(Outcome_clean) %>%
  #     summarise(
  #       all_na   = all(is.na(risk_key)),
  #       nuniq    = n_distinct(risk_key, na.rm = TRUE),
  #       rmin_key = if (all_na) NA_real_ else min(risk_key, na.rm = TRUE),
  #       rmax_key = if (all_na) NA_real_ else max(risk_key, na.rm = TRUE),
  #       .groups = "drop"
  #     )
  #   df <- df %>% left_join(stats, by = "Outcome_clean")
  # 
  #   # optional drops
  #   df <- df %>%
  #     group_by(Outcome_clean) %>%
  #     filter(
  #       !all_na,
  #       !(isTRUE(drop_mono_risk)        & nuniq <= 1),
  #       !(isTRUE(drop_single_threshold) & !is.na(n_thresh) & n_thresh <= 1)
  #     ) %>%
  #     ungroup()
  # 
  #   # label category
  #   if (nrow(df) > 0) {
  #     df <- df %>%
  #       mutate(
  #         is_max = !is.na(risk_key) & risk_key >= rmax_key,
  #         is_min = !is.na(risk_key) & risk_key <= rmin_key,
  #         fill_cat = case_when(
  #           is_max ~ "High",
  #           is_min ~ "Low",
  #           TRUE   ~ "Medium"
  #         ),
  #         fill_cat = factor(fill_cat, levels = c("High","Medium","Low"))
  #       )
  #   }
  # }
  if (nrow(df) > 0) {
    df <- df %>%
      mutate(
        risk_key = ifelse(is.na(risk), NA_real_, round(risk, risk_digits)),
        risk_fmt = ifelse(is.na(risk), "", sprintf(paste0("%.", risk_digits, "f"), risk))
      )
    
    # per-outcome stats
    stats <- df %>%
      group_by(Outcome_clean) %>%
      summarise(
        all_na   = all(is.na(risk_key)),
        nuniq    = n_distinct(risk_key, na.rm = TRUE),
        rmin_key = if (all_na) NA_real_ else min(risk_key, na.rm = TRUE),
        rmax_key = if (all_na) NA_real_ else max(risk_key, na.rm = TRUE),
        .groups = "drop"
      )
    df <- df %>% left_join(stats, by = "Outcome_clean")
    
    # optional drops
    df <- df %>%
      group_by(Outcome_clean) %>%
      filter(
        !all_na,
        !(isTRUE(drop_mono_risk)        & nuniq <= 1),
        !(isTRUE(drop_single_threshold) & !is.na(n_thresh) & n_thresh <= 1)
      ) %>%
      ungroup()
    
    if (nrow(df) > 0) {
      # Tail boundaries and near-max threshold
      low_tail_cut   <- 8     # hb < 8 considered "low tail"
      high_tail_cut  <- 13    # hb >= 13 considered "high tail"
      near_max_frac  <- 0.90  # ≥ 90% of per-outcome max
      
      df <- df %>%
        group_by(Outcome_clean) %>%
        mutate(
          is_global_min = !is.na(risk_key) & risk_key == rmin_key,
          is_global_max = !is.na(risk_key) & risk_key == rmax_key,
          in_low_tail   = hb <  low_tail_cut,
          in_high_tail  = hb >= high_tail_cut,
          
          # Near-max in tails only if we have >=3 unique risks
          is_near_max_tail = nuniq >= 3 &
            !is.na(risk_key) &
            risk_key >= near_max_frac * rmax_key &
            (in_low_tail | in_high_tail)
        ) %>%
        ungroup() %>%
        mutate(
          # IMPORTANT: global min stays green even if it meets the near-max test
          fill_cat = case_when(
            is_global_min                          ~ "Low",
            is_global_max | is_near_max_tail       ~ "High",
            TRUE                                   ~ "Medium"
          ),
          fill_cat = factor(fill_cat, levels = c("High","Medium","Low"))
        )
    }
  }
  
  
  # --- reference rows (WHO, HEALTHY) for this trimester plot ---
  if (isTRUE(add_reference_rows)) {
    hb_vals <- if (nrow(grid_df)) sort(unique(grid_df$hb)) else seq(5, 18, by = 0.5)
    
    # Which set of cuts to use
    trim_key <- if (trimester %in% c("Trim1","Trim2","Trim3")) trimester else "All"
    
    #EDIT THESE DUMMY CUTS 
    # WHO: <= low -> High; (low, mid] -> Medium; > mid -> Low
    who_cuts <- list(
      All   = c(low = 7.0, mid = 11.0),
      Trim1 = c(low = 7.0, mid = 11.0),
      Trim2 = c(low = 7.0, mid = 10.5),
      Trim3 = c(low = 7.0, mid = 11.0)
    )
    
    # Healthy Cohort:
    # < h1 -> High; [h1, h2) -> Medium; [h2, h3] -> Low; > h3 -> High
    healthy_cuts <- list(
      All   = c(h1 = 8.23, h2 = 9.0, h3 = 13.0),
      Trim1 = c(h1 = 8.95, h2 = 9.75, h3 = 13.5),
      Trim2 = c(h1 = 8.37, h2 = 9.15, h3 = 13.0),
      Trim3 = c(h1 = 8.11, h2 = 9.00, h3 = 13)
    )
    
    # Helpers to turn cuts into "High/Medium/Low" fills
    fill_who <- function(x, cuts) dplyr::case_when(
      x <= cuts["low"] ~ "High",
      x <= cuts["mid"] ~ "Medium",
      TRUE             ~ "Low"
    )
    fill_healthy <- function(x, cuts) dplyr::case_when(
      x <  cuts["h1"] ~ "High",
      x <  cuts["h2"] ~ "Medium",
      x <= cuts["h3"] ~ "Low",
      TRUE            ~ "High"
    )
    
    who_fill     <- fill_who(hb_vals,     who_cuts[[trim_key]])
    healthy_fill <- fill_healthy(hb_vals, healthy_cuts[[trim_key]])
    
    ref_df <- dplyr::bind_rows(
      tibble::tibble(
        Outcome_clean = "WHO Guideline",
        Group         = reference_group,
        Trimester     = trimester,
        hb            = hb_vals,
        risk_key      = NA_real_,
        risk_fmt      = "",
        fill_cat      = factor(who_fill, levels = c("High","Medium","Low"))
      ),
      tibble::tibble(
        Outcome_clean = "Healthy Cohort",
        Group         = reference_group,
        Trimester     = trimester,
        hb            = hb_vals,
        risk_key      = NA_real_,
        risk_fmt      = "",
        fill_cat      = factor(healthy_fill, levels = c("High","Medium","Low"))
      )
    )
    
    # combine (allow empty df on the left)
    if (nrow(df) == 0) {
      df <- ref_df
    } else {
      keep <- intersect(names(df), names(ref_df))
      df   <- dplyr::bind_rows(df[, keep], ref_df[, keep])
    }
    
  }
  
  if (is.null(title)) {
    title <- switch(trimester,
                    "Trim1" = "Heatmap of Hb Threshold Risk Rank — Trimester 1",
                    "Trim2" = "Heatmap of Hb Threshold Risk Rank — Trimester 2",
                    "Trim3" = "Heatmap of Hb Threshold Risk Rank — Trimester 3",
                    "All"   = "Heatmap of Hb Threshold Risk Rank — All Trimesters")
  }
  
  # facet order & row order
  df$Group <- fct_relevel(df$Group, group_order)
  ref_order <- c("WHO Guideline","Healthy Cohort")
  row_order <- df %>%
    distinct(Group, Outcome_clean) %>%
    arrange(
      factor(Group, levels = group_order),
      dplyr::case_when(
        Group == reference_group & Outcome_clean %in% ref_order ~ match(Outcome_clean, ref_order),
        TRUE ~ Inf
      ),
      Outcome_clean
    )
  df <- df %>% mutate(Outcome_label = factor(Outcome_clean, levels = row_order$Outcome_clean))
  
  # plot
  risk_colors <- c("High"="#e74c3c","Medium"="#f1c40f","Low"="#2ecc71")
  
  g <- ggplot(df, aes(x = hb, y = Outcome_label, fill = fill_cat)) +
    geom_tile(height = 1, width = 0.5, colour = mono_grey, linewidth = 0.25) +
    scale_fill_manual(values = risk_colors,
                      name = "Threshold Risk Rank (per outcome)",
                      na.value = "grey90") +
    { if (isTRUE(show_values)) geom_text(aes(label = risk_fmt), size = 3) else NULL } +
    scale_x_continuous(breaks = seq(5, 18, by = 1), expand = c(0,0)) +
    facet_grid(Group ~ ., scales = "free_y", space = "free_y") +
    labs(x = "Hemoglobin (g/dL)", y = NULL, title = title) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      legend.position = "right",
      strip.text.y = element_text(face = "bold"),
      axis.text.y = element_text(size = 8, hjust = 1, margin = margin(r = 6)),
      axis.text.x = element_text(size = 8),
      axis.title  = element_text(size = 9),
      plot.title  = element_text(size = 10, hjust = 0),
      legend.text = element_text(size = 7),
      legend.title= element_text(size = 8),
      legend.key.width  = grid::unit(1, "cm"),
      legend.key.height = grid::unit(0.3, "cm")
    )
  g
}


#All trimesters
# List all .rda files starting with 'out_' and ending with '_50.rda' 
files <- list.files("iso_results", pattern = "^out_.*_50\\.rda$", full.names = TRUE)

# Load each file and assign it to its original name (without .rda) ----
for (file in files) {
  obj_name <- load(file)  # Loads the object and returns its stored name
  file_name <- tools::file_path_sans_ext(basename(file))  # e.g., "out_anc_dpr_50"
  assign(file_name, get(obj_name), envir = .GlobalEnv)
}


##Trimester 1 ----

# List all .rda files starting with 'out_' and ending with '_50_trim1.rda' 
files <- list.files("iso_results", pattern = "^out_.*_50_trim1\\.rda$", full.names = TRUE)

# Load each file and assign it to its original name (without .rda) 
for (file in files) {
  obj_name <- load(file)  # Loads the object and returns its stored name
  file_name <- tools::file_path_sans_ext(basename(file))  # e.g., "out_anc_dpr_50"
  assign(file_name, get(obj_name), envir = .GlobalEnv)
}

##Trimester 2 ----

# List all .rda files starting with 'out_' and ending with '_50_trim2.rda' 
files <- list.files("iso_results", pattern = "^out_.*_50_trim2\\.rda$", full.names = TRUE)

# Load each file and assign it to its original name (without .rda) 
for (file in files) {
  obj_name <- load(file)  # Loads the object and returns its stored name
  file_name <- tools::file_path_sans_ext(basename(file))  # e.g., "out_anc_dpr_50"
  assign(file_name, get(obj_name), envir = .GlobalEnv)
}


##Trimester 3 ----
# List all .rda files starting with 'out_' and ending with '_50_trim3.rda' 
files <- list.files("iso_results", pattern = "^out_.*_50_trim3\\.rda$", full.names = TRUE)

# Load each file and assign it to its original name (without .rda) 
for (file in files) {
  obj_name <- load(file)  # Loads the object and returns its stored name
  file_name <- tools::file_path_sans_ext(basename(file))  # e.g., "out_anc_dpr_50"
  assign(file_name, get(obj_name), envir = .GlobalEnv)
}


#Create all Interval Datasets ----
intervals        <- collect_intervals_from_env("All")
intervals_trim1  <- collect_intervals_from_env("Trim1")
intervals_trim2  <- collect_intervals_from_env("Trim2")
intervals_trim3  <- collect_intervals_from_env("Trim3")

#Make Hb Grid Combined Datasets for All Trimesters ----
hb_grid   <- expand_to_grid(intervals,        hb_seq = seq(5, 18, by = 0.5))
hb_grid_t1 <- expand_to_grid(intervals_trim1,  hb_seq = seq(5, 18, by = 0.5))
hb_grid_t2 <- expand_to_grid(intervals_trim2,  hb_seq = seq(5, 18, by = 0.5))
hb_grid_t3 <- expand_to_grid(intervals_trim3,  hb_seq = seq(5, 18, by = 0.5))
 
#Plot all Plots By Trimester ----
p_all <- plot_risk_heatmap(hb_grid,   trimester = "All")
p_t1  <- plot_risk_heatmap(hb_grid_t1, trimester = "Trim1")
p_t2  <- plot_risk_heatmap(hb_grid_t2, trimester = "Trim2")
p_t3  <- plot_risk_heatmap(hb_grid_t3, trimester = "Trim3")

#Print all plots ----
print(p_all)
print(p_t1)
print(p_t2)
print(p_t3)
