
out <- out_sga3_50  # Your isotonic model output
out <- out_sga3_50[, colSums(!is.na(out_sga3_50)) > 0]

# Get variable names
n_vars     <- grep("^N[0-9]+$", names(out), value = TRUE)
risk_vars  <- grep("^Risk\\.Mscore[0-9]+$", names(out), value = TRUE)
thres_vars <- grep("^Thres[0-9]+$", names(out), value = TRUE)

# Extract group numbers safely
n_groups     <- as.integer(gsub("N", "", n_vars))
risk_groups  <- as.integer(gsub("Risk\\.Mscore", "", risk_vars))
thres_groups <- as.integer(gsub("Thres", "", thres_vars))

# Keep only group numbers common to all
common_groups <- intersect(intersect(n_groups, risk_groups), thres_groups)

# Build long-format dataframeS
long_df <- lapply(common_groups, function(i) {
  data.frame(
    Group = i,
    N = out[[paste0("N", i)]],
    Risk = round(out[[paste0("Risk.Mscore", i)]], 3),
    Threshold = out[[paste0("Thres", i)]]
  )
}) |> bind_rows()

# View as flextable
model_label <- "SGA 3rd (0.50 change model)"  # change as needed

flextable(long_df) |>
  set_caption(paste("Isotonic Threshold –", model_label)) |>
  autofit() |>
  theme_booktabs()


hist(data$hb)
