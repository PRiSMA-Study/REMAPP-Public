# outdata <- function(data, outcome) {
#   # Load required package
#   if (!requireNamespace("dplyr", quietly = TRUE)) {
#     stop("Please install the 'dplyr' package first")
#   }
#   
#   # 1. Safely process group statistics
#   df <- tryCatch({
#     data.frame(
#       group = data$groups,
#       estimate = data$estimates
#     ) %>%
#       dplyr::group_by(group) %>%
#       dplyr::summarise(
#         group_n = dplyr::n(),
#         risk_mscore = mean(estimate, na.rm = TRUE),
#         .groups = "drop"
#       ) %>%
#       dplyr::arrange(group)
#   }, error = function(e) {
#     stop("Error processing group statistics: ", e$message)
#   })
#   
#   # 2. Process breakpoints safely
#   breaks <- tryCatch({
#     sort(unique(data$brkPoints))
#   }, error = function(e) {
#     warning("Error processing breakpoints: ", e$message)
#     numeric(0)
#   })
#   
#   # 3. Create range labels with robust error handling
#   range_labels <- tryCatch({
#     if (length(breaks) > 0) {
#       c(
#         paste0("\u2264 ", breaks[1]),  # ≤ first break
#         if (length(breaks) > 1) {
#           sapply(1:(length(breaks)-1), function(i) {
#             paste0(">", breaks[i], " - ", "\u2264", breaks[i+1])  # middle ranges
#           })
#         },
#         paste0("> ", breaks[length(breaks)])  # > last break
#       )
#     } else {
#       "All values"
#     }
#   }, error = function(e) {
#     warning("Error creating range labels: ", e$message)
#     rep("Range error", nrow(df))
#   })
#   
#   # 4. Build output dataframe with safe assignments
#   result <- data.frame(Outcome = outcome, stringsAsFactors = FALSE)
#   
#   for (i in seq_len(nrow(df))) {
#     tryCatch({
#       result[[paste0("N", i)]] <- df$group_n[i]
#       result[[paste0("Risk/Mscore", i)]] <- round(df$risk_mscore[i], 4)
#       result[[paste0("Thres", i)]] <- if (i <= length(range_labels)) range_labels[i] else "N/A"
#     }, error = function(e) {
#       warning("Error assigning values for group ", i, ": ", e$message)
#       result[[paste0("N", i)]] <<- NA
#       result[[paste0("Risk/Mscore", i)]] <<- NA
#       result[[paste0("Thres", i)]] <<- "Error"
#     })
#   }
#   
#   return(result)
# }


outdata <- function(data, outcome) {
  # Load required package
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Please install the 'dplyr' package first")
  }
  
  # 1. Safely process group statistics
  df <- tryCatch({
    data.frame(
      group = data$groups,
      estimate = data$estimates
    ) %>%
      dplyr::group_by(group) %>%
      dplyr::summarise(
        group_n = dplyr::n(),
        risk_mscore = mean(estimate, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::arrange(group)
  }, error = function(e) {
    stop("Error processing group statistics: ", e$message)
  })
  
  # 2. Process breakpoints safely
  breaks <- tryCatch({
    sort(unique(data$brkPoints))
  }, error = function(e) {
    warning("Error processing breakpoints: ", e$message)
    numeric(0)
  })
  
  # 3. Create range labels with robust error handling
  # 3. Create range labels (including last group as "> last breakpoint")
  range_labels <- tryCatch({
    if (length(breaks) == 0) {
      "All values"
    } else {
      labels <- c(
        paste0("\u2264 ", breaks[1])  # ≤ first break
      )
      if (length(breaks) > 1) {
        for (i in 2:length(breaks)) {
          labels <- c(labels,
                      paste0("> ", breaks[i - 1], " - ", "\u2264 ", breaks[i]))  # middle ranges
        }
      }
      # Add the final range: > last breakpoint
      labels <- c(labels, paste0("> ", breaks[length(breaks)]))
      labels
    }
  }, error = function(e) {
    warning("Error creating range labels: ", e$message)
    rep("Range error", nrow(df))
  })
  
  # 4. Build output dataframe with safe assignments
  result <- data.frame(Outcome = outcome, stringsAsFactors = FALSE)
  
  for (i in seq_len(nrow(df))) {
    tryCatch({
      result[[paste0("N", i)]] <- df$group_n[i]
      result[[paste0("Risk/Mscore", i)]] <- round(df$risk_mscore[i], 4)
      result[[paste0("Thres", i)]] <- if (i <= length(range_labels)) range_labels[i] else "N/A"
    }, error = function(e) {
      warning("Error assigning values for group ", i, ": ", e$message)
      result[[paste0("N", i)]] <<- NA
      result[[paste0("Risk/Mscore", i)]] <<- NA
      result[[paste0("Thres", i)]] <<- "Error"
    })
  }
  
  return(result)
}