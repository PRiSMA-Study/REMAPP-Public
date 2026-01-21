

library(tidyverse)
library(rmarkdown)


UploadDate <- "2025-10-31"  # Define the upload date

output_path <- paste0("D:/Users/williams_pj/Documents/Analysis/ReMAPP/Aim2/",UploadDate , "/ReMAPP_Aim2_Plots_", UploadDate, "_v" ,".pdf")

rmarkdown::render(
  input = "~/REMAPP-Public/ReMAPP-Aim2/3_remapp_aim2_plots_report.Rmd",
  output_format = "pdf_document",
  output_file = output_path
)
