


library(tidyverse)
library(rmarkdown)


UploadDate <- "2026-01-30"  # Define the upload date

output_path <- paste0("D:/Users/williams_pj/Documents/Analysis/ReMAPP/Aim1/",UploadDate , "/ReMAPP_Aim1_Report_", UploadDate ,".pdf")

rmarkdown::render(
  input = "~/REMAPP-Public/ReMAPP-Aim1/2_remapp_aim_one_analysis.Rmd",
  output_format = "pdf_document",
  output_file = output_path
)

