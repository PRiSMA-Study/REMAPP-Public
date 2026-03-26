

library(tidyverse)
library(rmarkdown)


UploadDate <- "2026-01-30"  # Define the upload date

output_path <- paste0("D:/Users/williams_pj/Documents/Analysis/ReMAPP/Aim2/",UploadDate , "/ReMAPP_Aim2_Plots_", UploadDate, "v.1" ,".pdf")

rmarkdown::render(
  input = "~/REMAPP-Public/ReMAPP-Aim2/3_remapp_aim_two_plots.Rmd",
  output_format = "pdf_document",
  output_file = output_path
)
