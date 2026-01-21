

library(tidyverse)
library(rmarkdown)


UploadDate <- "2025-10-31"  # Define the upload date

output_path <- paste0("D:/Users/williams_pj/Documents/Analysis/ReMAPP/Manuscript/",UploadDate , "/ReMAPP_Aim1and2_Results", UploadDate ,".pdf")

rmarkdown::render(
  input = "~/REMAPP-Public/ReMAPP-Manuscript/ReMAPP_Aim1_Aim2_Results.Rmd",
  output_format = "pdf_document",
  output_file = output_path
)

