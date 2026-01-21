<<<<<<< HEAD


library(tidyverse)
library(rmarkdown)


UploadDate <- "2025-10-31"  # Define the upload date

output_path <- paste0("D:/Users/williams_pj/Documents/Analysis/ReMAPP/Aim1/",UploadDate , "/ReMAPP_Aim1_Report_", UploadDate ,".pdf")

rmarkdown::render(
  input = "~/REMAPP-Public/ReMAPP-Aim1/2_ReMAPP_Aim1_Reports_current.Rmd",
  output_format = "pdf_document",
  output_file = output_path
)

=======


library(tidyverse)
library(rmarkdown)


UploadDate <- "2025-10-31"  # Define the upload date

output_path <- paste0("D:/Users/williams_pj/Documents/Analysis/ReMAPP/Aim1/",UploadDate , "/ReMAPP_Aim1_Report_", UploadDate ,".pdf")

rmarkdown::render(
  input = "~/REMAPP-Public/ReMAPP-Aim1/2_ReMAPP_Aim1_Reports_current.Rmd",
  output_format = "pdf_document",
  output_file = output_path
)

>>>>>>> 045b86a432752d5e21150935b2d614559ea27c90
