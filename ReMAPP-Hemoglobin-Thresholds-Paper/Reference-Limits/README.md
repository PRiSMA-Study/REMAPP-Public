# Reference Limits

This directory contains the public code used to reproduce the **primary reference-limit analysis** for the ReMAPP maternal hemoglobin thresholds paper.

## Analysis

Maternal hemoglobin reference limits are estimated among pregnancies meeting the ReMAPP healthy cohort criteria.

The analysis uses **fractional polynomial regression (FPR)** to model the relationship between gestational age and hemoglobin. 

The fitted model is used to estimate the **5th, 50th, and 95th hemoglobin centiles** for each exact gestational week. 

## Folder structure

```text
Reference-Limits/
├── README.md
└──  reference_limits.R
