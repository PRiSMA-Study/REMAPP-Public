# ReMAPP Hemoglobin Thresholds Paper

## Overview

This directory contains the public code used to reproduce the **primary analyses** for the ReMAPP maternal hemoglobin thresholds paper.

The study evaluates maternal hemoglobin (Hb) thresholds during pregnancy using two complementary approaches:

1. **Reference limits** — describing the distribution of hemoglobin across gestation among pregnancies meeting the ReMAPP healthy cohort criteria.
2. **Decision limits** — evaluating maternal hemoglobin in relation to maternal and infant outcomes.


# Reference Limits

The `Reference-Limits` directory contains the primary analysis used to estimate gestational-age-specific maternal hemoglobin reference limits.

## Analysis population

Reference limits are estimated using hemoglobin measurements from pregnancies meeting the ReMAPP healthy cohort criteria.

The dataset contains participant-level healthy cohort indicators and available hemoglobin measurements. 
The analysis code reconstructs the longitudinal hemoglobin data required for modeling and restricts the primary analysis to the specified healthy cohort.

Hemoglobin measurements between **6 and 40 completed weeks of gestation** are included in the primary FPR analysis.

## Fractional polynomial regression

Fractional polynomial regression (FPR) is used to characterize the non-linear relationship between gestational age and maternal hemoglobin.
The fitted distributions are used to estimate the:
* **5th centile**
* **50th centile**
* **95th centile**
for each exact gestational week from 6 through 40 weeks.

## Reference-limit outputs

The primary Reference Limits script generates:

* FPR model parameters and fractional polynomial powers;
* predicted 5th, 50th, and 95th Hb centiles by exact gestational week;
* number of Hb observations contributing to each gestational week;
* overall and trimester-specific summaries of the predicted centiles;
* WHO-style anemia severity thresholds derived from the estimated 5th centile; and
* the FPR centile curve figure.

The overall and trimester-specific summaries represent the **median of the gestation-specific predicted centile values** within each specified gestational-age interval.
Numerical outputs are combined into a single Excel workbook with separate worksheets for the major FPR results. The primary figure is saved separately.

# Decision Limits
The `Decision-Limits` directory contains the primary analyses used to evaluate maternal hemoglobin concentrations in relation to maternal and infant outcomes.

## Analysis approach
The Decision Limits analysis uses two complementary modeling approaches:

* **Restricted cubic spline regression** to characterize potentially non-linear relationships between maternal Hb and each outcome.
* **Isotonic regression** to characterize monotonic changes in outcome risk across Hb concentrations and estimate candidate outcome-based decision limits.

For each analysis, the user specifies the:
* analysis dataset;
* outcome variable;
* corresponding Hb variable;
* outcome type; and
* outcome label.

## Analysis periods
Where applicable and supported by the outcome-specific dataset, analyses are conducted for:

* all pregnancy Hb measurements combined;
* first trimester;
* second trimester;
* third trimester; and
* postpartum periods for applicable maternal outcomes.

Outcome-specific Hb variables are used so that the selected hemoglobin measurement corresponds to the relevant analysis period and outcome.

## Outcome types
Most maternal and infant outcomes are analyzed as **binary outcomes**.
The fatigue analysis uses the corresponding **continuous fatigue score** and continuous-outcome modeling functions.

## Decision-limit outputs
For each eligible outcome and analysis period, the analysis generates:

* the analytic sample size;
* restricted cubic spline model results;
* isotonic regression results;
* standardized decision-limit hemoglobin thresholds; and
* graphical displays of the Hb–outcome relationship.

# Public analysis datasets
The scripts in this directory begin with **ReMAPP analysis datasets**.
These datasets contain the derived variables required to reproduce the primary manuscript analyses and therefore represent the starting point for the public code.

# Reproducibility
The code in this directory has been prepared specifically for public reproducibility of the primary manuscript analyses.
Users should begin with the corresponding paper-ready ReMAPP dataset and run the appropriate script within either:
* `Reference-Limits`, for the primary FPR analysis; or
* `Decision-Limits`, for the primary outcome-based threshold analyses.
The resulting outputs can be compared with the corresponding tables, estimates, and figures reported in the manuscript.

## Software
Analyses were conducted using **R**.
Package requirements are specified within the individual analysis scripts. 
