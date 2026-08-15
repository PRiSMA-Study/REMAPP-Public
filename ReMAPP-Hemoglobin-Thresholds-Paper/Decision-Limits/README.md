
# Decision Limits

This directory contains the public code used to reproduce the **decision-limit analyses** for the ReMAPP maternal hemoglobin thresholds paper.

## Analysis

Decision limits are evaluated by examining the relationship between maternal hemoglobin and maternal and infant outcomes.

The primary analysis uses:

- **restricted cubic spline regression** to characterize non-linear Hb–outcome relationships; and
- **isotonic regression** to characterize monotonic changes in outcome risk and estimate candidate decision limits.

The analysis framework can be applied to the maternal and infant outcomes included in the manuscript by specifying the corresponding analysis dataset, outcome variable, hemoglobin variable, and outcome type.

## Folder structure

```text
Decision-Limits/
├── README.md
├── decision_limits.R
└──  iso-code/
  ├── ...
  └── ...
