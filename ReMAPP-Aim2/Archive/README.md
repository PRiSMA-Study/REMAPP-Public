# ReMAPP Aim 2

## Description

Isotonic regression on hemoglobin and outcomes. 

#### :pushpin: Updated on 2025-01-21
#### :pushpin: Originally drafted by: Xiaoyan Hu (xyh@gwu.edu)

## File structure

**`1.Data prep.R`**
1. Load data: All the stacked and outcome data are loaded in this step. Also select the variables needed for analysis.
2. hb data: Prepare all hb data, including adjust hb with smoke and altitude, find extreme value at different conditions.
3. Maternal outcome data: merge maternal outcome data and corresponding hb data, remove NA rows.
4. Infant outcome data: merge infant outcome data and corresponding hb data, remove NA rows.
5. Save data for analysis.
6. Heatmap data: these codes are not related to isotonic regression. Mainly to visualize data collection status for each outcome.

**`2.run_b_asph.R to 2.run_b_stillbirth20.R`**
Run spline and isotonic for each binary outcome. 
Depending on the computer's capabilities, multiple codes can run simultaneously. 

**`2.run_c_dpr_score.R to 2.run_c_ftg_score.R`**
Run spline and isotonic for each continuous outcome. 
Depending on the computer's capabilities, multiple codes can run simultaneously.

**`3.plot_b_asph.R to 3.plot_b_stillbirth20.R`**
Plot spline and isotonic line for each binary outcome. 

**`3.plot_c_dpr_score.R to 3.plot_c_ftg_score.R`**
Plot spline and isotonic line for each continuous outcome. 

**`4.ReMAPP-Aim2-Plots-Trims.RMD`**
Generate isotonic panel plots.

**`5.ReMAPP-Aim2-Plots-Heatmap.RMD`**
Generate heatmap for current outcome collection status.

**_iso_code_**

**`spling_binary.R`**
Generate splines to model nonlinear relationships between the predictor and the binary response variable, using a Generalized Linear Model (GLM) with a logit link function, while selecting optimal knots for the spline.

**`spline_continuous.R`**
Generate splines to model nonlinear relationships between the predictor and the continuous response variable, using a Generalized Linear Model (GLM) with an identity link function, while selecting optimal knots for the spline.

**`iso_binary.R`** 
Implement a flexible isotonic regression approach using a Generalized Linear Model (GLM) with a logit link function, incorporating mixed effects (fixed and random), splines, and statistical significance testing to identify optimal partition points in the response variable, while accounting for covariates and random effects.

**`Iso_continuous.R`** 
Implement a flexible isotonic regression approach using mixed-effects models (linear mixed-effects with random effects), to detect multiple change points in a continuous response variable, while adjusting for covariates and random effects, and determining optimal partition points based on statistical significance.

**`plot_binary.R`** 
Plot function for binary outcome.

**`plot_continuous.R`**
Plot function for continuous outcome.

**`outdata.R`**
Function for isotonic regression results.



