# ReMAPP Aim 3

## Description

To describe the underlying contributing factors for anemia during pregnancy. 

#### :pushpin: Updated on 2024-12-27
#### :pushpin: Originally drafted by: Yipeng Wei (yipeng_wei@gwu.edu)

## Part A

**Sample** : Aim 2 cohort of 1200-2000 women at study enrollment per site

**Method** : 1. Relative risk: For non-distal risk factors, Relative Risk (RR) is derived from univariate GEE (log-binomial) model adjusted for site and distal risk factors (demographic). For distal risk factors, Relative Risk (RR) is derived from multivariate GEE model. If log-binomial GEE model fails to converge, modified Poisson GEE model will be implemented.
2. Population attributable risk (PAR) : PAR from the perspective of causal effect using the stdReg R package. 
Reference: Sjölander, A. (2018). Estimation of causal effect measures with the R-package stdReg. European journal of epidemiology, 33(9), 847-858.

## Part B
**Sample** : Aim 3 sample of 300 women per clinical site (100 from each trimester of pregnancy)

**Method** : Relative risk: Relative Risk (RR) is derived from log-binomial generalized linear mixed effects model (link=log, distribution=binomial) adjusted for site (random intercept). If the model fails to converge, modified Poisson generalized linear mixed effects model (link=log, distribution=poisson) will be implemented.
