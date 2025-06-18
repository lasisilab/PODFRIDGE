# Regression to predict CODIS database proportions

**Author:** Hannah Van Wyk  
**Date:** 2025-06-15

## Overview

To estimate the number of people in CODIS by state and race, we use several data sources:
- Direct data on the number of Black and White people in CODIS for 7 states (Murphy & Tong, 2020)
- The number of people of each racial group in prison in each state (Klein et al., 2023)
- The number of people in the State DNA Indexing System (SDIS) and the National DNA Indexing System (NDIS)

We separate states into three categories:
1. States with Murphy & Tong data (no estimation needed)
2. States with NDIS and SDIS data but not in Murphy & Tong
3. States with only NDIS data

## Data Preparation

- Census, prison, and CODIS data are merged by state.
- Proportions of Black and White populations and incarcerated individuals are calculated.
- For states without direct data, regression models are used to estimate racial proportions in CODIS.

## Regression Model for Racial Composition

A regression model is fit using the Murphy & Tong states to estimate the proportion of Black and White people in CODIS as a function of:
- Census proportion of each race
- Prison proportion of each race
- Race indicator (Black/White)
- Interactions between race and census/prison proportions

The model equation:

$$Proportion_{race} = \beta_0 + \beta_1census_{proportion} + \beta_2prison_{proportion} + \beta_3race + \beta_4race*census_{proportion} + \beta_5race*prison_{proportion}$$

## SDIS Regression for States with Only NDIS Data

For states with only NDIS data, additional regression models are used to estimate the number of people in SDIS:

- For arrestees:
  $$N_{arrestee} = \beta_0 + \beta_1census_{black} + \beta_2census_{white} + \beta_3prison_{black} + \beta_4prison_{white} + \beta_5NDIS_{arrestees} $$
- For offenders:
  $$N_{offender} = \beta_0 + \beta_1census_{black} + \beta_2census_{white} + \beta_3prison_{black} + \beta_4prison_{white} + \beta_5NDIS_{offenders} $$

## Final Estimates

- For each state, the number of Black and White people in CODIS is estimated using the best available data and regression predictions.
- Results are visualized with:
  - Plots comparing estimated and actual values for states with data
  - Maps showing the difference between CODIS and census racial proportions
  - Pie charts for each state comparing census and CODIS racial composition

## Key Findings

- The regression model fits the Murphy & Tong data well ($R^2 = 0.93$).
- Black Americans are overrepresented and White Americans are underrepresented in CODIS compared to census proportions in most states.
- The workflow provides state-level, data-driven estimates of racial composition in CODIS using the most accurate available methods.

## References

1. Murphy, Erin, and Jun H. Tong. "The racial composition of forensic DNA databases." Calif. L. Rev. 108 (2020): 1847.
2. Klein, Brennan, et al. "COVID-19 amplified racial disparities in the US criminal legal system." Nature 617.7960 (2023): 344-350.
3. https://le.fbi.gov/science-and-lab/biometrics-and-fingerprints/codis/codis-ndis-statistics
