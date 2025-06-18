# Detailed Sources and Summary for Forensic DNA Database Estimates

**Authors:** Stella BooydeGraaff & Hannah Van Wyk  
**Date:** 2025-06-15

## Murphy & Tong 2020 Freedom of Information Act (FOIA) Appendix Data

Murphy & Tong (2020) provided information on the racial composition of NDIS for each state. Their data includes:
- Proportions of racial groups in CODIS (NDIS) for each state
- Proportions of racial groups in the census for each state
- Over- or under-representation of each racial group in CODIS compared to census

Key steps:
- Data is read from `murphy_foia_cleaned.csv`
- CODIS and census data are separated and harmonized
- Over/under-representation is calculated as the difference between CODIS and census proportions

## State DNA Indexing System (SDIS)

SDIS is a state-level database of DNA profiles. Data is separated into offenders and arrestees, and is publicly available for some states. For each state, the following are reported:
- Total number of profiles
- Number of arrestees and offenders
- Source of the data

Key steps:
- Data is read from `SDIS.csv`
- Missing totals are imputed as the sum of arrestees and offenders
- Proportion of offenders is calculated

## National DNA Index System (NDIS)

NDIS is the national-level CODIS database. All NDIS database sizes are publicly available. For each state, the following are reported:
- Total number of profiles
- Number of arrestees and offenders

Key steps:
- Data is read from `NDIS.csv`
- Proportion of offenders is calculated

## Prison Data

The number of incarcerated people of each racial group for each state is compiled from various sources (Klein et al., 2023; Vera Institute for Michigan). For each state, the following are reported:
- Total incarcerated
- Incarcerated by race (White, Black, Hispanic, Native American, Asian)

Key steps:
- Data is read from `populations_states.csv`
- Data is filtered for 2022 and cleaned
- Michigan's Black incarcerated population is estimated from a separate source

## Data Visualizations and Summaries

- Tables and plots summarize the racial and gender breakdowns in CODIS and census data
- Pie charts and bar charts illustrate the distribution of offender types, gender, and race
- Summary statistics highlight states with the highest and lowest Black population proportions

## References

1. Murphy, Erin, and Jun H. Tong. "The racial composition of forensic DNA databases." Calif. L. Rev. 108 (2020): 1847.
2. https://le.fbi.gov/science-and-lab/biometrics-and-fingerprints/codis/codis-ndis-statistics
3. Klein, Brennan, et al. "COVID-19 amplified racial disparities in the US criminal legal system." Nature 617.7960 (2023): 344-350.
4. https://www.vera.org/downloads/pdfdownloads/state-incarceration-trends-michigan.pdf
