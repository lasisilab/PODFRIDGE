# Instructions: IPUMS Fertility & Sibling Distribution Mirrored Plots

These instructions describe how to reproduce two mirrored frequency bar plots showing (1) the distribution of children ever born and (2) the derived sibling distribution, by race, census year, and age group. The data comes from IPUMS USA Census microdata.

---

## 1. IPUMS Data Extract Specification

### Source
**IPUMS USA** — https://usa.ipums.org/usa/

**Citation:** Steven Ruggles, Sarah Flood, Matthew Sobek, Daniel Backman, Annie Chen, Grace Cooper, Stephanie Richards, Renae Rogers, and Megan Schouweiler. IPUMS USA: Version 14.0 [dataset]. Minneapolis, MN: IPUMS, 2023. https://doi.org/10.18128/D010.V14.0

### Samples to select
- 1960 Census (1% sample)
- 1970 Census (1% sample)
- 1980 Census (5% sample)
- 1990 Census (5% sample)

### Variables to include in the extract

| Variable  | IPUMS name | Description |
|-----------|-----------|-------------|
| `SEX`     | SEX       | Sex of respondent |
| `AGE`     | AGE       | Age of respondent |
| `BIRTHYR` | BIRTHYR   | Year of birth |
| `RACE`    | RACE      | Race of respondent |
| `CHBORN`  | CHBORN    | Children ever born (asked of women; coded categorically in IPUMS) |
| `NCHILD`  | NCHILD    | Number of own children in household |
| `NCHLT5`  | NCHLT5    | Number of own children under age 5 in household |

### Universe
Women (SEX == Female). The analysis further filters to women aged 40+.

---

## 2. Data Processing / Recoding

After downloading the IPUMS extract, apply these steps in R:

### 2a. Filter to women aged 40+

```r
df <- raw_data %>%
  filter(SEX == "Female", AGE >= 40)
```

### 2b. Recode CHBORN to numeric `chborn_num`

IPUMS codes `CHBORN` as a categorical/labeled variable. Convert it to a plain integer count of children ever born. The exact recoding depends on the IPUMS coding scheme for your extract — typically:
- 0 = "No children" or "N/A"
- 1 = "1 child"
- 2 = "2 children"
- ... etc.
- 12+ should be capped (see below)

The resulting numeric column should be called `chborn_num`.

### 2c. Recode RACE

Keep only two groups and recode to these exact labels:
- `"White"`
- `"Black/African American"`

```r
df <- df %>%
  filter(RACE %in% c("White", "Black/African American")) %>%
  mutate(RACE = factor(RACE, levels = c("White", "Black/African American")))
```

### 2d. Create age range bins

```r
df <- df %>%
  mutate(
    AGE_RANGE = case_when(
      AGE >= 70 ~ "70+",
      AGE >= 60 ~ "60-69",
      AGE >= 50 ~ "50-59",
      AGE >= 40 ~ "40-49"
    )
  )
```

### 2e. Cap children at 12+

For the CHBORN factor label (used in display), cap at "12+ children":

```r
df <- df %>%
  mutate(
    CHBORN = factor(case_when(
      chborn_num == 0  ~ "No children",
      chborn_num == 1  ~ "1 child",
      chborn_num == 2  ~ "2 children",
      chborn_num == 3  ~ "3 children",
      chborn_num == 4  ~ "4 children",
      chborn_num == 5  ~ "5 children",
      chborn_num == 6  ~ "6 children",
      chborn_num == 7  ~ "7 children",
      chborn_num == 8  ~ "8 children",
      chborn_num == 9  ~ "9 children",
      chborn_num == 10 ~ "10 children",
      chborn_num == 11 ~ "11 children",
      chborn_num >= 12 ~ "12+ children"
    ), levels = c("No children", "1 child", "2 children", "3 children",
                  "4 children", "5 children", "6 children", "7 children",
                  "8 children", "9 children", "10 children", "11 children",
                  "12+ children"), ordered = TRUE)
  )
```

### 2f. Final columns needed

The processed data frame should have at minimum: `YEAR`, `SEX`, `AGE`, `BIRTHYR`, `RACE`, `CHBORN`, `AGE_RANGE`, `chborn_num`

---

## 3. Plot 1: Distribution of Number of Children (Mirrored Bar Chart)

### Data preparation

Calculate proportions within each (YEAR, RACE, AGE_RANGE) group, then mirror the White proportions to the negative side:

```r
df_proportions <- df %>%
  group_by(YEAR, RACE, AGE_RANGE, chborn_num) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(YEAR, RACE, AGE_RANGE) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()

df_mirror <- df_proportions %>%
  mutate(proportion = if_else(RACE == "White", -proportion, proportion))
```

### Plot code

```r
my_colors <- colorRampPalette(c("#FFB000", "#F77A2E", "#DE3A8A", "#7253FF", "#5E8BFF"))(13)

child_plot <- ggplot(df_mirror, aes(x = chborn_num, y = proportion, fill = as.factor(chborn_num))) +
  geom_col(aes(alpha = RACE)) +
  geom_hline(yintercept = 0, color = "black", size = 0.5) +
  facet_grid(AGE_RANGE ~ YEAR, scales = "free_y") +
  coord_flip() +
  scale_y_continuous(
    labels = function(x) abs(x),
    limits = function(x) c(-max(abs(x)), max(abs(x)))
  ) +
  scale_x_continuous(breaks = 0:12, labels = c(0:11, "12+")) +
  scale_fill_manual(values = my_colors) +
  scale_alpha_manual(values = c("White" = 0.7, "Black/African American" = 1), guide = "none") +
  labs(
    title = "Distribution of Number of Children by Census Year, Race, and Age Range",
    x = "Number of Children",
    y = "Proportion",
    fill = "Number of Children",
    caption = "White population shown on left (negative values), Black/African American on right (positive values)\nProportions normalized within each age range, race, and census year\nThe category '12+' includes families with 12 or more children."
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, hjust = 0.5),
    axis.text.y = element_text(size = 8),
    strip.text = element_text(size = 10),
    legend.position = "none",
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  )
```

### Key design elements
- **Facets:** rows = AGE_RANGE, columns = YEAR
- **Orientation:** `coord_flip()` — children count on the y-axis, proportion on x-axis
- **Mirror:** White proportions are negated so they extend to the left; Black/AA extends to the right
- **Alpha:** White bars at 0.7 opacity, Black/AA bars at full opacity
- **Legend:** hidden (`legend.position = "none"`)

---

## 4. Plot 2: Distribution of Number of Siblings (Mirrored Bar Chart)

### Deriving sibling counts from the children-born data

The sibling distribution is derived from the fertility data. The key idea: if a mother has `chborn_num` children, each of those children has `chborn_num - 1` siblings. To get the correct *individual-level* sibling frequency, each mother's row is weighted by the number of children she had.

```r
df2 <- df %>%
  dplyr::select(RACE, YEAR, AGE_RANGE, chborn_num) %>%
  mutate(
    n_siblings = chborn_num - 1,
    sibling_freq = ifelse(chborn_num != 1, chborn_num * 1, 1)
  )
```

**Explanation of `sibling_freq`:** Each mother represents `chborn_num` children in the next generation. So a mother with 5 children contributes 5 individuals who each have 4 siblings. The weight is `chborn_num` (except mothers with 1 child contribute 1 individual with 0 siblings).

### Aggregate sibling data

```r
df_siblings <- df2 %>%
  group_by(YEAR, RACE, AGE_RANGE, n_siblings) %>%
  summarise(sibling_count = sum(sibling_freq), .groups = "drop")

df_sibling_proportions <- df_siblings %>%
  group_by(YEAR, RACE, AGE_RANGE) %>%
  mutate(proportion = sibling_count / sum(sibling_count)) %>%
  ungroup()

df_sibling_mirror <- df_sibling_proportions %>%
  mutate(proportion = if_else(RACE == "White", -proportion, proportion))
```

### Plot code

```r
my_colors <- colorRampPalette(c("#FFB000", "#F77A2E", "#DE3A8A", "#7253FF", "#5E8BFF"))(13)

sibling_plot <- ggplot(
  data = df_sibling_mirror %>% filter(n_siblings != -1),
  aes(x = n_siblings, y = proportion, fill = as.factor(n_siblings))
) +
  geom_col(aes(alpha = RACE)) +
  geom_hline(yintercept = 0, color = "black", size = 0.5) +
  facet_grid(AGE_RANGE ~ YEAR, scales = "free_y") +
  coord_flip() +
  scale_y_continuous(
    labels = function(x) abs(x),
    limits = function(x) c(-max(abs(x)), max(abs(x)))
  ) +
  scale_x_continuous(breaks = 0:11, labels = c(0:10, "11+")) +
  scale_fill_manual(values = my_colors) +
  scale_alpha_manual(values = c("White" = 0.7, "Black/African American" = 1)) +
  labs(
    title = "Distribution of Number of Siblings by Census Year, Race, and Age Range",
    x = "Number of Siblings",
    y = "Proportion",
    fill = "Number of Siblings",
    caption = "White population shown on left (negative values), Black/African American on right (positive values)\nProportions normalized within each age range, race, and census year\nThe category '11+' includes individuals with 11 or more siblings."
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, hjust = 0.5),
    axis.text.y = element_text(size = 8),
    strip.text = element_text(size = 10),
    legend.position = "none",
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  )
```

### Key differences from Plot 1
- X-axis is `n_siblings` (0 to 11+) instead of `chborn_num` (0 to 12+)
- Rows with `n_siblings == -1` (from mothers with 0 children) are filtered out
- Sibling counts are capped at "11+" (since 12+ children → 11+ siblings)

---

## 5. R Package Dependencies

```r
library(dplyr)
library(tidyverse)  # includes ggplot2, tidyr, readr, etc.
library(viridis)
library(scales)
```

---

## 6. Color Palette Reference

The 13-color gradient ramp used for both plots:

```r
my_colors <- colorRampPalette(c("#FFB000", "#F77A2E", "#DE3A8A", "#7253FF", "#5E8BFF"))(13)
```

This produces a warm-to-cool gradient: gold → orange → magenta → purple → blue.

---

## 7. Figure Dimensions

The original plots use `fig.width = 9, fig.height = 7` in the R Markdown chunk options.
