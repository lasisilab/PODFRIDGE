# Predict racial proportions in the M&T forensic database

**Author:** Junhui He  
**Date:** 2025-06-15

## 1 Objective

To analyse the comparison between the M&T forensic DNA database and census database, we have to attain the racial breakdown (e.g., the proportions of black, white people) of the state-level forensic database across the US. However, the breakdown is given in only 7 states, which are California, Florida, Indiana, Maine, Nevada, South Dakota and Texas. Thus, we need to develop a statistical model to estimate the proportions of black and white Americans in the remaining 43 states. Specifically, we focus on the differences of the racial breakdown between the forensic database and census database.

## 2 Binomial Logistic Regression Model

In the underlying work, we establish binomial regression models to make the estimations. To make the main findings interpretable, we will give a brief introduction of this common model. *(If you don't care about the details of the binomial regression model, just feel free to skip this part.)*

Firstly, we introduce the concept of binomial distributions. Suppose in each observation, an event has two possible states, success or failure, and the probability of success is defined as $p$. Then in $n$ observations, the number of success $y\in \{0,1,\ldots, n\}$ follows a binomial distribution $\text{B}(n,p)$. A special case of the binomial distribution is $n=1$, where the distribution is the simple Bernoulli distribution. The success probability $p$ determines the characteristics of the binomial distribution. An important property is that the expectation is given by $np$.

Now we can delve into the binomial regression. Let the predictor be $x$ and the response variable be $y$, we assume:

$$y|x\sim \text{B}(n, g(x^\top \beta)),$$

where $\beta$ is the coefficients, and $g$ is a link function taking values in $[0,1]$. The popular choices of $g$ include Logit and Probit functions. Here we choose the Logit link function for its interpretability, which is given by:

$$g:\mathbb{R}\rightarrow [0,1], \quad g(z)=\frac{\exp(z)}{\exp(z)+1}.$$

To fit this model, we will estimate the linear coefficients $\hat{\beta}$, thereby we can predict the success probability $\hat{p}=g(x_*^\top \hat{\beta})$ on a new point $x_*$, which can be considered as the proportion of success events given $x_*$.

## 3 Data and Model Setting

In this section, we demonstrate the response variable and predictors used in the binomial regression, and give the concrete model equation.

### 3.1 Response variable:

- The total number of people for each state in the M&T forensic database.
- The number of each race for each state in the M&T database.

### 3.2 Predictor:

- The proportion of black and white people for each state in the census database.
- The proportion of black and white people of the prison population for each state.

### 3.3 Stick-breaking:

We divide the people of the US into three categories, black + white + other. Then we need to estimate the complete racial breakdown for three categories in each state. To ensure the sum of these percents is equal to 1, we use a simple stick-breaking technique. That is, we separately estimate the percent of white people in all people $p_{1}$ and the percent of black people in non-white people $p_{2}$. Then the racial breakdown is given by:

$$\pi_{white}=p_1, \quad \pi_{black}=(1-p_1)p_2,\quad \pi_{other}=(1-p_1)(1-p_2).$$

We run 2 binomial regression models to predict $p_1$ and $p_2$ for each state.

### 3.4 Model equation

Compared to the predictors used by Hanna, we remove the racial indicator and the interaction between the racial indicator and the census/prison proportion to avoid colinearity, which leads to a singular problem for linear regression. Therefore, the model equations are defined as:

$$ \frac{white}{all} = g(\beta_{00} + \beta_{01}*white_{census}+ \beta_{02}*white_{prison}), $$

$$ \frac{black}{nonwhite} = g(\beta_{10} + \beta_{11}*black_{census} + \beta_{12}*black_{prison}). $$

### 3.5 Coefficient interpretability

To interpret those coefficients, we simply denote the binomial logistic regression model as:

$$p = g(\beta_0 + \sum_{j=1}^p\beta_j x_{j}),$$

where $g$ is the logit link function and $x_j$ for $1,\ldots,p$ are covariates. We consider an odds as $p/(1-p)$, which is the ratio of success probabilities versus failure probabilities. Thus, using the logit link function, the model equation can be written as:

$$\log(\frac{p}{1-p})=\beta_0 + \sum_{j=1}^p\beta_j x_{j}.$$

Therefore, we can interpret the coefficients $\beta_j$ as the increase in the log odds for every unit increase in $x_j$.

## 4 Model Evaluation

*(If you just want to read the answers to the key questions, please skip this part and go to the next section directly.)*

- The model is fit using binomial regression for both White and Black Americans.
- Confidence intervals and p-values are calculated for the predictions.
- Goodness-of-fit is assessed by comparing fitted values to ground truth for the 7 states with available data.

## 5 Main Findings

- Black Americans are significantly overrepresented while White Americans are underrepresented in the M&T forensic database compared to census representation.
- Side-by-side pie charts for each state show the racial composition according to the census (left) versus the estimated racial composition of CODIS (right) for each state.
- Absolute and relative differences of racial proportions are visualized, showing Black/African Americans are significantly overrepresented in all states and White Americans are underrepresented in most states in the M&T forensic database compared to census representation.

## 6 Statistical Inference

- The normal approximation for the log odds is used for hypothesis testing and confidence interval construction.
- For each state, a hypothesis test is performed for the difference of white proportions between Census and CODIS.
- Confidence intervals for the fitted probability are calculated using binomial regression.
- For Black Americans, a Bonferroni method is used to construct the confidence interval.
- As the number of population in each state with available data is very large, the estimated standard errors are very small, causing the interval widths to be almost zero compared to the point estimation. This also explains why the p-values for the above hypothesis testing are nearly zero.
