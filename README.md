
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bayeslincom

<!-- badges: start -->

<!-- badges: end -->

<!-- The goal of bayeslincom is to ... -->

## Installation

<!-- You can install the released version of bayeslincom from [CRAN](https://CRAN.R-project.org) with: -->

<!-- ``` r -->

<!-- install.packages("bayeslincom") -->

<!-- ``` -->

The development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("josue-rodriguez/bayeslincom")
```

## Example: **BGGM**

This is a basic example which shows you how to solve a common problem in
the `R` package **BGGM**:

``` r
library(bayeslincom)
## basic example code
```

## Example: **BBcor**

## Example: Data Frame

### **MCMCpack**

There are a variety of `R` packages that provide samples from the
posterior distribution. By placing the respective samples into a
`data.frame`, **bayeslincom** can be used to test linear combinations.
Here is an example using **MCMCpack**.

``` r
library(MCMCpack)

# data 
Y <- mtcars

# fit model
fit_bayes <- MCMCpack::MCMCregress(mpg ~ vs + hp, 
                                   data = Y, 
                                   mcmc = 100000)
                                   
# data frame
samps <- as.data.frame(fit_bayes)
                                   
# test hypothesis
test <- lin_comb(lin_comb = "vs - hp = 0", 
                 obj  = samps, 
                 cri_level = 0.95)                                 

# print results
fit_bayes

#> bayeslincom: linear combinations of posterior samples
#> ------ 
#> Call:
#> lin_comb.data.frame(lin_comb = lin_comb, obj = obj, cri_level = cri_level, 
#>    rope = rope)
#> ------ 
#> combination: vs + hp = 0 
#> ------ 
#> Posterior Summary:
#> 
#> Post.mean Post.sd Cred.lb Cred.ub Pr.less Pr.greater
#>      2.64    2.04   -1.38    6.65    0.09       0.91
#> ------ 
#> note:
#> Pr.less: Posterior probability less than zero
#> Pr.greater: Posterior probability less than zero
```

Note that the hypothesis could also be written as `vs = hp`.

#### Comparison to **multcomp**

The above can also be implemented with the `R` package **multcomp**.

``` r
library(multcomp)

# fit model
fit_lm <- lm(mpg ~ vs + hp, Y) 

# confidence interval
confint(
  multcomp::glht(fit_lm, linfct = "vs - hp == 0"), 
  level = 0.95
  )

#>   Simultaneous Confidence Intervals
#> 
#> Fit: lm(formula = mpg ~ vs + hp, data = Y)
#> 
#> Quantile = 2.0452
#> 95% family-wise confidence level
 

#> Linear Hypotheses:
#>              Estimate lwr     upr    
#> vs - hp == 0  2.6308  -1.3763  6.6378
```

Although the results are nearly identical, note that **bayeslincom** (1)
provides the posterior probability of a positive and negative
difference; and (2) is compatible with essentially all `R` packages for
Bayesian analysis.
