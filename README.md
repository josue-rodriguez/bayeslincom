
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bayeslincom

<!-- badges: start -->

<!-- badges: end -->

The goal of bayeslincom is to …

## Installation

You can install the released version of bayeslincom from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("bayeslincom")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("josue-rodriguez/bayeslincom")
```

## Example: BGGM

This is a basic example which shows you how to solve a common problem:

``` r
library(bayeslincom)
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

## Example: Data Frame

### MCMCpack

There are a variety of `R` packages that provide samples from the
posterior distribution. By placing the respective samples into a
`data.frame`, the **bayeslincom** can be used to test linear
combinations. Here is an example using **MCMCpack**.

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
#>      2.53    2.06   -1.53    6.58    0.11       0.89
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
