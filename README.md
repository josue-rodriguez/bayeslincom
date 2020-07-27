
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

The **BBcor** package provides Bayesian bootstrapped correlations and
partial correlations. In the following, the test is whether one
correlation is larger than the sum of two correlations. This is
implemented with

``` r
library(BBcor)

# data
Y <- mtcars[,c("mpg", "wt", "hp")]

# fit model
fit <- bbcor(Y, method = "spearman")

# print
fit

#>           mpg         wt         hp
#> mpg  1.0000000 -0.8842590 -0.8946683
#> wt  -0.8842590  1.0000000  0.7754185
#> hp  -0.8946683  0.7754185  1.0000000
```

Next the hypothesis is written out. In this case

``` r
hyp <- "-1*mpg--wt - (-1*mpg--hp + wt--hp) = 0"
```

Note that each relation is multiplied by the sign. This ensures the
magnitude is being compared. Next the hypothesis is tested as follows

``` r
lin_comb(hyp, 
         obj = fit, 
         cri_level = 0.95)

#> bayeslincom: linear combinations of posterior samples
#> ------ 
#> Call:
#> lin_comb.bbcor(lin_comb = lin_comb, obj = obj, cri_level = cri_level, 
#>     rope = rope)
#> ------ 
#> combination: -1*mpg--wt - (-1*mpg--hp + wt--hp) = 0 
#> ------ 
#> Posterior Summary:
#> 
#>  Post.mean Post.sd Cred.lb Cred.ub Pr.less Pr.greater
#>      -0.78    0.09   -0.93   -0.59       1          0
#> ------ 
#> note:
#> Pr.less: Posterior probability less than zero
#> Pr.greater: Posterior probability less than zero
```

In this case, the sum of the relations, `mpg--hp` and `wt-hp`, are
larger than the relation `mpg--wt`, with a posterior probability of 1.

This example also demonstrates a key contribution of **bayeslincom**,
i.e., testing linear combinations among dependent correlations. This is
the only Bayesian (and possibly in general) implementation in `R`.

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
fit_lm <- lm(mpg ~ vs + hp, 
             data = Y) 

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
