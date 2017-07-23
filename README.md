[![Travis-CI Build Status](https://travis-ci.org/nuffe/penfmbr.svg?branch=master)](https://travis-ci.org/nuffe/penfmbr)

## penfmbr
simple implementation of the penalized FamaMacBeth estimator*

## Installation
```r
install.packages("devtools")
devtools::install_github("nuffe/penfmbr")
```

## Example
```r
library(penfmbr)
library(stargazer)

# Load factors and returns data from the accompanying test folder
path = paste(.libPaths()[1], "/penfmbr/tests/testthat", sep = "")
factors = read.csv(paste(path, "/factors.csv", sep = ""))[, -1]
returns = read.csv(paste(path, "/returns.csv", sep = ""))[, -1]

models = list(
  fmb_mkt = fmb(factors[, 1, drop = FALSE], returns),
  fmb_ff = fmb(factors[, c(1, 2, 3)], returns),
  fmb_ch = fmb(factors, returns),
  penfmb_mkt = penfmb(factors[, 1, drop = FALSE], returns),
  penfmb_ff = penfmb(factors[, c(1, 2, 3)], returns),
  penfmb_ch = penfmb(factors, returns)
)

stargazer(models,
          # (generic methods are added for vcov - computes shrinakge rate and/or GMM/Shanken standard errors)
          se = lapply(models, function(x) sqrt(diag(vcov(x)))),
          keep.stat = 'adj.rsq',
          type = 'text',
          dep.var.caption = '',
          dep.var.labels = '',
          star.cutoffs = NA,
          column.labels = c('FamaMacBeth', 'Penalized FamaMacBeth'),
          column.separate = c(3, 3),
          notes = c('FamaMacBeth: Shanken standard errors in parenthesis',
                    'Penalized FamaMacBeth: bootstrapped shrinkage rate in parenthesis'),
          notes.append = FALSE,
          notes.align = 'l')

# ==================================================================================
#                                                                                   
#                          FamaMacBeth                   Penalized FamaMacBeth      
#                 (1)          (2)          (3)        (4)        (5)        (6)    
# ----------------------------------------------------------------------------------
# Mkt            -0.697       -0.795       0.025      -0.697     -0.795     -0.747  
#               (0.483)      (0.355)      (0.437)    (0.000)    (0.000)    (0.000)  
#                                                                                   
# SMB                         0.141        0.152                 0.141      0.143   
#                            (0.143)      (0.144)               (0.000)    (0.000)  
#                                                                                   
# HML                         0.426        0.447                 0.426      0.429   
#                            (0.135)      (0.137)               (0.000)    (0.000)  
#                                                                                   
# Mom                                      2.308                            0.000   
#                                         (0.792)                          (0.710)  
#                                                                                   
# Constant       1.474        1.350        0.574      1.474      1.350      1.299   
#               (0.442)      (0.290)      (0.376)    (0.000)    (0.000)    (0.000)  
#                                                                                   
# ----------------------------------------------------------------------------------
# Adjusted R2    0.165        0.674        0.723      0.165      0.673      0.657   
# ==================================================================================
# Note:       FamaMacBeth: Shanken standard errors in parenthesis                   
#             Penalized FamaMacBeth: bootstrapped shrinkage rate in parenthesis
```

> *S. Bryzgalova, Spurious Factors in Linear Asset Pricing Models (2016)
