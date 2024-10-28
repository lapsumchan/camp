# Censoring-based analysis of microbiome proportions (CAMP)

### Overview
This repository provides a demonstration on how to use the `R` package `camp`

# System Requirements

## Software Requirements

### OS Requirements

The source code has been tested on Microsoft's Windows 10 operating system and Linux (Ubuntu 18.04). The source code should be compatible with Windows, Mac, and Linux operating systems.

Before using `camp` package, users should have `R` version 4.3.0 or higher.

### Installation  

First, install `camp` from GitHub using `devtools`:  

    # Install devtools if not already installed
    # install.packages("devtools") 
    devtools::install_github("lapsumchan/camp")
    
Installation should complete within a couple of minutes on a standard machine.

# Demo

To get started, load the necessary packages:

```
library(camp)
library(survival)
library(MCMCpack)
```

Simulated data can be generated using the following code:

```
seed <- 1
n1 <- 100
n2 <- 100
p <- 1000
signal.prop <- 0.8
sum.const <- signal.prop
a <- 353

camp.dat <- GenData(seed, n1, n2, p, signal.prop, sum.const, a)
```

`camp.dat` is a list which contains an (`n x p`) OTU table and a length `n` vector of group labels.

To run camp, use the OTU table and group labels as inputs:
```
camp.res <- camp(camp.dat$OTU, camp.dat$Group)
head(camp.res)
```

# Output

Depending on the number of covariates included, the output will either be a length `p` vector of p-values or a (`p x q`) data frame of p-values, where `q` is the number of covariates. In this case, since the covariate `cov` only involves the group label (length `n` vector), the output of `camp` is a length `n` vector of p-values.

```
> head(camp.res)
[1] 3.710057e-03 1.110223e-16 3.561546e-05 3.406422e-03 1.208696e-05 1.957646e-04
```

# Citation

If you find `camp` useful, please cite:
> Chan, Lap Sum, and Gen Li. ["Zero is not absence: censoring-based differential abundance analysis for microbiome data."](https://academic.oup.com/bioinformatics/article/40/2/btae071/7603976) Bioinformatics
> 40.2 (2024): btae071.
