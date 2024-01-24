# CAMP

### Overview
This repository provides a demonstration on how to use the `R` package CAMP

# System Requirements

## Software Requirements

### OS Requirements

The source code has been tested on Microsoft's Windows 10 operating system and Linux (Ubuntu 18.04). The source code should be compatible with Windows, Mac, and Linux operating systems.

Before using the DrFARM source code, users should have `R` version 4.3.0 or higher, and several packages installed.

### Installation  
First, we need to install `camp`:
    install.packages("camp_1.0.tar.gz", repos = NULL)

Then, we need to install `survival`:  

    install.packages("survival")
    
which should install within a couple of minutes on a standard machine.
   
To run the simulated data example, you will also need to install `MCMCpack`:  

    install.packages("MCMCpack")

# Demo

We first load all the source code dependencies:

```
library(camp)
library(survival)
library(MCMCpack)
```

The data can be generated by running the following code:

```
seed <- 1
n1 <- 100
n2 <- 100
p <- 1000
signal_prop <- 0.8
signal <- p * signal_prop
null <- p - signal
sum_const <- signal_prop
a <- 353

#Function for generating OTU data
GenData <- function(seed, n1, n2, p, 
                    signal_prop, sum_const, a,
                    L1 = 25000, U1 = 40000,
                    L2 = 25000, U2 = 40000) {
  set.seed(seed)
  N1 <- sample(L1:U1, n1, replace = T) #Raw library size group 1
  N2 <- sample(L2:U2, n2, replace = T) #Raw library size group 2
  
  null_prop <- 1 - signal_prop
  
  signal = p * signal_prop
  null <- p - signal
  
  #Generate True relative abundance vector for group 1
  P1_signal <- runif(signal)
  P1_signal <- sum_const * P1_signal / sum(P1_signal)
  P1_null <- runif(null)
  P1_null <- (1 - sum_const) * P1_null / sum(P1_null)
  P1 <- c(P1_signal, P1_null)
  
  P2_signal <- runif(signal)
  P2_signal <- sum_const * P2_signal / sum(P2_signal)
  P2_null <- P1_null
  P2 <- c(P2_signal, P2_null)
  
  a1 <- rdirichlet(n1, a*P1)
  
  #Generate OTU matrix for group 1
  X1 <- matrix(0, n1, p)
  for (i in 1:n1) {
    X1[i,] <- rmultinom(n = 1, size = N1[i], prob = a1[i,])
  }
  
  a2 <- rdirichlet(n2, a*P2)
  
  X2 <- matrix(0, n2, p)
  for (i in 1:n2) {
    X2[i,] <- rmultinom(n = 1, size = N2[i], prob = a2[i,])
  }
  
  X <- rbind(X1, X2)
  
  cond <- c(rep("Group 1", n1), rep("Group 2", n2))
  return(list(OTU = X, Group = cond))
}

camp_dat <- GenData(seed, n1, n2, p, signal_prop, sum_const, a)
```

camp_dat is a list which contains a (n x p) OTU table and a length n vector of group labels.

Then, in order to run camp, we simply need to input the OTU table and corresponding group labels.
```
camp_res <- camp(camp_dat$OTU, camp_dat$Group)
head(camp_res)
```

# Output
```
> head(camp_res)
[1] 4.381610e-03 1.110223e-16 1.297482e-08 3.955734e-03 1.526363e-07
[6] 1.177519e-05

```
