#' Generate simulated micrboiome data
#'
#' @param seed An integer, seed used for reproducibility in random number generation
#' @param n1 An integer, the number of samples in Group 1
#' @param n2 An integer, the number of samples in Group 2
#' @param p An integer, the number of taxa to be generated
#' @param signal.prop A number between 0 to 1, indicating the proportion of taxa that are signals (non-null)
#' @param sum.const A number between 0 to 1, representing the signal taxa's cumulative relative abundance (i.e., how much should they add up to)
#' @param a A positive number for the Dirichlet distribution (controls the dispersion of the abundance)
#' @param L1 An integer, the minimum library size sampling in Group 1. Default is 25000.
#' @param U1 An integer, the maximum library size sampling in Group 1. Default is 40000.
#' @param L2 An integer, the minimum library size sampling in Group 2. Default is 25000.
#' @param U2 An integer, the maximum library size sampling in Group 2. Default is 40000.
#'
#' @return A list containing two elements:
#' \item{OTU}{A size ((n1 + n2) x p) matrix, representing the OTU table for both groups combined}
#' \item{Group}{A length (n1 + n2) vector, indicating the group labels ("Group 1" or "Group 2") for each sample}
#' @examples
#' \dontrun{
#' seed <- 1
#' n1 <- 100
#' n2 <- 100
#' p <- 1000
#' signal.prop <- 0.8
#' sum.const <- signal.prop
#' a <- 353
#'
#' camp.dat <- GenData(seed, n1, n2, p, signal.prop, sum.const, a)
#' head(camp.dat)
#' }
#' @export
GenData <- function(seed, n1, n2, p,
                    signal.prop, sum.const, a,
                    L1 = 25000, U1 = 40000,
                    L2 = 25000, U2 = 40000) {
  set.seed(seed)
  N1 <- sample(L1:U1, n1, replace = T) #Raw library size group 1
  N2 <- sample(L2:U2, n2, replace = T) #Raw library size group 2

  null_prop <- 1 - signal.prop

  signal = p * signal.prop
  null <- p - signal

  #Generate True relative abundance vector for group 1
  P1.signal <- runif(signal)
  P1.signal <- sum.const * P1.signal / sum(P1.signal)
  P1.null <- runif(null)
  P1.null <- (1 - sum.const) * P1.null / sum(P1.null)
  P1 <- c(P1.signal, P1.null)

  P2.signal <- runif(signal)
  P2.signal <- sum.const * P2.signal / sum(P2.signal)
  P2.null <- P1.null
  P2 <- c(P2.signal, P2.null)

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

#' Censoring-based analysis of microbiome proportions (CAMP)
#'
#' @param  X A size (n x p) matrix, p is the number of taxa
#' @param  cov A length n vector (in case of only group label) or a size (n x q) data frame, q is the number of covariates (including group)
#' @return A length p vector of p-values, or a size (p x q) data frame of p-values corresponding the covariates
#' @examples
#' \dontrun{
#' seed <- 1
#' n1 <- 100
#' n2 <- 100
#' p <- 1000
#' signal.prop <- 0.8
#' sum.const <- signal.prop
#' a <- 353
#'
#' camp.dat <- GenData(seed, n1, n2, p, signal.prop, sum.const, a)
#'
#' camp.res <- camp(camp.dat$OTU, camp.dat$Group)
#' head(camp.res)
#' }
#' @export
camp <- function(X, cov) {
  n <- dim(X)[1] # Sample size
  p <- dim(X)[2] # Number of taxa
  Delta <- (X != 0) + 0 # Creating censoring indicator for survival analysis
  T <- X
  T[which(T == 0, arr.ind = TRUE)] <- min(X[X != 0]) # Zero replacement

  # Normalize the data
  T <- T / matrix(rowSums(T), n, p)
  T <- -log(T) # "Time" data

  if (is.vector(cov)) {
    pval.vec <- rep(NA, p)
    for (i in 1:p) {
      fit <- survdiff(Surv(T[,i], Delta[,i]) ~ cov)
      pval.vec[i] <- 1 - pchisq(fit$chisq, 1)
    }
    return(pval.vec)
  } else if (is.data.frame(cov)) {
    q <- dim(cov)[2]
    pval.df <- data.frame(matrix(NA, p, q))
    colnames(pval.df) <- names(cov)
    for (i in 1:p) {
      formula_str <- paste0("Surv(T[,", i, "], Delta[,", i, "]) ~ ", paste(names(cov), collapse = " + "))
      fit <- coxph(as.formula(formula_str))
      pval.df[i, ] <- summary(fit)$coefficients[1:q, 5]
    }
    return(pval.df)
  } else {
    print("Cov is of the wrong data type!")
  }
}
