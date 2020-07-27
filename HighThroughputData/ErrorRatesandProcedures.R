## Exercises - Error Rates
N <- 8793
p <- 0.05
1-p^N
# prob of incorrectly rejecting at least one null = 1
1-(1-.05)^(1/N)
# using Sidak's procedure - necessary cutoff for 5% error rate
# but requires statistical independence

## Bonferroni Correction (alpha/m)
alphas <- seq(0,0.25,0.01)
m <- seq(2, 10,0.05)
length(alphas)
library(rafalib)
mypar(5,6)
sapply(alphas, function(alpha){
  plot(alpha/m, (1-(1-alpha)^(1-m)))
})
# Rejects less null hypotheses than Sidak

# MC simulation of p-values w/ bonferroni cutoff
mypar(1)
alpha <- 0.05
B <- 10000
cutoff <- alpha/N
set.seed(1)
pvalSample <- function() {
  pvals <- runif(N,0,1)
  sum(pvals < cutoff)
}
pvalsEst <- replicate(B, pvalSample())
mean(pvalsEst)

# For Sidak's procedure
set.seed(1)
pvalSample <- function() {
  pvals <- runif(N,0,1)
  sum(pvals < 1-(1-alpha)^(1/N))
}
pvalsEst <- replicate(B, pvalSample())
mean(pvalsEst)

## False Discovery Rate
library(devtools)
library(rafalib)
install_github("genomicsclass/GSE5859Subset")
library(BiocManager)
BiocManager::install(c("genefilter", "qvalue"))

library(GSE5859Subset)
data(GSE5859Subset)
library(genefilter)

groups <- factor(sampleInfo$group)
pvals <- rowttests(geneExpression, groups)$p.value
sum(pvals < 0.05)
sum(pvals < 0.05/N)
# FWER < 0.05 using Bonferroni

# FDR is defined on a list of features
# mth q-value is FDR of list of m most significant tests ordered by p-value

pvalsadjust <- p.adjust(pvals, method = "fdr")
sum(pvalsadjust < 0.05)

library(qvalue)
qvals <- qvalue(pvals)$qvalues
sum(qvals < 0.05)
pi0 <- qvalue(pvals)$pi0
pi0
# pi0 is estimate for proportion of features where null is true
# qvalue function provides less conservative estimate

ps <- seq(0, 0.5, 0.005)
pqratio <- sapply(ps, function(p) {
  pv <- sum(pvalsadjust < p)
  qv <- sum(qvals < p)
  pv/qv
})
plot(ps, pqratio)

# MC simulation where cases have +2 expression
n <- 24
m <- 8793
delta <- 2
positives <- 500
m0 <- m - positives
m1 <- positives

B <- 1000

set.seed(1)
groups <- factor(c(rep(1,n/2), rep(0,n/2)))

# compute significant p-values from each procedure
corrections <- function(mat) {
  mat <- matrix(rnorm(n*m),m,n)
  mat[1:positives,1:(n/2)] <- mat[1:positives,1:(n/2)] + delta
  
  pvals <- rowttests(mat, groups)$p.value
  pvalsadjust <- p.adjust(pvals, method = "fdr")
  qvals <- qvalue(pvals)$qvalues
  
  matrix(c(pvals < 0.05/m, pvalsadjust < 0.05, qvals < 0.05), nrow = length(pvals))
}
simData <- replicate(B, corrections(mat))

# true values for significant genes
trueVals <- c(rep(TRUE, m1), rep(FALSE, m0))
comp <- !xor(simData, trueVals)

# False positive and false negative rates
FPs <- 1 - rowMeans(colSums(comp[501:m,,]))/m0
FPs

FNs <- 1 - rowMeans(colSums(comp[1:500,,]))/m1
FNs