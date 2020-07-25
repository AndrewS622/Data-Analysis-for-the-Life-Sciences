## Monte Carlo Simulation
set.seed(1)
N <- 5
x <- rnorm(N)
tval <- sqrt(N)*mean(x)/sd((x))
tval

# simulating the t distribution
set.seed(1)
B <- 1000
tstat <- function(N) {
  x <- rnorm(N)
  tval <- sqrt(N)*mean(x)/sqrt(var(x))
  return(tval)
}
tstats <- replicate(B, tstat(N))
mean(tstats > 2)
1-pt(2,df=4)

# QQ plots for the t distribution
B <- 100
ps <- seq(1/(B+1), 1-1/(B+1), len=B)
Ns <- seq(5,30,5)
qs <- qt(ps, df=4)

set.seed(1)
B <- 1000
tstat <- function(N) {
  x <- rnorm(N)
  tval <- sqrt(N)*mean(x)/sqrt(var(x))
  return(tval)
}
tstats <- sapply(Ns, function(n) {
  replicate(B, tstat(n))
})
pvals <- sapply(ps, function(p) {
  colMeans(tstats < p)
})

library(rafalib)
mypar(3,2)
for (i in 1:6) {
  qqplot(qs, pvals[i,])
  qqline(pvals[i,])
}

# Repeat using t.test function
set.seed(1)
tstats <- function(N) {
  stat <- t.test(rt(N, df=N-1), rt(N, df=N-1), var.equal = TRUE)$statistic
  return(stat)
}
ts <- sapply(Ns, function(n) {
  replicate(B, tstats(n))
})
pvals <- sapply(ps, function(p) {
  colMeans(ts < p)
})

mypar(3,2)
for (i in 1:6) {
  qqplot(qs, pvals[i,])
  qqline(pvals[i,])
}

# Example: Binary data
par(mfrow = c(1,1))
set.seed(1)
N <- 15
B <- 10000
tstats <- replicate(B,{
  X <- sample(c(-1,1), N, replace=TRUE)
  sqrt(N)*mean(X)/sd(X)
})
ps=seq(1/(B+1), 1-1/(B+1), len=B) 
qqplot(qt(ps,N-1), tstats, xlim=range(tstats))
abline(0,1)
# Large tails confirm non-normality, cannot use t

# However, if we make N = 1000
set.seed(1)
N <- 1000
B <- 10000
tstats <- replicate(B,{
  X <- sample(c(-1,1), N, replace=TRUE)
  sqrt(N)*mean(X)/sd(X)
})
ps=seq(1/(B+1), 1-1/(B+1), len=B) 
qqplot(qt(ps,N-1), tstats, xlim=range(tstats))
abline(0,1)
# Because of CLT, average is ~normal and df=999 ~ normal

# Distributions of the median
B <- 10000
N <- 5
zs <- function(N) {
  X <- rnorm(N)
  return(median(X))
}
zsamples <- replicate(B, zs(N))
qqnorm(zsamples)
abline(0,1/sqrt(N))

B <- 10000
N <- 50
zs <- function(N) {
  X <- rnorm(N)
  return(median(X))
}
zsamples <- replicate(B, zs(N))
qqnorm(zsamples)
abline(0,1/sqrt(N))
# Approximately normal with mean 0 and sd > 1/sqrt(N)

## Permutations
library(downloader)
url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/babies.txt"
filename <- basename(url)
download(url, destfile=filename)
babies <- read.table("babies.txt", header=TRUE)

library(dplyr)
bwt.nonsmoke <- filter(babies, smoke==0) %>% select(bwt) %>% unlist 
bwt.smoke <- filter(babies, smoke==1) %>% select(bwt) %>% unlist

N <- 10

# Generate null distribution from permutations
set.seed(1)
B <- 1000
permute <- function() {
  dat <- c(smokers,nonsmokers)
  shuffle <- sample(dat)
  smokersstar <- shuffle[1:N]
  nonsmokersstar <- shuffle[(N+1):(2*N)]
  return(mean(smokersstar)-mean(nonsmokersstar))
}
permutations <- replicate(B, permute())
# Note: p-value for permutation test adds 1 in numerator and denominator
(sum(abs(permutations) > abs(obs)) + 1) / (length(permutations) + 1)
hist(permutations)
abline(v=obs)     

# Now for the median
set.seed(1)

permute <- function() {
  dat <- c(smokers,nonsmokers)
  shuffle <- sample(dat)
  smokersstar <- shuffle[1:N]
  nonsmokersstar <- shuffle[(N+1):(2*N)]
  return(median(smokersstar)-median(nonsmokersstar))
}
permutations <- replicate(B, permute())
(sum(abs(permutations) >= abs(obs)) + 1) / (length(permutations) + 1)
hist(permutations)
abline(v=obs) 

## Association Tests
library(dplyr)
d <- read.csv("assoctest.csv")
tab <- table(d)
chisq.test(tab)
fisher.test(tab)
# Fisher's exact test uses hypergeometric distn 
# and requires knowing the number of observations
# Chi-squared uses F distribution and does not require this