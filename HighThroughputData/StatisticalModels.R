## Statistical Models
# Binomial distribution
pgirl <- 0.49
dbinom(2,4,pgirl)
dbinom(4,10,pgirl)

pgCG <- 0.40
1 - pbinom(10,20,pgCG)

plotto <- 1/175223510
Ntickets <- 189000000
1 - pbinom(0,N,plotto)
1 - pbinom(1,N,plotto)

# Normal approximation
N <- 20
binomEst <- pbinom(0.45*N, N, pgCG) - pbinom(0.35*N, N, pgCG)
z45 <- (0.45 * N - pgCG*N)/sqrt(pgCG*N*(1-pgCG))
z35 <- (0.35 * N - pgCG*N)/sqrt(pgc*N*(1-pgCG))
pnorm(0.45*N,pgc*N,sqrt(pgc*N*(1-pgc))) - pnorm(0.35*N,pgc*N,sqrt(pgc*N*(1-pgc)))
normEst <- pnorm(z45)- pnorm(z35)
abs(normEst - binomEst)

N <- 1000
binomEst <- pbinom(0.45*N, N, pgCG) - pbinom(0.35*N, N, pgCG)
z45 <- (0.45 * N - pgCG*N)/sqrt(pgCG*N*(1-pgCG))
z35 <- (0.35 * N - pgCG*N)/sqrt(pgc*N*(1-pgCG))
pnorm(0.45*N,pgc*N,sqrt(pgc*N*(1-pgc))) - pnorm(0.35*N,pgc*N,sqrt(pgc*N*(1-pgc)))
normEst <- pnorm(z45)- pnorm(z35)
abs(normEst - binomEst)

Ns <- c(5,10,30,100)
ps <- c(0.01,0.10,0.5,0.9,0.99)
library(rafalib)
mypar(length(Ns),length(ps))
pMe <- function(i) {
  N <- Ns[ceiling(i/5)]
  p <- ps[i %% 5]
  if (identical(p, numeric(0))) p = ps[5]
  ks <- 1:(N-1)
  pexact <- mapply(dbinom, ks, size=N, prob=p)
  a <- (ks + 0.5 - N*p)/sqrt(N*p*(1-p))
  b <- (ks - 0.5 - N*p)/sqrt(N*p*(1-p))
  papprox <- pnorm(a) - pnorm(b)
  plot(pexact,papprox,xlab=paste("N = ", N, " and p = ", p))
  c(pexact, papprox)
}
is <- 1:(length(Ns)*length(ps))
pMes <- sapply(is, pMe)
# Approximations break down when p very small or large

mypar(1,1)

N <- 189000000
p <- 1/175223510
dbinom(2,N,p)
a <- (2.5 - N*p)/sqrt(N*p*(1-p))
b <- (1.5 - N*p)/sqrt(N*p*(1-p))
pnorm(a) - pnorm(b)
100*(dbinom(2,N,p) - pnorm(a) + pnorm(b))/dbinom(2,N,p)

lambda <- N*p
dpois(2,lambda) - dbinom(2,N,p)
1 - ppois(1, lambda)

## Maxmimum Likelihood Exercises
library(devtools)
install_github("genomicsclass/dagdata")
# palindromes in human cytomegalovirus genome
library(dagdata)
data(hcmv)
library(rafalib)
mypar()
plot(locations,rep(1,length(locations)),ylab="",yaxt="n")

# break genome into 4000-bp sequences so Np not so small
breaks <- seq(0, 4000*round(max(locations)/4000),4000)
tmp <- cut(locations,breaks)
counts=as.numeric(table(tmp))
hist(counts)

probs <- dpois(counts,4)
likelihood <- prod(probs)
likelihood

# use log-likelihood instead since likelihood is so small
logprobs <- dpois(counts,4,log=TRUE)
loglikelihood <- sum(logprobs)
loglikelihood

# optimize lambda
lambdas <- seq(0,15,len=300)
MLE <- function(lambda) {
  logprobs <- dpois(counts, lambda, log=TRUE)
  loglikelihood <- sum(logprobs)
}
LEs <- sapply(lambdas, MLE)
hist(LEs)
lambdas[which(LEs == max(LEs))]

# counts per location in genome
binLocation <- (breaks[-1]+breaks[-length(breaks)]) / 2
plot(binLocation,counts,type="l",xlab=)
binLocation[which(counts == max(counts))]
max(counts)

# is this bigger than expected by chance?
lambda <- mean(counts[-which.max(counts)])
pval <- 1 - ppois(14,lambda)
pval
0.05 / length(counts) #Bonferroni

# check QQ plots
ps <- (seq(along=counts) - 0.5)/length(counts)
lambda <- mean(counts[-which.max(counts)])
poisq <- qpois(ps,lambda)
qqplot(poisq,counts)
abline(0,1)
qqnorm(counts)

## Models for Variance
library(devtools)
install_github("genomicsclass/tissuesGeneExpression")
library(tissuesGeneExpression)
data("tissuesGeneExpression")
library(genefilter)
y <- e[,which(tissue=="endometrium")]
vars <- rowVars(y)
# clearly not normal
qqnorm(vars)
qqnorm(sqrt(vars))

library(limma)
Fdist <- fitFDist(vars, 14)
Fdist

ps <- seq(0,0.99,len=length(vars))
qs <- mapply(qf, ps, df1=14, df2=Fdist$df2)
qqplot(sqrt(qs),sqrt(vars))
qqline(sqrt(vars))
# decent approximation except for ~top 5% of variances