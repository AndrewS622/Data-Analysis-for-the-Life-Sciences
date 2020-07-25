## The Normal Distribution
library(downloader)
url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/femaleControlsPopulation.csv"
filename <- basename(url)
download(url, destfile=filename)
x <- unlist(read.csv(filename))

set.seed(1)
xm <- mean(x)
xsm5 <- vector("numeric", 1000)
for (i in 1:1000) {
  xs <- sample(x, 5)
  xsm5[i] <- mean(xs)
}

set.seed(1)
xsm50 <-vector("numeric", 1000)
for (i in 1:1000) {
  xs <- sample(x, 50)
  xsm50[i] <- mean(xs)
}

par(mfrow=c(1,2))
hist(xsm5)
hist(xsm50)
(sum(xsm50 < 25) - sum(xsm50 < 23))/1000

pnorm(25, 23.9, 0.43) - pnorm(23, 23.9, 0.43)

## Populations, Samples, and Estimates
url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/mice_pheno.csv"
filename <- basename(url)
download(url, destfile=filename)
dat <- read.csv(filename) 
dat <- na.omit(dat)

library(dplyr)
x <- filter(dat, Diet == "chow") %>% 
  filter(Sex == "M") %>% 
  select(Bodyweight) %>% 
  unlist()
mean(x)

library(rafalib)
popsd(x)
N <- length(x)
sd(x)*sqrt((N-1)/N) #true population standard deviation

set.seed(1)
xs <- sample(x, 25)
mean(xs)

y <- filter(dat, Diet == "hf") %>% 
  filter(Sex == "M") %>% 
  select(Bodyweight) %>% 
  unlist()
mean(y)
popsd(y)

set.seed(1)
ys <- sample(y, 25)
mean(ys)

abs(abs(mean(y)-mean(x))
    -abs(mean(ys)-mean(xs)))


w <- filter(dat, Diet == "chow") %>% 
  filter(Sex == "F") %>% 
  select(Bodyweight) %>% 
  unlist()
z <- filter(dat, Diet == "hf") %>% 
  filter(Sex == "F") %>% 
  select(Bodyweight) %>% 
  unlist()
popsd(w)
popsd(z)
set.seed(2)
ws <- sample(w, 25)    
set.seed(2)
zs <- sample(z, 25)
abs(abs(mean(w)-mean(z))-abs(mean(ws)-mean(zs)))
# sample less variable for women 
# possibly because population variance 
# smaller for females

## The Central Limit Theorem

# Properties of the normal
pnorm(1) - pnorm(-1)
pnorm(2) - pnorm(-2)
pnorm(3) - pnorm(-3)

library(downloader) 
url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/mice_pheno.csv"
filename <- basename(url)
download(url, destfile=filename)
dat <- na.omit( read.csv(filename) )

library(dplyr)
y <- filter(dat, Diet == "chow") %>% filter(Sex == "M") %>% select(Bodyweight) %>% unlist()

# Properties of the population
ym <- mean(y)
library(rafalib)
ysd <- popsd(y)

# Compare to theoretical
mean(abs(y-ym) < ysd)
mean(abs(y-ym) < 2 * ysd)
mean(abs(y-ym) < 3 * ysd)

# Confirm normality with QQ plot
qqnorm((y - ym)/ysd)
abline(0,1)

# Check for all four populations of sex and diet
par(mfrow = c(2,2))
y <- filter(dat, Sex=="M" & Diet=="chow") %>% select(Bodyweight) %>% unlist
z <- ( y - mean(y) ) / popsd(y)
qqnorm(z);abline(0,1)
y <- filter(dat, Sex=="F" & Diet=="chow") %>% select(Bodyweight) %>% unlist
z <- ( y - mean(y) ) / popsd(y)
qqnorm(z);abline(0,1)
y <- filter(dat, Sex=="M" & Diet=="hf") %>% select(Bodyweight) %>% unlist
z <- ( y - mean(y) ) / popsd(y)
qqnorm(z);abline(0,1)
y <- filter(dat, Sex=="F" & Diet=="hf") %>% select(Bodyweight) %>% unlist
z <- ( y - mean(y) ) / popsd(y)
qqnorm(z);abline(0,1)

# Distributions of averages
y <- filter(dat, Sex=="M" & Diet=="chow") %>% select(Bodyweight) %>% unlist
set.seed(1)
avgs <- replicate(10000, mean( sample(y, 25)))
mypar(1,2)
hist(avgs)
qqnorm(avgs)
qqline(avgs)
mean(avgs)
popsd(avgs)

## The Central Limit Theorem and the t Distribution
library(downloader)
url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/femaleMiceWeights.csv"
filename <- "femaleMiceWeights.csv"
if(!file.exists("femaleMiceWeights.csv")) download(url,destfile=filename)
dat <- read.csv(filename)

# For simple die with known probability:
# mu = p = 1/6
# sigma = sqrt(p*(1-p)/n)
set.seed(1)
n <- 100
ntimes <- 10000
x <- replicate(ntimes, sample(1:6, n, replace=TRUE))
p <- 1/6
z <- vector("numeric", ntimes)
for (i in 1:ntimes) {
  z[i] <- (mean(x[,i]==6) - p) / sqrt((p * (1-p) / n))
}
mean(abs(z) > 2)

par(mfrow = c(1,1))
qqnorm(z)
abline(0,1)

# Dependence of normality on p, n
par(mfrow = c(2,2))
ps <- c(0.5, 0.5, 0.01, 0.01)
ns <- c(5, 30, 30, 100)
is <- 1:4
zs <- sapply(is, function(i) {
  set.seed(1)
  n <- ns[i]
  ntimes <- 10000
  p <- ps[i]
  x <- replicate(ntimes, sample(1:1/p, n, replace=TRUE))
  z <- vector("numeric", ntimes)
  for (j in 1:ntimes) {
    z[j] <- (mean(x[,j]==1) - p) / sqrt((p * (1-p) / n))
  }
  mean(abs(z) > 2)
  qqnorm(z)
  qqline(z)
})

#Need to estimate parameters from non-binary distributions
library(dplyr)
X <- filter(dat, Diet=="chow") %>% 
  select(Bodyweight) %>% 
  unlist()
Y <- filter(dat, Diet=="hf") %>% 
  select(Bodyweight) %>% 
  unlist()
ux <- mean(X)
uy <- mean(Y)
sdx <- sd(X)
sdy <- sd(Y)

# Prob that ux is more than 2 off from true mean
2 * (1 - pnorm(sqrt(12) * 2 / (sdx)))

# Standard error
se <- sqrt(sdx^2/12 + sdy^2/12)

# t statistic
tstat <- (uy - ux) / se

# Properties of the t distribution
1 - pt(3, df=3)
1 - pt(3, df=15)
1 - pt(3, df=30)
1 - pnorm(3)

# Significance testing
pval <- 2 * (1 - pnorm(abs(tstat))) # equivalent to
pval <- 1 - pnorm(abs(tstat)) + pnorm(-abs(tstat))

# Check for normality of sample
par(mfrow = c(1,2))
qqnorm(X)
qqline(X)
qqnorm(Y)
qqline(Y)

# Use t-test
t.test(Y, X)

# Larger p-value due to accounting for 
# variability in estimating standard error