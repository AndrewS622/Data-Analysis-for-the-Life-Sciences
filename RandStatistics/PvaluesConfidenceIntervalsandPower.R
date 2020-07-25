# t-test Exercises Using babies.txt dataset
library(downloader)
url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/babies.txt"
filename <- basename(url)
download(url, destfile=filename)
babies <- read.table("babies.txt", header=TRUE)

# Differences in birthweight for smoking and non-smoking mothers
library(dplyr)
bwt.nonsmoke <- filter(babies, smoke==0) %>% select(bwt) %>% unlist()
bwt.smoke <- filter(babies, smoke==1) %>% select(bwt) %>% unlist()

# Assume large dataset represents entire population
library(rafalib)
mean(bwt.nonsmoke)-mean(bwt.smoke)
popsd(bwt.nonsmoke)
popsd(bwt.smoke)

# Take random samples and find t statistic
set.seed(1)
N <- 25
dat.ns <- sample(bwt.nonsmoke, N)
dat.s <- sample(bwt.smoke, N)
tval <- abs(mean(dat.ns) - mean(dat.s)) / sqrt(var(dat.ns)/N + var(dat.s)/N)
tval

# Determine significance
pval <- 1-(pnorm(abs(tval))-pnorm(-abs(tval)))
pval

pval <- 2*pnorm(-abs(tval))
pval

## Confidence Intervals
# continued from above

Qn <- qnorm(1-0.01/2)
Qt <- qt(1-0.01/2, df = 2*N-2)
se <- sqrt(var(dat.ns)/N + var(dat.s)/N)
diff <- abs(mean(dat.ns)-mean(dat.s))
tval <- diff/se

# For significance level alpha, prob of a Type I error
# i.e. failing to reject null, is alpha

interval <- c(diff-Qt*se, diff+Qt*se)
interval

interval <- c(diff-Qn*se, diff+Qn*se)
interval

t.test(dat.ns, dat.s)

set.seed(1)
N <- 5
t.test(sample(bwt.smoke, N), sample(bwt.nonsmoke, N))

## Power Calculations
# Power = 1 - prob of type II error, 
# TII error is failing to accept null when it is true

N <- 5
set.seed(1)
Ns <- 10000
rejects <- replicate(Ns, t.test(sample(bwt.smoke, N), sample(bwt.nonsmoke, N))$p.value)
mean(abs(rejects) < 0.05)
# Power is very low for N = 5
# Repeat for 30, 60

N <- 30
set.seed(1)
Ns <- 10000
rejects <- replicate(Ns, t.test(sample(bwt.smoke, N), sample(bwt.nonsmoke, N))$p.value)
mean(abs(rejects) < 0.05)

N <- 60
set.seed(1)
Ns <- 10000
rejects <- replicate(Ns, t.test(sample(bwt.smoke, N), sample(bwt.nonsmoke, N))$p.value)
mean(abs(rejects) < 0.05)

# Decreasing alpha requires more samples for same power
N <- 90
set.seed(1)
Ns <- 10000
rejects <- replicate(Ns, t.test(sample(bwt.smoke, N), sample(bwt.nonsmoke, N))$p.value)
mean(abs(rejects) < 0.01)

