## Contrasts
species <- factor(c("A", "A", "B", "B"))
condition <- factor(c("control", "treated", "control", "treated"))
model.matrix(~ species + condition)
library(contrast)

# Using leg-pair/frictional coefficient dataset
url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/spider_wolff_gorb_2013.csv"
filename <- "spider_wolff_gorb_2013.csv"
library(downloader)
if (!file.exists(filename)) download(url, filename)
spider <- read.csv(filename, skip=1)
friction <- spider$friction
type <- spider$type
leg <- spider$leg

fit <- lm(friction ~ type + leg)
contrast(fit, list(leg = "L4", type = "push"), list(leg = "L2", type = "push"))

# SE of contrast C*Beta-hat = sqrt(CSigmaCt)
# Sigma is the covariance matrix of Beta-hat
# diag(Sigma) gives variance of each coefficient
X <- model.matrix(~ type + leg, data=spider)
(Sigma <- sum(fit$residuals^2)/(nrow(X) - ncol(X)) * solve(t(X) %*% X))
C <- matrix(c(0,0,-1,0,1),1,5)
Sigma
Sigma["legL4", "legL2"]

# Var(betaL4 - betaL2) = var(betaL4) + var(betaL2) - 2Cov(betaL4, betaL2)
sqrt(Sigma["legL4", "legL4"] + Sigma["legL2", "legL2"] - 2*Sigma["legL4", "legL2"])
sqrt(C %*% Sigma %*% t(C))

## Interactions
url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/spider_wolff_gorb_2013.csv"
filename <- "spider_wolff_gorb_2013.csv"
library(downloader)
if (!file.exists(filename)) download(url, filename)
spider <- read.csv(filename, skip=1)

# Apply log2 to make within-group variances comparable
spider$log2friction <- log2(spider$friction)

boxplot(friction ~ type*leg, data=spider)
boxplot(log2friction ~ type*leg, data=spider)

fit <- lm(log2friction ~ type*leg, data=spider)
summary(fit)
summary(fit)$coefficients['typepush:legL4', 't value']

# Are all leg pairs affected by push vs. pull the same?
anova(fit)
anova(fit)["type:leg", "F value"]

library(contrast)
contrast(fit, list(type = "pull", leg = "L1"), list(type = "pull", leg = "L2"))
contrast(fit, list(type = "push", leg = "L1"), list(type = "push", leg = "L2"))

## F-test
# Explanatory power of model
# mean SS explained / mean SSR
# Generate null F distribution
N <- 40
p <- 4
B <- 1000
group <- factor(rep(1:p,each=N/p))
X <- model.matrix(~ group)

fdistn <- function() {
  Y <- rnorm(N,mean=42,sd=7)
  mu0 <- mean(Y)
  initial.ss <- sum((Y-mu0)^2)
  
  s <- split(Y, group)
  after.group.ss <- sum(sapply(s, function(x) sum((x-mean(x))^2)))
  (group.ss <- initial.ss - after.group.ss)
  
  group.ms <- group.ss / (p-1)
  after.group.ms <- after.group.ss / (N-p)
  f.value <- group.ms / after.group.ms
}

set.seed(1)
fval <- replicate(B, fdistn())
mean(fval)

hist(fval, col="grey", border="white", breaks=50, freq=FALSE)
xs <- seq(from=0,to=6,length=100)
lines(xs, df(xs, df1 = p - 1, df2 = N - p), col="red")
