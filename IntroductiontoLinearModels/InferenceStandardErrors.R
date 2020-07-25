## Inference Review
g <- 9.8
h0 <- 56.67
v0 <- 0
n <- 25
tt <- seq(0,3.4,len=n)
y <- h0 + v0 * tt - 0.5 * g * tt^2 + rnorm(n,sd=1)
# y are observed values, need to estimate parameters

X <- cbind(1,tt,tt^2)
A <- solve(crossprod(X)) %*% t(X)
-2 * (A %*% y) [3] #estimate for g

# Estimating standard error
set.seed(1)
B <- 100000
fallobj <- function() {
  y <- h0 + v0 * tt - 0.5 * g * tt^2 + rnorm(n,sd=1)
  X <- cbind(1,tt,tt^2)
  A <- solve(crossprod(X)) %*% t(X)
  -2 * (A %*% y) [3]
}
gest <- replicate(B, fallobj())
sd(gest)

# Father-son height dataset
library(UsingR)
x <- father.son$fheight
y <- father.son$sheight
n <- length(y)

library(dplyr)
N <- 50
B <- 10000
set.seed(1)
fatherson <- function() {
  index <- sample(n,N)
  sampledat <- father.son[index,]
  x <- sampledat$fheight
  y <- sampledat$sheight
  betahat <- lm(y~x)$coef
}
fsEst <- replicate(B, fatherson())
sd(fsEst[2,])

# Covariance
mean((y - mean(y)) * (x - mean(x)))

## Standard Errors
set.seed(1)
index <- sample(n, N)
sampledat <- father.son[index,]
x <- sampledat$fheight
y <- sampledat$sheight
fit <- lm(y~x)
betahat <- fit$coef

# variance is sigma^2 (Xt X)^-1

r <- y - fit$fitted.values
SSR <- sum(r^2)
sigma2 <- SSR/(50 - 2)

X <- cbind(1, x)
inv <- solve(t(X) %*% X)
inv[1,1]
varEst <- sigma2 * diag(inv)
sqrt(varEst)