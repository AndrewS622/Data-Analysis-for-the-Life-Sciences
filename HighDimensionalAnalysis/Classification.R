## Conditional Expectation
# for binary data, probability = expectation
n <- 1000
y <- rbinom(n,1,0.25)
sum(y==1)/length(y)
mean(y)

# random data simulating heights for men and women
n <- 10000
set.seed(1)
men <- rnorm(n,176,7)
women <- rnorm(n,162,7)
y <- c(rep(0,n),rep(1,n))
x <- round(c(men,women))

# mix up data
ind <- sample(seq(along=y))
y <- y[ind]
x <- x[ind]

# probability of being a woman given height = 176 cm
library(tidyverse)
df <- data.frame(height = x, woman = y)
women176 <- df %>% filter(height == 176 & woman == 1)
people176 <- df %>% filter(height == 176)
dim(women176)[1]/dim(people176)[1]

# probability of being a woman vs. height
womanProb <- function(i) {
  tallWomen <- df %>% filter(height == i & woman == 1)
  tall <- df %>% filter(height == i)
  sum(tallWomen)/sum(tall)
}

hs <- seq(160,178)
w <- sapply(hs, womanProb)

plot(hs, w)

# max height where probability > 0.5
hs[which(w > 0.5)[length(which(w > 0.5))]]

## Smoothing and loess
RNGkind(sample.kind = "Rounding")

n = 10000
set.seed(1)
men = rnorm(n,176,7) #height in centimeters
women = rnorm(n,162,7) #height in centimeters
y = c(rep(0,n),rep(1,n))
x = round(c(men,women))
ind = sample(seq(along=y))
y = y[ind]
x = x[ind]

# sample 250 from population
set.seed(5)
N <- 250
ind <- sample(length(y), N)
Y <- y[ind]
X <- x[ind]

# fit local regression
l <- loess(Y~X)
predict(l, data.frame(X=168))

# find standard error of this prediction
set.seed(5)
B <- 1000
ls <- replicate(B, {
  ind <- sample(length(y), N)
  Y <- y[ind]
  X <- x[ind]
  
  l <- loess(Y~X)
  predict(l, data.frame(X=168))
})
sd(ls)

## kNN and Cross-Validation
RNGkind(sample.kind = "Rounding")
library(GSE5859Subset)
data(GSE5859Subset)
y <- factor(sampleInfo$group)
X <- t(geneExpression)

# only consider autosomal genes
out <- which(geneAnnotation$CHR %in% c("chrX", "chrY"))
X <- X[,-out]
library(caret)
library(class)

# create folds for 10X cross-validation
set.seed(1)
idx <- createFolds(y, 10)
idx$Fold03[2]

# find 8 genes with smallest p-values
# use fold 2 as the validation set
m <- 8
k <- 5
library(genefilter)
train <- X[-idx$Fold02,]
ytrain <- y[-idx$Fold02]
validate <- X[idx$Fold02,]
yvalidate <- y[idx$Fold02]

pvals <- colttests(train, ytrain)$p.value
pvalsrank <- rank(pvals)
smallp <- which(pvalsrank %in% 1:m)

# train k-nearest neighbors model with 5 centers
pred <- knn(train=train[,smallp], test = validate[,smallp], cl = ytrain, k = k)
sum(pred != y[idx$Fold02])

# repeat for all 10 folds and determine error rate
folds <- 1:10
predictions <- sapply(folds, function(fold) {
  train <- X[-idx[[fold]],]
  ytrain <- y[-idx[[fold]]]
  validate <- X[idx[[fold]],]
  yvalidate <- y[idx[[fold]]]
  
  pvals <- colttests(train, ytrain)$p.value
  pvalsrank <- rank(pvals)
  smallp <- which(pvalsrank %in% 1:m)
  
  pred <- knn(train=train[,smallp], test = validate[,smallp], cl = ytrain, k = k)
  c(sum(yvalidate != pred), length(yvalidate))
})
predictions
tot_predict <- sum(predictions[2,])
errors <- sum(predictions[1,])
errors/tot_predict

# optimize k and m
set.seed(1)
ms <- 2^c(1:11)
ks <- seq(1,9,2)
params <- expand.grid(k=ks, m=ms)

# 10X cross-validation for each combination of m and k
error_rate <- sapply(1:length(params[,1]), function(i) {
  k <- params[i,1]
  m <- params[i,2]
  predictions <- sapply(folds, function(fold) {
    train <- X[-idx[[fold]],]
    ytrain <- y[-idx[[fold]]]
    validate <- X[idx[[fold]],]
    yvalidate <- y[idx[[fold]]]
    
    pvals <- colttests(train, ytrain)$p.value
    pvalsrank <- rank(pvals)
    smallp <- which(pvalsrank %in% 1:m)
    
    pred <- knn(train=train[,smallp], test = validate[,smallp], cl = ytrain, k = k)
    c(sum(yvalidate != pred), length(yvalidate))
  })
  tot_predict <- sum(predictions[2,])
  errors <- sum(predictions[1,])
  c(k,m,errors/tot_predict)
})
error_rate
plot(error_rate[3,])
min_rate <- which(error_rate[3,] == min(error_rate[3,]))
min_combo <- params[min_rate,]
paste("k = ", min_combo[1], " m =", min_combo[2])

# don't remove cross-validation fold before filtering
error_rate <- sapply(1:length(params[,1]), function(i) {
  k <- params[i,1]
  m <- params[i,2]
  predictions <- sapply(folds, function(fold) {
    train <- X[-idx[[fold]],]
    ytrain <- y[-idx[[fold]]]
    validate <- X[idx[[fold]],]
    yvalidate <- y[idx[[fold]]]
    
    pvals <- colttests(X, y)$p.value
    pvalsrank <- rank(pvals)
    smallp <- which(pvalsrank %in% 1:m)
    
    pred <- knn(train=train[,smallp], test = validate[,smallp], cl = ytrain, k = k)
    c(sum(yvalidate != pred), length(yvalidate))
  })
  tot_predict <- sum(predictions[2,])
  errors <- sum(predictions[1,])
  c(k,m,errors/tot_predict)
})
error_rate
which(error_rate[3,] == min(error_rate[3,]))
# result is biased with much lower error rates

# use month instead of group
y <- factor(as.numeric(format(sampleInfo$date, "%m") == "06"))
error_rate <- sapply(1:length(params[,1]), function(i) {
  k <- params[i,1]
  m <- params[i,2]
  predictions <- sapply(folds, function(fold) {
    train <- X[-idx[[fold]],]
    ytrain <- y[-idx[[fold]]]
    validate <- X[idx[[fold]],]
    yvalidate <- y[idx[[fold]]]
    
    pvals <- colttests(train, ytrain)$p.value
    pvalsrank <- rank(pvals)
    smallp <- which(pvalsrank %in% 1:m)
    
    pred <- knn(train=train[,smallp], test = validate[,smallp], cl = ytrain, k = k)
    c(sum(yvalidate != pred), length(yvalidate))
  })
  tot_predict <- sum(predictions[2,])
  errors <- sum(predictions[1,])
  c(k,m,errors/tot_predict)
})
min(error_rate[3,])
