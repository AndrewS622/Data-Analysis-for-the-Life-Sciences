url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/spider_wolff_gorb_2013.csv"
filename <- "spider_wolff_gorb_2013.csv"
library(downloader)
if (!file.exists(filename)) download(url, filename)
spider <- read.csv(filename, skip=1)

spider$log2friction <- log2(spider$friction)

boxplot(friction ~ type*leg, data=spider)
boxplot(log2friction ~ type*leg, data=spider)

fit <- lm(log2friction ~ type*leg, data=spider)
summary(fit)

anova(fit)

library(contrast)
contrast(fit, list(type = "push", leg = "L1"), list(type = "push", leg = "L2"))

## F-test
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
