## Collinearity
sex <- factor(rep(c("female","male"), each=4))
trt <- factor(c("A", "A", "B", "B", "C", "C", "D", "D"))
X <- model.matrix(~ sex + trt)
# matrix is not full rank --> collinear
c(paste("rank = ", qr(X)$rank), paste("nrow = ", nrow(X)))

# outcome observed
Y <- 1:8

# Fix two coefficients by subtracting off their contributions
makeYstar <- function(a,b) Y-X[,2] * a - X[,5] * b
# Once fixed, optimize the rest by minimizing the SSE
fitTheRest <- function(a,b) {
  Ystar <- makeYstar(a,b)
  Xrest <- X[,-c(2,5)]
  betarest <- solve(t(Xrest) %*% Xrest) %*% t(Xrest) %*% Ystar
  residuals <- Ystar - Xrest %*% betarest
  sum(residuals^2)
}

fitTheRest(1,2)

betas <- expand.grid(-2:8, -2:8)
rss <- apply(betas,1,function(x) fitTheRest(x[1],x[2]))
rss
min(rss)
index <- which (rss %in% min(rss))
betas[index,]
# multiple beta combinations show minimum due to collinearity

library(rafalib)
themin <- min(rss)
plot(betas[which(rss==themin),])

## QR Decomposition
url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/spider_wolff_gorb_2013.csv"
filename <- "spider_wolff_gorb_2013.csv"
library(downloader)
if (!file.exists(filename)) download(url, filename)
spider <- read.csv(filename, skip=1)

fit <- lm(friction ~ type + leg, data=spider)
betahat <- coef(fit)

# Use of QR decomposition is less prone to numerical instability
Y <- matrix(spider$friction, ncol=1)
X <- model.matrix(~ type + leg, data=spider)
QR <- qr(X)
Q <- qr.Q(QR)
R <- qr.R(QR)

Q[1,1]
R[1,1]
(t(Q) %*% Y)[1,1]
solve(R) %*% t(Q) %*% Y
fit$coefficients

# Instead of using solution beta-hat = (XtX)^-1 * XtY
# decompose X to QR,and solution is R^-1 (QtY)
