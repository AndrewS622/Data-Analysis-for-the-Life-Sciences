## Course Intro
library(UsingR)
data("father.son", package="UsingR")
mean(father.son$sheight)

fheight <- father.son$fheight
sheight <- father.son$sheight

library(dplyr)
f71 <- filter(father.son, round(father.son$fheight) == 71)
mean(f71$sheight)

## Matrix Notation
X <- matrix(1:1000, 100, 10)
X[25,3]

x <- 1:10
Y <- cbind(x, 2*x, 3*x, 4*x, 5*x)
sum(Y[7,])
matrix(1:60,20,3,byrow = TRUE)
Y %*% matrix(1,ncol(Y))

## Matrix Operations/Solving Linear Systems
data <- c(3, 4, -5, 1, 2, 2, 2, -1, 1, -1, 5, -5, 5, 0, 0, 1)
X <- t(matrix(data, nrow=4, ncol = 4))
b <- c(10, 5, 7, 4)
coef <- solve(X) %*% b
#solve() computes the inverse
coef
X %*% coef

a <- matrix(1:12, nrow=4)
b <- matrix(1:15, nrow=3)
c <- a %*% b
c
c[3,2]
sum(a[3,] * b[,2])

## Matrix Algebra Examples
# Linear models of the form X*beta
X <- matrix(c(1,1,1,1,0,0,1,1), nrow=4)
rownames(X) <- c("a", "a", "b", "b")
X
beta <- c(5,2)
X %*% beta

X <- matrix(c(1,1,1,1,1,1,0,0,1,1,0,0,0,0,0,0,1,1), nrow=6)
rownames(X) <- c("a", "a", "b", "b", "c", "c")
X
beta <- c(10,3,-3)
X %*% beta
