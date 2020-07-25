## Linear Models as a t-test
# Numerator of both is difference in means
# Need to show sigma^2 (Xt X)^-1 = 
# (1/nx + 1/ny) * (sum {(xi - mu_x)^2} + 
# sum{(yi - mu-y)^2})/(nx + ny -2)
# assuming equal variance
# sigma^2 = SSR / (N-2) = (sum {(xi - mu_x)^2} + 
# sum{(yi - mu-y)^2})/(nx + ny -2)
# Show (Xt X)^-1 [2,2] = (1/nx + 1/ny)
nx <- 5
ny <- 7
X <- cbind(rep(1,nx + ny), 
           rep(c(0,1),
           c(nx, ny)))

t(X) %*% X
# From formula for 2x2 determinant, [2,2] = 
# 12/(12*7-7*7) = 1/5 + 1/7