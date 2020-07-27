library(SpikeInSubset)
data(mas133)

e <- exprs(mas133)
plot(e[,1],e[,2],main=paste0("corr=",signif(cor(e[,1],e[,2]),3)),cex=0.5)
k <- 3000
b <- 1000 #a buffer
# correlation between first two samples
polygon(c(-b,k,k,-b),c(-b,-b,k,k),col="red",density=0,border="red")

# points in red box
sum(e[,1] < k & e[,2] < k)/length(e[,1])

# plot log of values so tails do not dominate
plot(log2(e[,1]),log2(e[,2]),main=paste0("corr=",signif(cor(log2(e[,1]),log2(e[,2])),2)),cex=0.5)
k <- log2(k)
b <- log2(0.5)
polygon(c(b,k,k,b),c(b,b,k,k),col="red",density=0,border="red")

sum(log2(e[,1]) < k & log2(e[,2]) < k)/length(e[,1])

# MA plot
e <- log2(exprs(mas133))
plot((e[,1]+e[,2])/2,e[,2]-e[,1],cex=0.5)

# standard deviation of log ratios
sd(e[,2]-e[,1])

# > 2-fold changes
sum(abs(e[,2]-e[,1]) > 1)
