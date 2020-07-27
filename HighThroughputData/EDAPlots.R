BiocManager::install("SpikeInSubset")
library(SpikeInSubset)
data(mas133)

e <- exprs(mas133)
plot(e[,1],e[,2],main=paste0("corr=",signif(cor(e[,1],e[,2]),3)),cex=0.5)
k <- 3000
b <- 1000 #a buffer
polygon(c(-b,k,k,-b),c(-b,-b,k,k),col="red",density=0,border="red")

sum(e[,1] < 3000 & e[,2] < 3000)/length(e[,1])

plot(log2(e[,1]),log2(e[,2]),main=paste0("corr=",signif(cor(log2(e[,1]),log2(e[,2])),2)),cex=0.5)
k <- log2(3000)
b <- log2(0.5)
polygon(c(b,k,k,b),c(b,b,k,k),col="red",density=0,border="red")

sum(log2(e[,1]) < log2(3000) & log2(e[,2]) < log2(3000))/length(e[,1])


e <- log2(exprs(mas133))
plot((e[,1]+e[,2])/2,e[,2]-e[,1],cex=0.5)
sd(e[,2]-e[,1])
sum(abs(e[,2]-e[,1]) > 1)
