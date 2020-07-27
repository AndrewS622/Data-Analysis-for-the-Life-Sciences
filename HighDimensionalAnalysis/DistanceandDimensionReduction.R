## Distances
# RNA expression levels for samples from seven tissues
library(tissuesGeneExpression)
data(tissuesGeneExpression)
head(e)
head(tissue)
sum(tissue == "hippocampus")

# compute distance between samples
d <- dist(t(e))
as.matrix(d)[3,45]

# distance between genes
g1 <- e[which(rownames(e) == "210486_at"),]
g2 <- e[which(rownames(e) == "200805_at"),]
sqrt(crossprod(g1 - g2))

# matrix of gene-gene distances would be too large
dim(e)
dim(e)[1]^2

length(d)
# dist() utilizes symmetry, storing only lower triangular matrix
# length is then ncol(e) * (ncol(e)-1)/2

## Projections
library(Biobase)
library(GSE5859Subset)
data(GSE5859Subset)

y <- geneExpression[,1:2]
A <- matrix(c(1,1,1,-1), nrow = 2)
MA <- y %*% A

## Singular Value Decomposition
library(tissuesGeneExpression)
data(tissuesGeneExpression)
s <- svd(e)

# can flip sign of any combination of columns and SVD will not change
# as long as both U and V columns signs are flipped
signflips <- sample(c(-1, 1), ncol(e), replace=TRUE)
signflips
newu <- sweep(s$u, 2, signflips, FUN="*")
newv <- sweep(s$v, 2, signflips, FUN="*")
all.equal(s$u %*% diag(s$d) %*% t(s$v), newu %*% diag(s$d) %*% t(newv))

# first column relates to mean of rows of e
m <- rowMeans(e)
cor(s$u[,1], m)

# changing means of rows does not change distances
newmeans <- rnorm(nrow(e))
newe <- e + newmeans
sqrt(crossprod(e[,3]-e[,45]))
sqrt(crossprod(newe[,3]-newe[,45]))

# thus, can subtract means off of each row
y <- e - rowMeans(e)
s <- svd(y)
resid <- y - s$u %*% diag(s$d) %*% t(s$v)
max(abs(resid))

# x*a = x[i,] * a[i] for matrix x, vector a
x <- matrix(rep(c(1,2),each=5),5,2)
x
x*c(1:5)
sweep(x,1,1:5,"*")

# creation of diagonal matrix for d is unnecessary
diag(s$d) %*% t(s$v)
s$d * t(s$v)

vd <- t(s$d * t(s$v))

# the following are equivalent to Y (= UDVt)
s$u %*% (s$d * t(s$v))
tcrossprod(s$u, vd)
s$u %*% t(vd)

z <- s$d * t(s$v)
identical(z, t(vd))

# all are equivalent
sqrt(crossprod(e[,3]-e[,45]))
sqrt(crossprod(y[,3]-y[,45]))
sqrt(crossprod(z[,3]-z[,45]))

# approximation using only two dimensions of z
sqrt(crossprod(z[1:2,3]-z[1:2,45])) - sqrt(crossprod(z[,3]-z[,45]))

# find min number of dimensions for approximation to be >90% accurate
r <- -100*(sqrt(crossprod(z[1:2,3]-z[1:2,45])) - sqrt(crossprod(z[,3]-z[,45])))/sqrt(crossprod(z[,3]-z[,45]))
i = 2
while(r > 10) {
  i = i + 1
  r <- -100*(sqrt(crossprod(z[1:i,3]-z[1:i,45])) - sqrt(crossprod(z[,3]-z[,45])))/sqrt(crossprod(z[,3]-z[,45]))
}
paste(i, "dimensions yields error of", r)

distances <- sqrt(apply(z[,-3]-z[,3],2,crossprod))
distances2 <- sqrt(apply(z[1:2,-3]-z[1:2,3],2,crossprod))
cor(distances2,distances, method="spearman")

## Multi-dimensional Scaling
library(rafalib)
ftissue <- factor(tissue)
mypar(1,1)
plot(z[1,],z[2,],col=as.numeric(ftissue))
legend("topleft",levels(ftissue),col=seq_along(ftissue),pch=1)

mds <- cmdscale(dist(t(e)))
cor(z[1,], mds[,1])
cor(z[2,], mds[,2])
# first two rows of z correlate (inversely) w/ first two columns of mds

mypar(1,2)
plot(z[1,],z[2,],col=as.numeric(ftissue))
legend("topleft",levels(ftissue),col=seq_along(ftissue),pch=1)
plot(mds[,1],mds[,2],col=as.numeric(ftissue))

# inverting sign of z produces identical plot
plot(-z[1,],-z[2,],col=as.numeric(ftissue))
plot(mds[,1],mds[,2],col=as.numeric(ftissue))

s <- svd(geneExpression - rowMeans(geneExpression))
z <- s$d * t(s$v)

# dimension 1 correlates most with group
which(cor(t(z),sampleInfo$group) == max(abs(cor(t(z),sampleInfo$group))))

#second-highest dimension is 6
which(rank(-abs(cor(t(z), sampleInfo$group))) == 2)

sampleInfo$date
month <- format(sampleInfo$date, "%m")
month <- factor(month)

# also correlates most with date
which(cor(t(z),as.numeric(month)) == max(abs(cor(t(z),as.numeric(month)))))

# stratify U[,6] by chromosome, remove chrUn, and make boxplot
mypar(1,1)
table(geneAnnotation$CHR)
boxplot(s$u[-which(geneAnnotation$CHR == "chrUn"),6] ~ geneAnnotation$CHR[-which(geneAnnotation$CHR == "chrUn")])

# Suggests that sampleInfo$group represents sex