## Hierarchical Clustering
set.seed(1)
m <- 10000
n <- 24

# create random matrix and cluster
x <- matrix(rnorm(m*n), m, n)
colnames(x) <- 1:n

d <- dist(t(x))
hc <- hclust(d)
plot(hc)

# repeat 100 times
# cluster number at a given height is a random variable
set.seed(1)
B <- 100
clustRep <- function() {
  x <- matrix(rnorm(m*n), m, n)
  d <- dist(t(x))
  hc <- hclust(d)
  cutree(hc, h = 143)
}
clusts <- replicate(B, clustRep())
maxs <- apply(clusts, 2, max)
sd(maxs)

## K-means Clustering
library(rafalib)
library(GSE5859Subset)
data(GSE5859Subset)
set.seed(10)
km <- kmeans(t(geneExpression), 5)
km$cluster

# clusters mostly in groups with date
table(cluster = km$cluster, true = sampleInfo$date)
table(cluster = km$cluster, true = sampleInfo$group)
table(cluster = km$cluster, true = sampleInfo$ethnicity)

## Heat Maps
# find genes with highest across-sample variance (using MAD)
library(matrixStats)
idx <- order(rowMads(geneExpression), decreasing = TRUE)[1:25]
library(gplots)
heatmap.2(scale(geneExpression[idx,]), trace = "none", ColSideColors = as.character(sampleInfo$group), labCol =  sampleInfo$date, labRow =  geneAnnotation[idx,]$CHR)

# random data independent of group
set.seed(17)
m <- nrow(geneExpression)
n <- ncol(geneExpression)
x <- matrix(rnorm(m*n), m, n)
g <- factor(sampleInfo$group)

library(rafalib)
mypar(1,2)
library(genefilter)

# examine genes with smallest p-values and largest variances
tt <- rowttests(x)
idx1 <- order(tt$p.value)[1:50]
sds <- sqrt(rowVars(x))
idx2 <- order(sds, decreasing=TRUE)[1:50]
heatmap.2(x[idx1,], ColSideColors = as.character(g), trace = "none")
heatmap.2(x[idx2,], ColSideColors = as.character(g), trace = "none")

# results show no trend with group (as expected)
# suggests that some p-values are just significant by chance