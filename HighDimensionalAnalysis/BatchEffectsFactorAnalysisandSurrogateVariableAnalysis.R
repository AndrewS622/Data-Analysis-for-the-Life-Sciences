## Batch Effects
library(GSE5859Subset)
data(GSE5859Subset)

# month and sex are partially confounded
sex = sampleInfo$group
month = factor(format(sampleInfo$date,"%m"))
table(sex, month)

library(genefilter)
library(qvalue)

# compare q-values with FDR cutoff of 0.10 (due to small sample size)
tt <- rowttests(geneExpression, as.factor(sex))
qs <- qvalue(tt$p.value, fdr.level = 0.10)
sum(qs$qvalues < 0.1)

# what proportion of these genes are on X or Y chromosomes
sigIdx <- which(qs$qvalues < 0.1)
chr <- geneAnnotation$CHR[sigIdx]
mean(chr == "chrX" | chr == "chrY")

# compare autosomal genes in two different months
FPs <- sigIdx[which(chr != "chrX" & chr != "chrY")]
ttmonth <- rowttests(geneExpression[FPs,], month)
mean(ttmonth$p.value < 0.05)
# majority of autosomal genes show differences due to month

# use linear model with sex and month for each gene
X = model.matrix(~sex+month)
is = 1:dim(geneExpression)[1]
ps <- sapply(is, function(i) {
  y = geneExpression[i,]
  fit = lm(y~X-1)
  monthval <- summary(fit)$coef["Xsex", "Pr(>|t|)"]
  sexval <- summary(fit)$coef["Xmonth10", "Pr(>|t|)"]
  c(monthval, sexval)
})
qs <- qvalue(ps[1,], fdr.level = 0.10)
sum(qs$qvalues < 0.1)
sigs <- which(qs$qvalues < 0.1)
chr <- geneAnnotation$CHR[sigs]
mean(chr == "chrX" | chr == "chrY")

# now check month influence
qs <- qvalue(ps[2,], fdr.level = 0.10)
sum(qs$qvalues < 0.1)

## Factor Analysis
library(Biobase)
library(GSE5859Subset)
data(GSE5859Subset)

y <- geneExpression - rowMeans(geneExpression)
corr <- cor(y)

# plot correlation of samples
image(corr)

date <- sampleInfo$date

library(dplyr)
df <- data.frame(t(y), date = date) %>% arrange(date)

# plot ordered by date
corrOrder <- cor(t(df[,-ncol(df)]))
image(corrOrder)
# ordering by date shows clear groups

# estimate hidden factors based on first two principal components
month <- as.factor(format(date, '%m'))

pcs = svd(y)$v[,1:2]
o = order(sampleInfo$date)
cols = as.numeric(month)[o]
mypar(2,1)
for(i in 1:2){
  plot(pcs[o,i],col=cols,xaxt="n",xlab="")
  label = gsub("2005-","",sampleInfo$date[o])
  axis(1,1:ncol(y),label,las=2)
}

# which samples are most different by factor 1
pcs <- data.frame(pcs, sampleInfo$date)
colnames(pcs) <- c("PC1", "PC2", "date")
means <- pcs %>% group_by(date) %>% summarize(m = mean(PC1))

# Use SVD to obtain PCs, identify how many explain > 10% of variability
s <- svd(y)
sum(100*s$d^2/sum(s$d^2) > 10)

# which correlates most with month?
which(cor(as.numeric(month), s$v) == max(cor(as.numeric(month), s$v)))
max(cor(as.numeric(month), s$v))

# and with sex?
sex <- sampleInfo$group
which(cor(as.numeric(sex), s$v) == max(cor(as.numeric(sex), s$v)))
max(cor(as.numeric(sex), s$v))

# use PCs as basis for linear model instead of month
X <- model.matrix(~sex + s$v[,1:2])
fit <- lm(t(y)~X-1)
is = 1:dim(geneExpression)[1]
ps <- sapply(is, function(i) {
  fit = lm(y[i,]~X-1)
  summary(fit)$coef["Xsex", "Pr(>|t|)"]
})

library(qvalue)
qs <- qvalue(ps, fdr.level = 0.10)
sum(qs$qvalues < 0.1)
sigs <- which(qs$qvalues < 0.1)
chr <- geneAnnotation$CHR[sigs]
sum(chr == "chrX" | chr == "chrY")
# all detected genes are on chrX or chrY

## Surrogate Variable Analysis
library(sva)
library(Biobase)
library(GSE5859Subset)
data(GSE5859Subset)

# first factor correlates with outcome of interest
s <- svd(geneExpression-rowMeans(geneExpression))
cor(sampleInfo$group,s$v[,1])

# sva() downweights genes correlating with outcome of interest
sex <- sampleInfo$group
mod <- model.matrix(~sex)
svafit <- sva(geneExpression,mod)
head(svafit)

# still highly correlated with PCs
for (i in 1:ncol(svafit$sv)) {
  print(cor(s$v[,i], svafit$sv[,i]))
}

head(svafit$sv)

# use these variables to fit a linear model
X <- model.matrix(~sex + svafit$sv[,1:5])

is = 1:dim(geneExpression)[1]
ps <- sapply(is, function(i) {
  fit = lm(y[i,]~X-1)
  summary(fit)$coef["Xsex", "Pr(>|t|)"]
})

library(qvalue)
qs <- qvalue(ps, fdr.level = 0.10)
sum(qs$qvalues < 0.1)
sigs <- which(qs$qvalues < 0.1)
chr <- geneAnnotation$CHR[sigs]
sum(chr == "chrX" | chr == "chrY")
#12/13 genes detected on expected chromosomes