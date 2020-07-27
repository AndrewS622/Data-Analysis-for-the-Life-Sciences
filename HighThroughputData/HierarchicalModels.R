## Bayes Rule
# cystic fibrosis test
pposD <- 0.99
pnegnoD <- 0.99
pD <- 0.00025
pnoD <- 1 - pD
pposnoD <- 1 - pnegnoD

ppos <- pposD*pD + pposnoD*pnoD
pDpos <- pposD*pD/ppos
pDpos

## Bayes Rule in Practice
tmpfile <- tempfile()
tmpdir <- tempdir()
download.file("http://seanlahman.com/files/database/lahman-csv_2014-02-14.zip",tmpfile)
filenames <- unzip(tmpfile,list=TRUE)
players <- read.csv(unzip(tmpfile,files="Batting.csv",exdir=tmpdir),as.is=TRUE)
unlink(tmpdir)
file.remove(tmpfile)

library(dplyr)
avgs <- filter(players,yearID==2012 | yearID==2010 | yearID == 2011) %>% 
  mutate(AVG=H/AB) %>% 
  filter(AB>=500) %>% 
  select(AVG) %>% 
  unlist()
mean(avgs)
sd(avgs)
qqnorm(avgs)
qqline(avgs)

# Jose Iglesias
p <- 0.450
N <- 20
sdEst <- sqrt(N*p*(1-p))/N
sdEst

# parameters for prior
mu <- 0.275
tau <- 0.027

# hierarchical model
beta <- sdEst^2/(sdEst^2 + tau^2)
posteriorMean <- beta*mu + (1-beta)*p
posteriorMean

## Hierarchical Models in Practice
library(Biobase)
library(SpikeInSubset)
library(genefilter)
data("rma95")
y <- exprs(rma95)

# 6 replicate samples of RNA w/ 16 genes artificially added in diff quantities
# artificial quantities specified in pData(rma95)
pData(rma95)
g <- factor(rep(0:1,each=3))
spike <- rownames(y) %in% colnames(pData(rma95))
ttests <- rowttests(y, g)

sig <- ttests[which(ttests$p.value < 0.01),]
positives <- colnames(pData(rma95))
TPs <- sig[which(rownames(sig) %in% positives),]
FPR <- (dim(sig)[1]-dim(TPs)[1])/dim(sig)[1]
FPR
FPs <- sig[which(!(rownames(sig) %in% positives)),]

nosig <- ttests[which(ttests$p.value >= 0.01),]
FNs <- nosig[which(rownames(nosig) %in% positives),]
TNs <- nosig[which(!(rownames(nosig) %in% positives)),]

sds <- rowSds(y[,g==1])
sdsFP <- sds[names(sds) %in% rownames(FPs)]
sdsTP <- sds[names(sds) %in% rownames(TPs)]
sdsFN <- sds[names(sds) %in% rownames(FNs)]
sdsTN <- sds[names(sds) %in% rownames(TNs)]
boxplot(sdsTP, sdsFP, sdsTN, sdsFN)

# alternative method to obtain same plot
sds <- rowSds(y[,g==1])
index <- paste0(as.numeric(spike), as.numeric(ttests$p.value<0.01))
index <- factor(index,levels=c("11","01","00","10"),labels=c("TP","FP","TN","FN"))
boxplot(split(sds,index))
# p-values can be small by chance

# use empirical Bayes hierarchical model
library(limma)
fit <- lmFit(y, design=model.matrix(~g))
colnames(coef(fit))
fit <- eBayes(fit)
sampleSD <- fit$sigma
posteriorSD <- sqrt(fit$s2.post)

# hierarchical adjusted SD vs. sample-based estimate
plot(sampleSD, posteriorSD)

# compute adjusted p-values
pvals <- fit$p.value[,2]
sigNew <- pvals[which(pvals < 0.01)]
TPs <- sigNew[which(names(sigNew) %in% positives)]
FPR <- (length(sigNew)-length(TPs))/length(sigNew)
FPR
