## Differential expression and multiple comparisons
library(Biobase)
library(maPooling)
data(maPooling)
pd=pData(maPooling)

# pooled entries contain RNA from all 12 samples
# individuals contain just one
pooled=which(rowSums(pd)==12)
individuals=which(rowSums(pd)==1)

# remove technical replicates
individuals=individuals[-grep("tr",names(individuals))]

# extract data
pool = exprs(maPooling)[,pooled] 
indiv = exprs(maPooling)[,individuals]
strain = ifelse(grepl("a",rownames(pData(maPooling))),0,1)

g_pool = strain[pooled]
g_indiv = strain[individuals]

# technical and biologlical variability
library(genefilter)
bio_var <- rowSds(indiv[,g_indiv == 1])
tech_var <- rowSds(pool[,g_pool == 1])
mean(bio_var > tech_var)

# t-tests and false-discovery rate
tt <- rowttests(pool, as.factor(g_pool))
library(qvalue)
q <- qvalue(tt$p.value)
sum(q$qvalues < 0.05)

# repeat using biological replicates on same population of genes
genes <- which(q$qvalues < 0.05)
tt_bio <- rowttests(indiv[genes,], as.factor(g_indiv))
q_bio <- qvalue(tt_bio$p.value)
mean(q_bio$qvalues > 0.05)

# now use moderated t-tests 
library(limma)
fit <- lmFit(indiv[genes,], design=model.matrix(~as.factor(g_indiv)))
fit <- eBayes(fit)
tt <- topTable(fit, coef=2, number=length(genes))
q_mod <- qvalue(tt$adj.P.Val)
mean(q_mod$qvalues < 0.05)

## Gene set analysis: Inference on coordinated expression changes
# variance of average of N standard normal variables = 1/N
var(rowMeans(matrix(rnorm(10000 * 10, 0, 1), ncol=10)))

# if a correlation is present
library(MASS)
Sigma = matrix(0.7, ncol=10, nrow=10)
diag(Sigma) = 1
num <- mvrnorm(n=10000,mu=rep(0,10),Sigma=Sigma)
var(rowMeans(num))
# variance increases with correlation, so null must be widened

# ROAST algorithm for differential gene set expression
library(limma)
library(qvalue)
library(GEOquery)
e = getGEO("GSE34313")[[1]]

# data is from human airway smooth muscle treated with glucocorticoid to reduce inflammation
pData(e)$condition = factor(pData(e)$characteristics_ch1.2)
levels(pData(e)$condition) = c("dex24", "dex4", "control")

# gene ontology terms for each gene
names(fData(e))
fData(e)$GO_ID[1:4]

# compare only control and 4-hr dexamethasone treatments
lvls = c("control", "dex4")
es = e[, e$condition %in% lvls]
es$condition = factor(es$condition, levels = lvls)

# perform differential analysis
design = model.matrix(~ es$condition)
fit <- lmFit(es, design = design)
fit <- eBayes(fit)
topTable(fit)[,c(6,7,18,22)]
# top genes are immune response (e.g. CSF2, IL6, CCL2, LIF)

# use ROAST for gene set testing based on GO_ID
# uses rotation of residuals to generate null distribution for summary of scores from each gene
# algorithm is self-contained (only relies on scores from the gene set being tested)
set.seed(1)
idx <- grep("GO:0006955", fData(es)$GO_ID)
length(idx)
r1 <- roast(es, idx, design)
r1

set.seed(1)
idx <- grep("GO:0045454", fData(es)$GO_ID)
length(idx)
r2 <- roast(es, idx, design)
r2

# mroast performs multiple tests for different gene sets
library(org.Hs.eg.db)
org.Hs.egGO2EG
go2eg <- as.list(org.Hs.egGO2EG)
head(go2eg)

# match Entrez ID to index in ExpressionSet
govector <- unlist(go2eg)
golengths <- sapply(go2eg, length)
head(fData(es)$GENE)

idxvector <- match(govector, fData(es)$GENE) 
table(is.na(idxvector))

# organize list of indexes for genes per GO term
idx <- split(idxvector, rep(names(go2eg), golengths))

# to access genes
go2eg[[1]]
fData(es)$GENE[idx[[1]]]

# remove NAs and small gene sets
idxclean <- lapply(idx, function(x) x[!is.na(x)])
idxlengths <- sapply(idxclean, length)
idxsub <- idxclean[idxlengths > 10]
length(idxsub)

# run mroast test
set.seed(1)
r3 <- mroast(es, idxsub, design)
head(r3)
r3[which.max(r3$PropUp),]

# rerun using sets with at least 50 genes
idxsub <- idxclean[idxlengths >= 50]
length(idxsub)
set.seed(1)
r4 <- mroast(es, idxsub, design)
head(r4)
r4[which.max(r4$PropUp),]
GID <- rownames(r4[which.max(r4$PropUp),])

# select this annotation to find the term
library(GO.db)
library(ensembldb)
ensembldb::select(GO.db, key = GID, keytype = "GOID", columns = "TERM")
