## Introduction to DNA methylation
library(BSgenome.Hsapiens.UCSC.hg19)
library(Biostrings)

# subset human genome
chr22 <- Hsapiens[["chr22"]]
s <- subseq(chr22, start = 23456789, width=1000)
print(as.character(s))

# GC content
letterFrequency(s, "GC", as.prob = TRUE)

# CpGs
countPattern("CG", s)

# GpCs
countPattern("GC", s)

# CpG islands
library(AnnotationHub)
ah <- AnnotationHub()
ah <- subset(ah, ah$genome=="hg19")
query(ah, "genes")
cgi_id <- names(query(ah, "CpG Islands"))
cgi <- ah[[cgi_id]]
class(cgi)
length(cgi)

# sequences of each CpG island
library(BSgenome.Hsapiens.UCSC.hg19)
cgiseq <- getSeq(Hsapiens, cgi)
identical(genome(cgi)[1:24], genome(Hsapiens)[1:24])

# Median letter content for islands
pcs <- letterFrequency(cgiseq, "C", as.prob = TRUE)
pgs <- letterFrequency(cgiseq, "G", as.prob = TRUE)
median(pcs)
median(pgs)
expected <- pcs*pgs*width(cgiseq)
observed <- vcountPattern("CG", cgiseq)
cpgoe <- observed/expected
median(cpgoe)
observedGpC <- vcountPattern("GC", cgiseq)
gpcoe <- observedGpC/expected
median(gpcoe)

# in other parts of the genome, restricted to mapped chromosomes
chr2use <- seqlevels(cgi)[1:24]
index <- which(seqnames(cgi) %in% chr2use)
noncgi <- shift(cgi[index], 20000)
noncgiseq <- getSeq(Hsapiens, noncgi)
nullres <- alphabetFrequency(noncgiseq)
keepIndex <- nullres[,"G"]>0 &  nullres[,"C"]>0 & nullres[,"N"]==0
nullres <- nullres[keepIndex,]
noncgiseq <- noncgiseq[keepIndex]

# compute CpG o/e
pcs <- letterFrequency(noncgiseq, "C", as.prob = TRUE)
pgs <- letterFrequency(noncgiseq, "G", as.prob = TRUE)
expected2 <- pcs*pgs*width(noncgiseq)
observed2 <- vcountPattern("CG", noncgiseq)
noncgioe <- observed2/expected2
median(noncgioe)

## DNA Methylation Measurement
# Finding differentially methylated regions
library(devtools)
install_github("genomicsclass/coloncancermeth")
library(coloncancermeth)
data(coloncancermeth)

# dimensions
dim(meth)
dim(pd)
print(gr)

# data structure
sum(pd$Status == "cancer")
which(pd$Status == "cancer" & pd$bcr_patient_barcode == "TCGA-A6-4107")

# Euclidean distance
d <- dist(t(meth))
mds <- cmdscale(d)
plot(mds[,1], mds[,2], bg = as.numeric(pd$Status), pch = 21)
legend("topleft", levels(pd$Status),col=seq(along=levels(pd$Status)), pch=15)

# compute p-values from limma
library(limma)
X <- model.matrix(~pd$Status)
fit <- lmFit(meth, X)
eb <- eBayes(fit)
pvals <- eb$p.value[,2]

# compute q-values
library(qvalue)
qvals <- qvalue(pvals)
mean(qvals$qvalues < 0.05)

# hypermethylated regions in cancer
idx <- which(qvals$qvalues < 0.05)
mean(fit$coefficients[idx,2] > 0)

# which differentially methylated CpGs are CpG islands?
library(AnnotationHub)
ah <- AnnotationHub()
cgi <- ah[["AH5086"]]
mean(gr[idx] %over% cgi)

# using bumphunter to separate CpGs into groups
library(bumphunter)
X <- model.matrix(~pd$Status)
chr <- as.character(seqnames(gr))
res <- bumphunter(meth, X, chr=chr, pos = start(gr), cutoff=0.1)
head(res$table)
# B argument can be used to assess uncertainty

# filter by region size instead due to computational intensity
dmrs <- res$table[res$table$L >= 3,]
dmrs <- makeGRangesFromDataFrame(dmrs)

# distances to nearest CpG island
d <- distanceToNearest(dmrs, cgi)

# overlapping CpG
mean(mcols(d)$distance == 0)

# within 2 kbp but not overlapping
mean(mcols(d)$ distance <= 2000 & mcols(d)$distance > 0)

# Reading Raw 450K Array Data
path <- "D:/Harvard/Other/CaseStudiesinFunctionalGenomics/rawdata-master/idats/"
list.files(path)

# read in targets CSV
targets <- read.csv(file.path(path,"targets.csv"), as.is = TRUE)
names(targets)
targets$Basename
sum(targets$Status == "cancer")
targets$Basename <- file.path(path,targets$Basename)

# read raw data
library(minfi)
rgset <- read.metharray(targets$Basename, verbose=TRUE)
rownames(targets) <- sampleNames(rgset)
targets <- as(targets, "DataFrame")
pData(rgset) <- targets
dim(getRed(rgset))
dim(getGreen(rgset))

# preprocessing and map to genome
mset <- preprocessIllumina(rgset)
mset <- mapToGenome(mset)

# obtain methylation values and CpG locations
dim(getBeta(mset, type = "Illumina"))
head(granges(mset))
methyl <- getBeta(mset, type="Illumina")
gr <- granges(mset)
idxgr <- names(gr)[which(seqnames(gr) == "chr4" & start(gr) == 153807318)]
methyl[idxgr, "5775041068_R04C01"]

# estimating CpG islands
library(bumphunter)
class(mset)
showMethods("bumphunter")
grset <- ratioConvert(mset, what="beta", type="Illumina")

# run bumphunter
dmrs <- bumphunter(grset, model.matrix(~pData(grset)$Status), cutoff=0.1)
head(dmrs$table)
dmrs$table$area[1]

# bumphunter with smoothing
index <- which(seqnames(grset) == "chr22")
grset2 <- grset[index,]
X <- model.matrix(~pData(grset2)$Status)
res <- bumphunter(grset2, X, cutoff = 0.25)
res2 <- bumphunter(grset2, X, cutoff = 0.25, smooth = TRUE)
tab1 <- res$table
tab2 <- res2$table
mean(tab1$L)
mean(tab2$L)
length(tab1$chr)
length(tab2$chr)

## Data Analysis and Integration
# CpG Island Shores
library(minfi)
library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# path to data
path <- "D:/Harvard/Other/CaseStudiesinFunctionalGenomics/tcgaMethylationSubset"
targets <- read.delim(file.path(path, "targets.txt"), as.is=TRUE)

dim(targets)
head(targets)
table(targets$Tissue, targets$Status)

# Find DMRs only from normal breast and colon samples
index <- which(targets$Status=="normal" & targets$Tissue%in%c("colon", "breast"))
targets <- targets[index,]

# read in data
dat <- read.metharray.exp(base=path, targets=targets, verbose=TRUE)

# convert to find methylation values
class(dat)
dat <- preprocessIllumina(dat)
class(dat)

# assign locations to CpGs
dat <- mapToGenome(dat)
class(dat)

# compute methylation values
dat <- ratioConvert(dat, type="Illumina")
class(dat)

# distribution of samples
library(rafalib)
mypar(1,1)
y <- getBeta(dat)
shist(y)
# smooth histogram does not show substantial differences

# MDS plot to search for outliers
mds <- cmdscale(dist(t(y)))
tissue <- as.factor(pData(dat)$Tissue)
cols <- ifelse(tissue == "breast", "red", "blue")
plot(mds, col=cols)

# find DMRs using inference
library(limma)
X <- model.matrix(~tissue)
fit <- lmFit(y, X)
head(fit$coefficients)

# which CpG has largest effect size
idx <- which.max(abs(fit$coefficients[,2]))
cpg <- rownames(fit$coefficients)[idx]
granges(dat)[cpg,]

# determine q-value
library(qvalue)
eb <- eBayes(fit)
qvals <- qvalue(eb$p.value[,2])$qvalue
-log10(qvals[cpg])

# find all CpGs within 5 kbp
gr <- granges(dat)[cpg,] + 5000
cpgs <- granges(dat)[which(granges(dat)%over%gr)]
cpgnames <- names(cpgs)
pos <- start(cpgs)

# plot methylation values, effect sizes and -log10(qvals)
matplot(pos, y[cpgnames,], col=cols)
plot(pos, fit$coefficients[cpgnames,2])
plot(pos, -log10(qvals[cpgnames]))

# repeat for top 10 CpGs ranked by effect size
o <- order(abs(fit$coefficients[,2]), decreasing=TRUE)[1:10]
mypar(1,3)
for(i in 1:10) {
  cpgi <- rownames(dat)[o[i]]
  gr <- granges(dat)[cpgi,] + 5000
  cpgs <- granges(dat)[which(granges(dat)%over%gr)]
  cpgnames <- names(cpgs)
  pos <- start(cpgs)
  matplot(pos, y[cpgnames,], col=cols)
  plot(pos, fit$coefficients[cpgnames,2])
  plot(pos, -log10(qvals[cpgnames]))
  invisible(readline(prompt="Press [enter] to continue"))
}

# use bumphunter to search for regions
# restrict to chr15 for speed
index <- which(seqnames(dat)=="chr15")
dat2 <- dat[index,]
library(doParallel)
ncores <- detectCores()
registerDoParallel(cores = ncores)
library(bumphunter)
set.seed(1)
res <- bumphunter(dat2, X, cutoff=0.1, B=100)
head(res$tab)
sum(res$tab$fwer < 0.05)

# compare results from bumphunter and CpG analysis
index <- which(qvals < 0.05 & abs(fit$coefficients[,2]) > 0.5 & seqnames(dat) == "chr15")
tab <- res$tab[res$tab$L >= 3,]
tab <- makeGRangesFromDataFrame(tab)
head(tab)
cpgranges <- granges(dat)[index,]
mean(cpgranges%over%tab)

# use annotations to find locations from true CpG islands
library(AnnotationHub)
cgi <- AnnotationHub()[["AH5086"]]
tab <- res$tab[res$tab$fwer < 0.05,]
tab <- makeGRangesFromDataFrame(tab)
dists <- distanceToNearest(tab, cgi)
mean(mcols(dists)$distance > 0 & mcols(dists)$distance < 2000)

# Cell Composition
grset <- getGenomicRatioSetFromGEO("GSE32148")
save(grset,file="grset.rda")
#load("grset.rda")

class(grset)
names(pData(grset))
head(pData(grset))

# group information (Crohn's vs. controls)
group <- rep("normal", nrow(pData(grset)))
group[grepl("ulcerative", pData(grset)[,1])]="ulcerative"
group[grepl("Crohn", pData(grset)[,1])]="crohn"
group <- factor(group, levels=c("normal", "ulcerative", "crohn"))

# remove CpGs with NAs and sex chromosomes
keep <- which(rowSums(is.na(getBeta(grset)))==0 &
                !seqnames(grset)%in% c("chrX", "chrY"))
grset2 <- grset[keep,]

# extract age
age <- pData(grset2)$'age (y):ch1'
age <- as.character(age)
age[age == "N/A"] <- NA
age <- as.numeric(age)

# MDS plot
d <- dist(t(getBeta(grset2)))
mds <- cmdscale(d)
cols <- ifelse(age > 40, "red", "blue")
plot(mds[,1], mds[,2], col=cols, pch = as.numeric(group))
legend("bottomleft", legend = c("normal", "ulcerative", "crohn"), pch = 1:3)

# Multi-resolution Analysis
library(minfi)
# use same path as CpG Islands
targets <- read.delim(file.path(path, "targets.txt"), as.is = TRUE)
index <- which(targets$Tissue == "colon")
targets <- targets[index,]
dat <- read.metharray.exp(base=path, targets=targets, verbose=TRUE)
dat <- preprocessIllumina(dat)
dat <- mapToGenome(dat)

# collapse the data
cdat <- cpgCollapse(dat)
class(cdat)
names(cdat)
nrow(cdat$object)

# types of regions represented
head(granges(cdat$object))
mean(granges(cdat$object)$type == "OpenSea")

# find DMRs between cancer and normal
status <- factor(pData(cdat$object)$Status,
                 level = c("normal", "cancer"))
X <- model.matrix(~status)
res <- blockFinder(cdat$object, X, cutoff=0.05)
head(res$table)
mean(res$table$value < 0)

## Whole Genome Analysis
# Measuring Methylation from Sequencing
library(BSgenome.Hsapiens.UCSC.hg19)
chr22 <- Hsapiens[["chr22"]]

# predict sizes of fragments for reduced representation bisulfite sequencing (RRBS)
CCGG <- matchPattern("CCGG", chr22)
l <- length(CCGG)

# fragments (add 3 to account for length of cut sequence)
fragments <- start(CCGG)[2:l] - end(CCGG)[1:l-1] + 3
hist(fragments, breaks = 100)
hist(fragments[fragments < 10000])
hist(log10(fragments))

# how many between 40 and 220 bp
mean(fragments >= 40 & fragments <= 220)
sum(fragments[fragments >= 40 & fragments <= 220])

# Whole-Genome Bisulfite Sequencing (WGBS)
# note: the file downloaded from the GitHub repo had space-separated columns 
# but spaces were also present in some string entries
# the data was pasted into Excel and exported as a .csv instead
path <- "D:/Harvard/Other/CaseStudiesinFunctionalGenomics/colonCancerWGBS"
targets <- read.table(file.path(path, "targets.txt"), header=TRUE, sep=",")

library(bsseq)

# specify coverage files and read them using read.bismark()
targets <- DataFrame(targets, row.names = as.character(targets$Run))
covfiles <- file.path(path, paste0(rownames(targets), ".chr22.cov"))
colonCancerWGBS <- read.bismark(files = covfiles, rmZeroCov = TRUE, 
                                colData = targets)

# examine dataset object
colonCancerWGBS
pData(colonCancerWGBS)
granges(colonCancerWGBS)

# extract coverage and methylation
cov <- getCoverage(colonCancerWGBS, type="Cov")
m <- getCoverage(colonCancerWGBS, type="M")

# proportion of CpGs with coverage in all samples
mean(rowSums(cov > 0) == 6)

# plot total coverage for each site vs. location
loc = start(granges(colonCancerWGBS))
totCov = rowSums(cov)
hist(totCov)
plot(loc, totCov, ylim = c(0, 2000))

# selected region of coverage
gr = GRanges(seqnames="22",ranges=IRanges(start=43793678,end= 45022550))
index=granges(colonCancerWGBS)%over%gr
library(rafalib)
i=1
index2=which(index & cov[,i]>=5 & cov[,i]<=50)
x=start(colonCancerWGBS)[index2]
y=m[index2,i]/cov[index2,i]
w=sqrt(cov[index2,i])/7
plot(x,y,cex=w)
