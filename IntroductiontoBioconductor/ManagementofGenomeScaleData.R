## ExpressionSet
library(GSE5859Subset)
library(Biobase)
library(GEOquery)
data(GSE5859Subset)
dim(geneExpression)
dim(sampleInfo)
dim(geneAnnotation)
?ExpressionSet

# check to ensure compatibility between data frames
identical(colnames(geneExpression), sampleInfo$filename)
identical(rownames(geneExpression), geneAnnotation$PROBEID)

# in ExpressionSet, rownmaes of assayData must match rownames of featureData
# rownames of phenoData must match colnames of assayData

# make AnnotatedDataFrame for phenoData
pd <- AnnotatedDataFrame(sampleInfo)
rownames(pd) <- colnames(geneExpression)

# find date for given sample
pData(pd)["GSM136530.CEL.gz", "date"]

# list variable names
varLabels(pd)

# repeat for featureData
fd <- AnnotatedDataFrame(geneAnnotation)
rownames(fd) <- rownames(geneExpression)
pData(fd)["204810_s_at", "CHR"]

# create object
eset <- ExpressionSet(geneExpression, pd, fd)

# explore expression differences on Y chromosome in different sexes
ind1 <- which(featureData(eset)$CHR == "chrY")
ind2 <- pData(eset)$group == 1

femaleY <- colMeans(exprs(eset)[ind1, ind2])
maleY <- colMeans(exprs(eset)[ind1, !ind2])

boxplot(maleY, femaleY)
median(maleY) - median(femaleY)

# subset first 10 features and 5 samples
eset[1:10, 1:5]

## Reading microarray data
# importing Affymetrix data from CEL files
# CEL files downloaded from course repo on Github
library("hgu95acdf")
library("affy")

# read files
wd <- getwd()
datadir <- paste0(wd, "/rawdata-master")
basedir <- paste0(datadir, "/celfiles")
setwd(basedir)
f <- read.delim("sampleinfo.txt")
rownames(f) <- f$filenames

# read CEL files 
fns <- list.celfiles(basedir)
fns %in% f[,1]
ab <- ReadAffy(phenoData = f)

# extract feature IDs
sum(probeNames(ab) == "36311_at")

# subset two samples and extract probe-level intensities
inds <- colnames(ab) == "1532a99hpp_av04.CEL.gz" | colnames(ab) == "1532b99hpp_av04.CEL.gz"
ab_sub <- ab[,inds]

# read probe-level intensities and subset
row <- probeNames(ab_sub) == "36085_at"
mat <- pm(ab_sub)[row,]
dim(mat)

# log2 ratios of intensities for each probe between two samples
pData(ab_sub)
expected <- log2(as.numeric(pData(ab_sub)[2,"X36085_at"])) - log2(as.numeric(pData(ab_sub)[1,"X36085_at"]))
observed <- log(mat[,2],2)-log(mat[,1],2)
boxplot(observed, ylim = c(-0.5, 1.5))
abline(h = expected, lty = 2, col = 'blue')
abline(h = 0)

## SummarizedExperiment
# RNA-seq data from untreated and steroid-treated airway smooth muscle cells
library(airway)
data(airway)
airway

# basic details of data
metadata(airway)
dim(airway)

# sample info
colData(airway)

# feature info
rowRanges(airway)
# each entry in a given row contains # of exons
rowRanges(airway)[100]
seqnames(rowRanges(airway)[100])
1-min(start(rowRanges(airway)[100])) + max(end(rowRanges(airway)[100]))

# analyzing a given gene expression across samples
sample <- airway["ENSG00000103196",]
mean(assay(sample))
dex <- sample[,which(colData(sample)$dex == "trt")]
mean(assay(dex))
nodex <- sample[,which(colData(sample)$dex == "untrt")]
mean(assay(nodex))
log2(mean(assay(dex))/mean(assay(nodex)))

## Importing NGS Data
library("Rsamtools")
library("GenomicAlignments")
# Rsamtools uses lower-level functions for NGS in standard formats
# both use BAM file formats
library("pasillaBamSubset")

# creating BamFile objects and extracting information on sequence length and summary
filename <- untreated1_chr4()
bf <- BamFile(filename)
seqinfo(bf)
sl <- seqlengths(bf)
quickBamFlagSummary(bf)

# what and which - reads from a particular genomic range
gr <- GRanges("chr4", IRanges(1, sl["chr4"]))
countBam(bf, param = ScanBamParam(which=gr))

# read 5 at a time to limit amount of data entered for memory purposes (or to debug)
reads <- scanBam(BamFile(filename, yieldSize = 5))

# analyzing the reads object
class(reads)
names(reads[[1]])
reads[[1]]$pos # start position
reads[[1]]$rname # chromosome
reads[[1]]$strand
reads[[1]]$qwidth
reads[[1]]$seq

# combine what and which to extract specific values over ranges
gr <- GRanges("chr4", IRanges(500000, 700000))
reads <- scanBam(bf, param = ScanBamParam(what = c("pos","strand"), which = gr))

# GenomicAlignments
ga <- readGAlignments(bf)
length(ga)
granges(ga[1])

# using GRanges functions
gr <- GRanges("chr4", IRanges(700000,800000))
fo <- findOverlaps(ga, gr)
countOverlaps(gr, ga)
table(ga %over% gr)

# Exercises
library(pasillaBamSubset)
library(Rsamtools)
library(Biostrings)
filename <- untreated1_chr4()

# count number of reads within range
gr <- GRanges("chr4", IRanges(440000, 470000))
bf <- BamFile(filename)
ga <- readGAlignments(bf)
countOverlaps(gr, ga)

# analyze sequence content - average GC content
reads <- scanBam(bf, param = ScanBamParam(what = "seq", which = gr))
reads <- reads$`chr4:440000-470000`$seq
match <- vcountPattern("C", reads) + vcountPattern("G", reads)
w <- width(reads)
mean(match/w)

# now use GenomicAlignments for analysis
library(GenomicAlignments)
ga <- readGAlignments(bf)
hist(start(ga), breaks = 100)

# load genes for Drosophila genome as GRanges
library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
g <- genes(TxDb.Dmelanogaster.UCSC.dm3.ensGene) 
g2 <- g[g %over% GRanges("chr4", IRanges(200000, 300000))]

# count overlaps for a given gene
gr <- g2["FBgn0039890", ]
countOverlaps(gr, ga)

## Count table creation
# load data and limit genes to chr4
library(pasillaBamSubset)
library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
txdb <- TxDb.Dmelanogaster.UCSC.dm3.ensGene
g <- genes(txdb)
g <- g[seqnames(g) == "chr4"]

# create list of exons for each gene
grl <- exonsBy(txdb, by="gene")
grl <- grl[names(g)]
all.equal(names(g), names(grl))

# get Bam file
library(Rsamtools)
bf <- BamFile(untreated1_chr4())
library(GenomicAlignments)

# summarize overlaps without considering strand information
so1 <- summarizeOverlaps(g, bf, ignore.strand = TRUE)
so2 <- summarizeOverlaps(grl, bf, ignore.strand = TRUE)

# compare results
plot(assay(so1) + 1, assay(so2) + 1, log="xy")
abline(0, 1)

# mean ratio of non-zero counts
a1 <- assay(so1)
a2 <- assay(so2)
mean(a2[a1 > 0]/a1[a1 > 0])
# GRangesList object shows lower counts due to excluding introns
# extra reads from introns may reflect precursor mRNA or other source (noisier data)

# proportion of reads aligned to each read per million reads
# i.e. RPM or FPM (reads/fragments per million)
FPM <- 1e6 * a2/sum(a2)
FPM[1]

# count should be proportional (roughly) to length of exons
widths <- sum(width(reduce(grl)))
summary(widths)

# fragments per kilobase of exonic basepairs per million mapped reads
FPKM <- 1000 * FPM/widths
FPKM[1]

## Public data repositories
library("hgu133plus2.db")
IDs <- select(hgu133plus2.db, keytype="SYMBOL", 
       columns=c("PROBEID"), keys="EGFR")$PROBEID
unique(IDs)

# using GO.db to find a tag for a biological process
library("GO.db")
key <- select(GO.db, keytype = "TERM", columns = "GOID", keys = "glial cell proliferation")
IDs <- select(hgu133plus2.db, keytype = "GO", columns = "PROBEID", keys = key$GOID)$PROBEID
unique(IDs)
