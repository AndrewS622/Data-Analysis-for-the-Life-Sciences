## Visualizing NGS
library(rtracklayer)

# after downloading ESRRA binding events bigWig
h2bw <- import("wgEncodeSydhTfbsHepg2ErraForsklnStdSig.bigWig")
h2bw
# score column contains coverage info
scores <- score(h2bw)

# median width of ranges
median(width(h2bw))

# compare to reported peaks from ERBS
library(ERBS) 
data(HepG2)

# overlapping regions
fo = findOverlaps(h2bw, HepG2)
inpeak = queryHits(fo)

# median score for regions identified as peaks in HepG2
median(scores[inpeak])

# median score for other regions
median(scores[-inpeak])

# find ESRRA gene 
library(Homo.sapiens)
esrra <- select(Homo.sapiens, keytype = "SYMBOL", 
                keys = "ESRRA",
                columns = c("TXSTART", "TXCHROM"))
esrra <- esrra[which.min(esrra$TXSTART),]

# find which peak in HepG2 overlaps with TXSTART
esrra_range <- makeGRangesFromDataFrame(esrra, start.field = "TXSTART", end.field = "TXSTART", seqnames.field = "TXCHROM")
esrra_peak <- HepG2[which(HepG2 %over% esrra_range),]

# find scores of regions in h2bw which overlap with this peak
overlaps <- findOverlaps(h2bw, esrra_peak)
inpeak <- queryHits(overlaps)
max(scores[inpeak])
peakcov <- h2bw[queryHits(fo)[subjectHits(fo) == 5]]

# plot coverage vs. middle point of each range
plot(0.5 * (start(peakcov) + end(peakcov)), peakcov$score)

# basepair with highest coverage
ranges(peakcov[which.max(peakcov$score),])

# Sushi package targets development of multipanel figures for genomics
library(Sushi)
data("Sushi_ChIPSeq_severalfactors.bed")
?Sushi_ChIPSeq_severalfactors.bed

# example from plotBed() man page
chrom = "chr15"
chromstart = 72800000
chromend = 73100000
plotBed(beddata = Sushi_ChIPSeq_severalfactors.bed,
        chrom = chrom,chromstart = chromstart,chromend =chromend,
        rownumber  = Sushi_ChIPSeq_severalfactors.bed$row, 
        type = "region", color=Sushi_ChIPSeq_severalfactors.bed$color, row="given", plotbg="grey95", 
        rowlabels=unique(Sushi_ChIPSeq_severalfactors.bed$name),
        rowlabelcol=unique(Sushi_ChIPSeq_severalfactors.bed$color),rowlabelcex=0.75)

# converting to GRanges
ss <- Sushi_ChIPSeq_severalfactors.bed
library(GenomicRanges)
k562gr <- GRanges(ss$chrom, IRanges(ss$start, ss$end))
k562grl <- split(k562gr, ss$name)

# how many distinct regions of genome
reduce(k562gr)
length(reduce(k562gr))

# plot binding display using Gviz
library(Gviz)
plotTracks(lapply(k562grl, AnnotationTrack))
plotTracks(lapply(1:length(k562grl), function(i) AnnotationTrack(reduce(k562grl[[i]]), name=names(k562grl)[i])))

## ggbio
# visualize GM12878 ESRRA peaks
library(ggbio)
library(GenomeInfoDb)
data(GM12878)
seqlevels(GM12878) = seqlevels(GM12878)[1:24]
autoplot(GM12878, layout="karyogram", aes(color=log(peak)))

# combine info on two cell lines
HepG2$cell = "HepG2"
GM12878$cell = "Bcell"
tot = c(GM12878, HepG2)
tot$peak10 = tot$peak/10 # for y-axis scale
seqlevels(tot, pruning.mode="coarse") = paste0("chr", 1:22)
library(scales)

# plot
p <- autoplot(seqinfo(tot))
p <- p + layout_karyogram(tot, aes(fill = cell, color = cell), geom="rect") + 
  scale_color_manual(values = alpha(c("green", "red"), 0.1)) + 
  scale_fill_manual(values = alpha(c("green", "red"), 0.1)) 
p + layout_karyogram(tot, aes(x = start, y = peak10), ylim=c(15,30), 
                    geom="point", color="blue", size=0.8)

# compute binding sites per chromosome
stot <- split(tot, as.character(seqnames(tot)))
w <- sapply(stot, function(x) sum(width(x)))
sort(w/seqlengths(tot)[names(w)], decreasing = TRUE)

# multitrack visualization
showz = function (sym = "ESRRA", radius = 1e+05)
{
  require(ggbio)
  require(erma)
  require(ERBS)
  es <- genemodel(sym)
  data(HepG2)
  data(GM12878)
  hsub <- subsetByOverlaps(HepG2, es+radius)
  gsub <- subsetByOverlaps(GM12878, es+radius)
  tracks(gene = GRangesList(es), hepNarr = autoplot(hsub), gmNarr = autoplot(gsub), 
         title = sym)
}

# use default options
p <- showz()
p

# zooming
p+zoom(2)

# debugging
debug(showz)
showz("OCM", radius = 2e6)
undebug(showz)

# Gviz
library(ph525x)
esrraScan

# Interactive visualization with shiny
library(Biobase)
library(hgu133a.db)

# hierarchical clustering function for expression set
esHclust = function(es) {
  emat <- t(exprs(es))
  rownames(emat) = sampleNames(es)
  dd <- data.frame(emat)
  dfHclust(dd)
}

# load in tissue expression data
library(tissuesGeneExpression)
data(tissuesGeneExpression)
tgeES <- ExpressionSet(e)
annotation(tgeES) = "hgu133a.db"
pData(tgeES) = tab
featureNames(tgeES) = 
  make.names(mapIds(hgu133a.db, keys=featureNames(tgeES), 
                    keytype="PROBEID", column="SYMBOL"), unique=TRUE)
sampleNames(tgeES) = make.names(tgeES$Tissue, unique=TRUE)
sum(is.na(rownames(tgeES)))

# plot
esHclust(tgeES[1:50,1:50])

# identify significant genes for targeting
library(limma)
mm <- model.matrix(~Tissue, data=pData(tgeES))
f1 <- lmFit(tgeES, mm)
ef1 <- eBayes(f1)
sig50 <- rownames(ef1$coefficients[order(ef1$F, decreasing = TRUE)[1:50],])
esHclust(tgeES[sig50,])

# using ML to assess set of 5 genes with greatest F-values for tissue discrimination
sig5 <- c("IL26", "ZNF674", "UBC.1", "C7orf25.1", "RPS13")
library(randomForest)
library(MLInterfaces)
set.seed(1234)
rf2 <- MLearn(Tissue~., tgeES[sig5,], randomForestI, xvalSpec("NOTEST"))
RObject(rf2)
