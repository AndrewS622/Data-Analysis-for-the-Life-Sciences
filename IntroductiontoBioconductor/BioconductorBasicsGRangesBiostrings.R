## IRanges and GRanges
# genomic ranges are reported ESRRA binding sites from ENCODE ChIP-seq on HepG2 and GM12878 cell lines
library(ERBS)
data("HepG2")
class(HepG2)
length(HepG2$signalValue)

# exploring the data set - median signal
median(HepG2$signalValue)

# which chromosome has the largest signal
chr <- as.character(seqnames(HepG2))
chr[which.max(HepG2$signalValue)]

# how many regions from chr16
sum(chr == "chr16")

# histogram of region widths
hist(width(HepG2))
median(width(HepG2))

# IRanges
library(IRanges)
ir <- IRanges(101, 200)

# multiplication zooms in
ir * 2

# narrow shortens by specified parameter
narrow(ir, start = 20)

# + expands width
ir + 25

# declaring multiple IRanges
ir <- IRanges(c(1, 11, 21), c(3, 15, 27))
sum(width(ir))

# plotting ranges
start <- c(101,106,201,211,221,301,306,311,351,361,401,411,501)
end <- c(150,160,210,270,225,310,310,330,390,380,415,470,510)
x <- IRanges(start, end)
library(ph525x)
plotRanges(x)

# gaps in ranges
gaps(x)
sum(width(gaps(x)))

# disjoint ranges
disjoin(x)

# resize(x,1) gives only the starting points of each range
par(mfrow=c(2,1))
plotRanges(x, xlim=c(0,600))
plotRanges(resize(x,1), xlim=c(0,600))

# GRanges
x <- GRanges("chr1", IRanges(c(1,101), c(50,150)), strand=c("+","-"))

# obtain the IRanges from a GRanges object
ranges(x)

# resize takes strand into account
plotGRanges = function(x) plotRanges(ranges(x))
par(mfrow=c(2,1))
plotGRanges(x)
plotGRanges(resize(x,1))

# example partially overlapping ranges
x = GRanges("chr1", IRanges(c(101,201,401,501),c(150,250,450,550)), strand="+")
y = GRanges("chr1", IRanges(c(101,221,301,401,541),c(150,250,350,470,550)), strand="+")
plotGRanges(x)
plotGRanges(y)

# can combine multiple GRanges into one list object
GRangesList(x,y)

# or into a single GRanges object
c(x,y)

# total width covered by x and y
both <- disjoin(c(x,y))
x %over% y
sum(width(intersect(x,y)))

# width in just x or just y
sum(width(x)) + sum(width(y)) - 2 * sum(width(intersect(x,y)))

# opposite strand: %over% is strand-sensitive
z <- GRanges("chr1", range(ranges(x)), strand="-")
x %over% z

## Using GRanges for genomic analysis
library(ERBS)
data("HepG2")
data("GM12878")
start(HepG2)[17]

# finding closest regions to given locations
d <- distanceToNearest(HepG2, GM12878)
GMloc <- subjectHits(d)[17]
start(GM12878)[GMloc]
values(d)$distance[17]

# proportion of distances < 2000 bp
mean(values(d)$distance < 2000)
mean(mcols(d)$distance < 2000)

# Genes as GRanges
par(mfrow=c(1,1))
library(Homo.sapiens)
ghs <- genes(Homo.sapiens)
ghs

# genes represented
length(ghs$GENEID)

# chromosome with the most genes
sort(table(seqnames(ghs)),decreasing=TRUE)

# widths
h <- hist(width(ghs))
h <- hist(width(ghs), nc=1000, xlim=c(0,1000000))

median(width(ghs))

# finding closest genes
library(ERBS)
data("HepG2")
data("GM12878")
res <- findOverlaps(HepG2,GM12878)
erbs <- HepG2[queryHits(res)]
erbs <- granges(erbs)
erbs2 <- intersect(HepG2, GM12878)

# exploring the differences between these methods
identical(erbs,erbs2)
mean(erbs %in% erbs2)
sum(width(erbs))
sum(width(erbs2))

# transcription start sites
library(Homo.sapiens)
ghs <- genes(Homo.sapiens)
tss <- resize(ghs, 1)
start(tss["100113402"])

# gene with TSS nearest a binding site
gene <- nearest(erbs, tss)
ghs$GENEID[gene[4]]
keys = as.character(values(ghs[gene])$GENEID)
res <- select(Homo.sapiens, keys=keys, columns = "SYMBOL", keytype="GENEID")
res[4,]

# Getting sequences
library(ERBS)
library(Biostrings)
library(GenomicRanges)
data("HepG2")
data("GM12878")
res <- findOverlaps(HepG2, GM12878)
erbs <- HepG2[queryHits(res)]
erbs <- granges(erbs)

# human genome data
library(BSgenome.Hsapiens.UCSC.hg19)

# extract sequence of each region in erbs
seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg19, erbs)

# GC content
f <- letterFrequency(seqs, "GC")
w <-width(seqs)
GCcontent <- f/w
median(GCcontent)

# create control set of regions
ctrl <- shift(erbs, 10000)
seqsctrl <- getSeq(BSgenome.Hsapiens.UCSC.hg19, ctrl)
fctrl <- letterFrequency(seqsctrl, "GC")
wctrl <-width(seqsctrl)
GCcontentctrl <- fctrl/wctrl
median(GCcontentctrl)
