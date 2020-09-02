## Relational Databases: SQLite
library(RSQLite)
library(GO.db)
con <- GO.db$conn

# issue SQL queries using connection
dbGetQuery(con, "select count(*) from go_term")

# benchmarking SQL vs. select()
library(microbenchmark)
m1 <- microbenchmark(
  dbGetQuery(GO.db$conn, "select term from go_term"), times = 10L, unit="ms"
)
m2 <- microbenchmark(
  AnnotationDbi::keys(GO.db, keytype = "TERM"), times=10L, unit="ms")
median(m2$time)/median(m1$time)

## HDF5 and GenomicFiles
# HDF5 with SummarizedExperiment
library(HDF5Array)
library(airway)
data("airway")
td <- tempfile()
saveHDF5SummarizedExperiment(airway, td)

# files in directory
dir(td, full=TRUE)

# load RDS object
X <- readRDS(dir(td, full=TRUE)[2])
class(X)
X
assay(X[1,1])
# need to use loadHDF5SummarizedExperiment() to access assay()

# Find exons of RNA-seq HNRNPC knockdown experiment
library(erma)
hn <- genemodel("HNRNPC")
e1 <- hn[1] # first exon

# load dataset using GenomicFiles interface
library(GenomicFiles)
library(RNAseqData.HNRNPC.bam.chr14)
gf <- GenomicFiles(files=RNAseqData.HNRNPC.bam.chr14_BAMFILES)
rowRanges(gf) = e1
gf


# Map alignments overlapping region of interest
library(GenomicAlignments)
MAP <- function(r, f) readGAlignmentPairs(f, param=ScanBamParam(which = r))
ali <- reduceByRange(gf, MAP=MAP)

elementNROWS(ali[[1]])

# visualize regulatory states of chromatin near TSS of CD28
ermaset = makeErmaSet()
stateProfile(ermaset[,c(4,6,30,31)], "CD28", short=FALSE)

## Tabix-indexed genomic data
# Tabix-indexed BAM: regulatory and transcribed sites
library(RNAseqData.HNRNPC.bam.chr14)
library(GenomicAlignments)
library(ERBS)
data("GM12878")
seqlevels(GM12878, pruning.mode="coarse") = "chr14"
library(Rsamtools)
param = ScanBamParam(which=GM12878)
tab <- summarizeOverlaps(GM12878, RNAseqData.HNRNPC.bam.chr14_BAMFILES, param=param)

# how many ESRRA binding peaks are not subject to transcription in any cells
sum(rowSums(assay(tab))==0)

# look at 5th region of chr14
mm <- ScanBamParam(which=rowRanges(tab)[5], what="mapq")
bf <- RNAseqData.HNRNPC.bam.chr14_BAMFILES
kk <- scanBam(bf[1], param=mm)

# how many reads aligned here
length(kk$`chr14:93552286-93553668`$mapq)

# mean quality score
mean(kk$`chr14:93552286-93553668`$mapq)

# repeat for 30th interval
mm <- ScanBamParam(which=rowRanges(tab)[30], what="mapq")
bf <- RNAseqData.HNRNPC.bam.chr14_BAMFILES
kk <- scanBam(bf[1], param=mm)
length(kk$`chr14:21560446-21561239`$mapq)
mean(kk$`chr14:21560446-21561239`$mapq)
