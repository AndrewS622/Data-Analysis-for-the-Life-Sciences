## ENCODE ATAC-seq data
library(AnnotationHub)
ah = AnnotationHub()
AnnotationHub::query(ah, "ENCODExplorerData")

fm = ah[["AH75132"]]

#filter to Homo sapiens, ATAC-seq
hfm = fm[which(fm$organism == "Homo sapiens"),]
tail(sort(table(hfm$assay)))
hfma = hfm[which(hfm$assay == "ATAC-seq"),]

# filter to A549 lung cancer BED files, keep only relevant columns
library(dplyr)
fhfma = as_tibble(hfma) %>%
  dplyr::select(biosample_name, output_type, file_format,
                treatment, treatment_duration, treatment_duration_unit, 
                cloud_metadata.url, biological_replicate_number) %>%
  dplyr::filter(biosample_name == "A549",
                file_format == "bed")

# make summary table of sample treatment status
fhfma %>%
  dplyr::count(biosample_name, treatment, 
               treatment_duration, treatment_duration_unit)

# extract URL filenames
library(GenomicFiles)
dexgf = GenomicFiles(files=fhfma$cloud_metadata.url)
colData(dexgf) = DataFrame(fhfma)

# Downloading ATAC-seq data and modeling treatment response to dexamethasone
# use Ensembl v79 genes
library(EnsDb.Hsapiens.v79)
g79 = genes(EnsDb.Hsapiens.v79)

# promoters returns interval 2000 bp upstream to 200 bp downstream of TSS
# change style of seqname to UCSC for BED file compatibility
p = promoters(g79)
GenomeInfoDb::seqlevelsStyle(p) = "UCSC"

# filter just to genes of interest to examine chromatin accessibility, remove unused seqnames
pint = p[p$symbol %in% c("VDR", "TREM1", "CEBPB")]
pint = GenomeInfoDb::keepStandardChromosomes(pint)
pint
rowRanges(dexgf) = pint

# run the query, specifying extra columns with cc since BED file format is non-standard
# request SummarizedExperiment with summarize = TRUE
rr = reduceByFile(dexgf, MAP=function(range, file) {
  cc = c(exta = "numeric", extb = "numeric", 
         extc = "numeric", extd = "numeric")
  cur = rtracklayer::import(file, extraCols=cc, which=range)
  cat(".")
  sum(GenomicRanges::width(cur))
}, summarize=TRUE)

# value = total width of accessible chromatin in promoter of genes of interest in each sample
assay(rr)
colData(rr) = colData(dexgf)

# replace treatment_time at untreated time point from NA to 0
rr$treatment_duration[which(is.na(rr$treatment_duration))] = 0

# model for relationship between dexamethasone exposure and chromatin accessibility
# approximate width of accessible regions as response of quadratic regression with duration as predictor
# fit only for second gene (row)
lm1 = lm(as.numeric(assay(rr[2, ])) ~ poly(rr$treatment_duration, 2))

# plot results
plot(as.numeric(assay(rr[2,])) ~ rr$treatment_duration,
     xlab="duration of dexamethasone exposure",
     ylab="avg. width accessible region")
lines(rr$treatment_duration[order(rr$treatment_duration)],
      predict(lm1)[order(rr$treatment_duration)])

# Assessment
pint$gene_name[2]

# which duration had the largest accessible width
rr$treatment_duration[which.max(assay(rr[2,]))]

# compute average score over intervals as well as total width
rr2 = reduceByFile(dexgf, MAP=function(range,file) {
  cc = c(exta = "numeric", extb = "numeric", 
         extc = "numeric", extd = "numeric")
  cur = rtracklayer::import(file, extraCols=cc, which=range)
  cat(".")
  c(swid = sum(GenomicRanges::width(cur)),
  msco = mean(cur$score, na.rm = TRUE))
}, summarize=TRUE)

# relationship between ATAC-seq scores and width of regions on second promoter
scos = apply(assay(rr2), 1:2, function(x) x[[1]][2])
ws = apply(assay(rr2), 1:2, function(x) x[[1]][1])
plot(scos["ENSG00000172216",],
     ws["ENSG00000172216",],
     xlab="average score over accessible region",
     ylab="accessible width over promoter")

# use linear model for score vs. width
df <- data.frame(scos = scos["ENSG00000172216",], ws = ws["ENSG00000172216",])
fit <- lm(ws ~ scos, data = df)
summary(fit)
# presence of outlier makes p large and model ineffective

# using robust regression to account for outlier presence
library(MASS)

# using M and MM methods (MM is more computationally intensive)
coefs_m <- summary(rlm(ws["ENSG00000172216",] ~ scos["ENSG00000172216",],
            method = "M"))
coefs_mm <- summary(rlm(ws["ENSG00000172216",] ~ scos["ENSG00000172216",],
            method = "MM"))
# thus, there is a significant inverse relationship

# with approximately Gaussian t-statistics, two-tailed p-value is:
p_m <- 2*pnorm(-abs(coefs_m$coef[2,"t value"]))
p_mm <- 2*pnorm(-abs(coefs_mm$coef[2,"t value"]))
-log10(p_m)
-log10(p_mm)

# incorporate score/width relationship into model for epigenetic response 
colData(rr2) = colData(dexgf)
rr2$treatment_duration[which(is.na(rr2$treatment_duration))] = 0
rr2df <- data.frame(
  wid = as.numeric(sapply(assay(rr2[2,]), "[", 1)),
  sco = as.numeric(sapply(assay(rr2[2,]), "[", 2)), 
  dur = rr2$treatment_duration)
lm1redo = lm(wid~poly(dur,2), data=rr2df)
lm1rev = lm(wid~poly(dur,2) + sco, data = rr2df)
summary(lm1redo)
summary(lm1rev)

# remove outlier
lm1rev2 <- lm(wid~poly(dur,2) + sco, data=rr2df, subset=-13)
summary(lm1rev2)

# Analyzing all promoter regions on chr17q
# q arm of chr17 associated with asthma
query(ah, c("cytoband", "hg38"))

# prepackaged promoter regions in following:
library(abEncoTools)
head(abEncoTools::p17q, 3)
data(p17q, package = "abEncoTools")
length(p17q)

# dexamethasone time course data at 100 promoter regions 
# FOR REFERENCE ONLY - DO NOT RUN
# library(GenomicFiles)
# library(abEncoTools)
# data(dexgf)
# data(p17q)
# rowRanges(dexgf) = p17q[1:100]
# library(BiocParallel)
# register(MulticoreParam(50))
# demo = reduceByRange(dexgf, MAP=function(range,file) {
#     cc = c(exta = "numeric", extb = "numeric",   # better for parallelization on Windows
#          extc = "numeric", extd = "numeric")     # to have this in the loop
#     cur = rtracklayer::import(file, extraCols=cc, which=range)
#     cat(".")
#     if (length(cur)==0) return(c(wid=NA, sco=NA))
#     c(wid=sum(width(cur)), msco=mean(cur$score, na.rm=TRUE))
#   }, summarize=TRUE)

# load results of above analysis
data("widScoSE")
widScoSE
which(rowData(widScoSE)$gene_name == "ORMDL3")

# fit quadratic mixed-effects models to genes 
# FOR REFERENCE ONLY - DO NOT RUN ON YOUR MACHINE
# library(parallel)
# options(mc.cores=50)
# library(abEncoTools)
# data(dexgf)
# len5 = length(widScoSE)
# mm = mclapply(1:len5, function(x) 
#      do.call(rbind, assay(widScoSE)[x,])) 
# ok = sapply(mm, function(x) !any(is.na(x[,1]))) # some genes have no accessibility recorded in their promoter regions
# mmm = mm[ok]  # we drop those genes
# dur = colData(dexgf)$treatment_duration
# dur[is.na(dur)] = 0
# now we assemble a list of data frames with columns width, score, replicate, and treatment duration
# mmmm = lapply(mmm, function(x) cbind(x, 
#     rep=colData(dexgf)$biological_replicate_number, dur=dur))
# library(nlme)
# NRpoly = mclapply(mmmm, function(x) summary(lme(wid~poly(dur,2), random=~1|rep, data=data.frame(x))))
# names(NRpoly) = make.names( rowData(widScoSE)$symbol[ok], unique=TRUE)

# load results
data("NRpoly")

# genes with accessible chromatin
length(NRpoly)

# lme fits
names(NRpoly)[1:5]
NRpoly[["WSB1"]]
anova(NRpoly[["WSB1"]])

# summarize p-values over all promoter regions with accessible chromatin
allap <- sapply(NRpoly, function(x) anova(x)[2,4])
head(sort(allap))

# dynamic response of accessibility
pp = function(ml) {
  d = matrix(unlist(ml), nr = 2)
  dur = rr2$treatment_duration
  plot(dur, d[1,])
}

# example
rowData(widScoSE)[147,"gene_name"]
pp(assay(widScoSE)[147,])

# computer FDRs with Benjamini-Hochberg
padj <- p.adjust(allap, method = "BH")
sum(padj < 0.05)
