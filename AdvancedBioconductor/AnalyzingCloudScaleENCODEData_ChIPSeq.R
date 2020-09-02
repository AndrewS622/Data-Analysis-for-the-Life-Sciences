## ENCODE ChIP-seq data
# Retrieving ChIP-seq metadata
library(AnnotationHub)
ah = AnnotationHub()
AnnotationHub::query(ah, "ENCODExplorerData")

# retrieve most recent full resource from above results
fm = ah[["AH75132"]] # implicitly retrieves if needed, or loads
dim(fm)
table(fm$organism)

# filter to Homo sapiens and ChIP-seq
hfm = fm[which(fm$organism == "Homo sapiens"),]
tail(sort(table(hfm$assay)))
hfmt = hfm[which(hfm$assay == "TF ChIP-seq"),]

# proteins most frequently assayed
tail(sort(table(hfmt$antibody_target)))

# restrict to BigWig file type
hfmtbw = hfmt[hfmt$file_type == "bigWig",]
dim(hfmtbw)
table(hfmtbw$output_type)

# use GenomicFiles to manage URL references
library(GenomicFiles)
gf1 = GenomicFiles(files=hfmtbw$cloud_metadata.url)
colData(gf1) = as(as.data.frame(hfmtbw), "DataFrame")
gf1

# to access full URL
GenomicFiles::files(gf1)[1]

# human genome assembly used, restrict to 38
table(gf1$assembly)

gf1 = gf1[, which(gf1$assembly == "GRCh38")]
gf1

# Exploring metadata
tls = function(x) t(t(head(sort(table(x), decreasing=TRUE))))

# sample cell lines and ontology
tls(gf1$biosample_name)
tls(gf1$biosample_ontology)

# library to create hierarchy of ontological labels
library(ontoProc)

# EFO = experimental factor ontology tags
ee = getEFOOnto()
ioi = tail(sort(table(hfmt$biosample_ontology)),12)

# simplify format of tags
nn = names(ioi[grep("EFO", names(ioi))])
tails = gsub(".*_", "", nn)
tailss = gsub("/", "", tails)
tt = paste0("EFO:", tailss)

# plot showing which tags map to which cell lines
onto_plot2(ee, tt, cex=.6)

# Downloading ChIP-seq data and sketching peak scores
# restrict to CREB1
gf1_creb1 = gf1[, which(gf1$target == "CREB1")]
table(gf1_creb1$treatment, gf1_creb1$biosample_name, exclude=NULL)

# available output types for cell lines; look at signal p-value
table(gf1_creb1$output_type, gf1_creb1$biosample_name, exclude=NULL)
gf1_creb1_use = gf1_creb1[, which(gf1_creb1$output_type == "signal p-value")]

# download data in short region on chr17
myr = GRanges("chr17", IRanges(66e6,67e6))
genome(myr) = "GRCh38"
rowRanges(gf1_creb1_use) = myr

# data provided in R package for Windows since BigWig functions do not work on this OS
library(rtracklayer)
if (.Platform$OS.type == "windows") {
  sels = edxAdvBioc::creb1_sels
} else {
  sels = reduceByFile(gf1_creb1_use, MAP=function(range, file, ...) {
    sel = rtracklayer::BigWigSelection(range)
    rtracklayer::import.bw(file, selection=sel, genome=genome(myr)[1])
  })
}

# pull out signal p-values for first two samples
lk1 = sels[[1]][[1]]
lk2 = sels[[2]][[1]]

# plot CREB1 peal -log10(p-values)
plot(start(lk1)+.5*(width(lk1)), 
     lk1$score, pch=19,
     xlab="midpoint of scored interval, chr17", 
     ylab="-log10 CREB1 signal p-value",
     cex=.5, 
     col=lava::Col("black", alpha=.3))
points(start(lk2)+.5*width(lk2), 
       lk2$score, 
       pch=19, 
       col=lava::Col("blue", alpha=.3), cex=.5)

# Assessment
# file size for bigWig files
median(as.numeric(colData(gf1)$cloud_metadata.file_size))
sum(as.numeric(colData(gf1)$cloud_metadata.file_size))/(1e9)

# most common target protein
tls(gf1$target)

# most commonly studied biological material
tls(gf1$biosample_name)

# filter metadata to ATAC-seq in Homo sapiens
hfma = hfm[which(hfm$assay == "ATAC-seq"),]
tls(hfma$biosample_name)

# timeline of assay types
library(dplyr)
library(lubridate)
newt = as_tibble(hfm) %>%
  dplyr::select(assay, date_created) %>%
  dplyr::mutate(ldate = lubridate::as_date(date_created))

par(las = 2, mar = c(12,4,2,2))
ameds = sapply(split(newt$ldate, newt$assay), median)
with(newt, boxplot(split(ldate, assay)[order(ameds)]))

# how many of most recent type
sum(newt$assay == "long read RNA-seq", na.rm = TRUE)

# look at which genes knocked down with CRISPRi
ii = as_tibble(hfm) %>%
  dplyr::filter(assay == "CRISPRi RNA-seq", 
                file_format == "tsv") %>%
  dplyr::mutate(targ = gsub("RNA-seq on K562 cells treated by CRISPR interference targeting (..*).", 
                            "\\1", dataset_description)) %>%
  dplyr::distinct(investigated_as, targ)
ii
sum(ii$investigated_as == "transcription factor")

# import bigWig files in specific range and output type
import_enc_bw = function(gf, target = "CREB1", 
                         biosample_names = c("HepG2", "MCF-7", "A549"),
                         output_type = "signal p-value", 
                         selection = GRanges("chr17", IRanges(66e6, 67e6), 
                                             genom = "GRCh38")) {
  stopifnot(length(selection) == 1)
  gf_use = gf[, which(gf$target == target)]
  gf_use = gf_use[, which(gf_use$output_type == output_type)]
  gf_use = gf_use[, which(gf_use$biosample_name %in% biosample_names)]
  rowRanges(gf_use) = selection
  ans = reduceByFile(gf_use, MAP=function(range, file, ...){
    sel = BigWigSelection(range)
    import.bw(file, selection=sel, genome=genome(myr)[1])
  })
  ans = GRangesList(unlist(ans, recursive = FALSE))
  mcols(ans) = colData(gf_use)
  names(ans) = make.unique(mcols(ans)$biosample_name)
  ans
}

# plot log(p-values) 
plot_pair = function(impeb, ylim=c(1,500), 
                     xlim=c(66.1e6,66.3e6),
                     logy=TRUE, leg_frac_x=0.015) {
  stopifnot(length(impeb) == 2)
  lk1 = impeb[[1]]
  lk2 = impeb[[2]]
  if (logy) {
    lk1 = lk1[which(lk1$score>0)]
    lk2 = lk2[which(lk2$score>0)]
  }
  seqn = as.character(seqnames(impeb[[1]]))[1]
  g = grep("chr", seqn)
  if (length(g)==0) seqn=paste0("chr", seqn)
  targ = mcols(impeb)$target[1]
  oty = mcols(impeb)$output_type[1]
  pspec = list(x=start(lk1)+.5*(width(lk1)),
               y=lk1$score, pch=19, 
               xlab = sprintf("midpoint of scored interval, %s", seqn),
               ylab = sprintf("bigWig '%s' for %s", oty, targ),
               cex=.5, col=lava::Col("orange", alpha=.1),
               ylim=ylim, xlim=xlim, log=ifelse(logy,"y",""))
  do.call(plot, pspec)
  points(start(lk2)+.5*width(lk2), lk2$score, pch=19,
         col=lava::Col("blue", alpha=.1), cex=.5)
  legend(xlim[1]+leg_frac_x*diff(xlim), ylim[2], pch=19,
         col=lava::Col(c("orange", "blue"), alpha=.4),
         legend=mcols(impeb)$biosample_name)
}

# download data (bigWig data cannot be imported on Windows due to reduceByFile function)
if (.Platform$OS.type == "windows") {
  skd = edxAdvBioc::skd
} else {
  skd = import_enc_bw(gf1)
}
table(mcols(skd)$biosample_name)

# select pairs of samples from different cell types and visualize
par(mfrow=c(2,1))
plot_pair(skd[c("HepG2.2", "MCF-7")], ylim=c(1,1500), xlim=c(66e6, 67e6), leg_frac_x=0.04)
plot_pair(skd[c("HepG2.1", "MCF-7.1")], ylim=c(1,1500), xlim=c(66e6, 67e6), leg_frac_x=0.04)

# simplify data on one HepG2 sample
s5 = skd[["HepG2.2"]]
GenomeInfoDb::seqlevelsStyle(s5) = "NCBI"
seqlevels(s5) = "17"

# find gene near strongest peak for CREB1
gr = s5[which.max(s5$score),]
library(EnsDb.Hsapiens.v79)
genes = genes(EnsDb.Hsapiens.v79)
ind = subjectHits(distanceToNearest(gr, genes))
genes[ind,"gene_name"]

# FOXA1 target in breast cancer cells
if (.Platform$OS.type == "windows") {
  uu = edxAdvBioc::uu
} else {
  uu = import_enc_bw(gf1, target="FOXA1", 
                     selection=GRanges("chr6", IRanges(151.5e6,152.5e6)), 
                     biosample_names=c("A549", "MCF-7"))
}
par(mfrow=c(1,1))
plot_pair(uu[c("A549", "MCF-7")], xlim=c(151.5e6,152.5e6))

# what gene is closest to strongest peak in MCF-7
uu = uu[[2]]
genome(uu) = "GRCh38"
seqlevels(uu) = "chr6"
gr = uu[which.max(uu$score)]
seqlevels(genes, pruning.mode = "coarse") = "6"
gran = ranges(gr)
genes_ran = ranges(genes)
ind = subjectHits(distanceToNearest(gran, genes_ran))
genes[ind,]