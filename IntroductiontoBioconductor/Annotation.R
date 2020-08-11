## Annotation of Genes and Transcripts
# Reference genomes
library(BSgenome)
library(Biostrings)
ag = available.genomes()

# how many non-masked zebrafish genomes
zebrafish <- grep("Drerio", ag, value=TRUE)
masked <- grep(".masked", ag, value=TRUE)
length(zebrafish) - length(intersect(zebrafish, masked))

# masked genomes isolate ambiguous, low-complexity, or uninformative segments
library(BSgenome.Hsapiens.UCSC.hg19.masked)
c17m = BSgenome.Hsapiens.UCSC.hg19.masked$chr17
c17m
c22m <- BSgenome.Hsapiens.UCSC.hg19.masked$chr22
c22m

## Packages for gene and transcript catalogs tutorial
# import TxDb transcript database
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb = TxDb.Hsapiens.UCSC.hg19.knownGene  # abbreviate
class(txdb)
methods(class="TxDb")

# extract and inspect genes from TxDb
genes(txdb)
table(strand(genes(txdb)))
summary(width(genes(txdb)))

# inspect largest gene in genome
id = which.max(width(genes(txdb)))
genes(txdb)[id]
library(org.Hs.eg.db)
select(org.Hs.eg.db, keys="286297", keytype = "ENTREZID", columns = c("SYMBOL", "GENENAME"))

# compare total size of exons to total size of genes
ex = exons(txdb)
rex = reduce(ex)
ex_width = sum(width(rex))    # bases in exons
gene_width = sum(width(genes(txdb)))    # bases in genes
ex_width/gene_width

# Gene and transcript model
library(devtools)
install_github("genomicsclass/ph525x")
library(ph525x)
stopifnot(packageVersion("ph525x") >= "0.0.16") 

# use modPlot to visualize ESR1 pathway
modPlot("ESR1", useGeneSym=FALSE, collapse=FALSE) 

# number of transcripts comprising ESR1 model
library(ensembldb)
edb <- EnsDb.Hsapiens.v75
library(EnsDb.Hsapiens.v75)

# extract ID for ESR1
idnum <- select(edb, keys = "ESR1", keytype = "SYMBOL", columns = "ENTREZID")$ENTREZID
txs <- transcripts(txdb, filter = list(gene_id = idnum))
length(txs)

## AnnotationHUb
library(AnnotationHub)
ah <- AnnotationHub()
ah

length(unique(ah$species))

ah_human <- subset(ah, species == "Homo sapiens")
ah_human

query(ah, "HepG2")
query(ah, c("HepG2", "H3K4me3"))
hepg2 <- query(ah, "HepG2")
hepg2_h3k4me3 <- query(hepg2, c("H3K4me3"))
hepg2_h3k4me3
hepg2_h3k4me3$tags

display(query(ah, "HepG2"))

e118_broadpeak <- query(hepg2_h3k4me3, c("E118", "broadPeak"))
id <- e118_broadpeak$ah_id
id

hepg2_h3k4me3_broad <- ah[["AH29728"]]
hepg2_h3k4me3_broad
alt_format <- ah[[id]]
identical(hepg2_h3k4me3_broad, alt_format)

# AnnotationHub Exercises
library(AnnotationHub)
ah = AnnotationHub()
mah = mcols(ah)
names(mah)

# sort number of entries per species
sort(table(mah$species), decreasing=TRUE)[1:10]

# nested queries to select CTCF binding in HepG2 genome
names(query(query(ah, "HepG2"), "CTCF"))

## liftOver: translating between reference builds
# learn about liftOver from rtracklayer
library(rtracklayer)
?liftOver

# chromosome 1 gene locations in hg38
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
tx38 <- TxDb.Hsapiens.UCSC.hg38.knownGene
seqlevels(tx38, pruning.mode="coarse") = "chr1"
g1_38 <- genes(tx38)

# download the hg38 to hg19 chain file
library(AnnotationHub)
ah <- AnnotationHub()
ah.chain <- subset(ah, rdataclass == "ChainFile" & species == "Homo sapiens")
query(ah.chain, c("hg19", "hg38"))
ch <- ah[["AH14108"]]

# perform the liftOver
g1_19L <- liftOver(g1_38, ch)
g1_19L

# liftOver Exercises
# hg19 to hg38 lifover file
download.file("http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz", "hg19ToHg38.over.chain.gz")
library(R.utils)
gunzip("hg19ToHg38.over.chain.gz")

# convert HepG2 binding addresses to hg38
library(ERBS)
data(HepG2)
library(rtracklayer)
ch = import.chain("hg19ToHg38.over.chain") 
nHepG2 = liftOver(HepG2, ch)

# shifting of start site of first range
ss19 <- 20335378
ss19 - start(nHepG2)[[1]]

## rtracklayer for data import and export
library(devtools)
install_github("genomicsclass/ERBS")    # install ERBS package

f1 = dir(system.file("extdata",package="ERBS"), full=TRUE)[1]    # access data
readLines(f1, 4)    # look at a few lines

library(rtracklayer)
imp = import(f1, format="bedGraph")    # import as bedGraph format
imp

genome(imp)    # genome identifier tag not set, but you should set it
genome(imp) = "hg19"    # set genome
genome(imp)

export(imp, "demoex.bed")    # export as BED format  
cat(readLines("demoex.bed", n=5), sep="\n")    # check output file

# rtracklayer Exercises
library(rtracklayer)
data(targets)
class(targets)

# create GRanges instance from targets data frame
library(GenomicRanges)
mtar <- with(targets,
             GRanges(chrom, IRanges(start,end), strand=strand,
                     targets=target, mirname=name))

# look at exported versions of data
cat(export(mtar[1:5], format="bed"), sep="\n")
cat("\n")
cat(export(mtar[1:5], format="gff3"), sep="\n")

## Annotation tools for systems biology
# OrgDb
library(org.Hs.eg.db)
org <- org.Hs.eg.db

# counting genes in a cytoband
genes <- select(org, key="17q21.1", keytype = "MAP", columns = c("GO", "GENENAME"))
length(unique(genes$GENENAME))

# most common annotations on this cytoband
sort(table(genes$GO), decreasing = TRUE)[1:5]

# how many annotations for ORMDL3 have traceable author statement as evidence
genes <- select(org, keytype = "SYMBOL", keys = "ORMDL3", columns = c("EVIDENCE", "GO"))
sum(genes$EVIDENCE == "TAS")

# enumerating genes where ER binds at TSS
library(Homo.sapiens)
g = genes(Homo.sapiens)
library(ERBS)
data(HepG2)

# set of genes where TSS lies at HepG2 ER binding site
kp = g[resize(g,1) %over% HepG2]

# interactive HTML5 report on gene annotation
nn = names(kp)
m = select(Homo.sapiens, keys=nn, keytype="ENTREZID",
           columns=c("SYMBOL", "GENENAME", "TERM", "GO"))
library(DT)
datatable(m)

# KEGG
library(KEGGREST)
obj <- unlist(keggGet("hsa:3845"))
obj$NAME

# number of pathways
length(keggLink("pathway", "hsa:3845"))

# Folate Biosynthesis pathway
library(png)
oo = keggGet("hsa00790", "image")
writePNG(oo, "im.png")

# Ontology lookup using rols
library(rols)
diab <- OlsSearch("diabetes")
olsRows(allRows(diab))

# look at all results
fulld <- olsSearch(allRows(diab))
adf <- as(fulld, "data.frame")
sort(table(adf$ontology_name), decreasing = TRUE)[1:10]

# create data table
library(DT)
datatable(adf)