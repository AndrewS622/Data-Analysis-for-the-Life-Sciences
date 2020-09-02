## ENCODE CRISPRi RNA-seq
library(AnnotationHub)
ah = AnnotationHub()
AnnotationHub::query(ah, "ENCODExplorerData")

# retrieve and filter
fm = ah[["AH75132"]]
hfm = fm[which(fm$organism == "Homo sapiens"),]
tail(sort(table(hfm$assay)))
hfmc = hfm[which(hfm$assay %in% c("CRISPR RNA-seq", "CRISPRi RNA-seq")),]
hfmc[1:4, 1:5]

# use gene quantifications with GRCh38 assembly
gq = hfmc[hfmc$output_type == "gene quantifications",]
table(gq$assembly)
gq = gq[which(gq$assembly == "GRCh38"),]

# store URLs and add metadata table
library(GenomicFiles)
gf2 = GenomicFiles(files = gq$cloud_metadata.url)
colData(gf2) = as(gq, "DataFrame")
tail(table(colData(gf2)$target))

# look at first few lines of first file
poke = read.delim(files(gf2)[1], sep = "\t", nrow=3)
poke

# keep only gene_id and expected count
col_cl = c(gene_id = "character", 
           `transcript_id(s)` = "NULL",
           length = "NULL", 
           effective_length = "NULL", 
           expected_count = "numeric",
           TPM = "NULL", 
           FPKM = "NULL",
           posterior_mean_count = "NULL",
           posterior_standard_deviation_of_count = "NULL",
           pme_TPM = "NULL",
           pme_FPKM = "NULL",
           TPM_ci_lower_bound = "NULL",
           TPM_ci_upper_bound = "NULL",
           FPKM_ci_lower_bound = "NULL",
           FPKM_ci_upper_bound = "NULL")

# helper function to read in GenomicFiles
get_EC = function(range, file, colClasses=NULL, ...) {
  x = as.data.frame(data.table::fread(file, colClasses = colClasses, 
                                      showProgress = FALSE))
  rownames(x) = x[,1]
  x = x[, 2, drop = FALSE]
  colnames(x) = gsub(".tsv", "", basename(file))
  x
}

# read in first 6 files
library(BiocParallel)
register(SerialParam())
rowRanges(gf2) = GenomicRanges::GRanges("chr0", IRanges(1, 2))
lk6 = reduceByFile(gf2[,1:6], MAP = function(range, file, ...){
  get_EC(range, file, colClasses = col_cl)
})

# combine results into single data frame where column anames are sample ids
fulldf = do.call(cbind, unlist(lk6, recursive = FALSE))

# load data from full 266-file data frame
data(criEC, package = "abEncoTools")

# to load this data manually (on multi-core machine)
# DESIGNED TO BE RUN ON A MULTI-CORE MACHINE - do not run on a personal computer
# lkall = reduceByFile(gf2, MAP=function(range, file, ...)
#   get_EC(range, file, colClasses=col_cl), BPPARAM=MulticoreParam(50))
# fulldf = do.call(cbind, unlist(lkall, recursive=FALSE))

# Exploring CRISPRi RNA-seq data
table(criEC$biosample_name)

# unique RNA binding proteins
length(unique(colData(criEC)$target[colData(criEC)$investigated_as == "RNA binding protein"]))
unique(colData(criEC)$target[colData(criEC)$investigated_as == "RNA binding protein"])

# most commonly targeted gene
sort(table(colData(criEC)$target), decreasing = TRUE)

# unique chromatin remodeler targets
length(unique(colData(criEC)$target[colData(criEC)$investigated_as == "chromatin remodeler"]))
unique(colData(criEC)$target[colData(criEC)$investigated_as == "chromatin remodeler"])

# target genes present among assay features
idx <- which(rowData(criEC)$gene_name %in% colData(criEC)$target)
length(idx)
rowData(criEC)$gene_name[idx]
length(unique(rowData(criEC)$gene_name[idx]))

# Examining CRISPRi Effects
sind = which(rowRanges(criEC)$gene_name == "STAT1")
contrids = grep("control", criEC$target)
expids = which(criEC$target == "STAT1")
limSE = criEC[sind, c(contrids, expids)]
boxplot(split(as.numeric(assay(limSE)), limSE$target), ylab = "STAT1 expected count")

# function to do the above for any gene
cribox = function(target, SE) {
  sind = which(rowRanges(SE)$gene_name == target)
  contrids = grep("control", SE$target)
  expids = which(criEC$target == target)
  limSE = SE[sind, c(contrids, expids)]
  boxplot(split(as.numeric(assay(limSE)), limSE$target), ylab = paste(target, "expected count"))
}

cribox("STAT5A", criEC)

# sources of variation in overall expression
# subset to genes expressed w/ expected count > 1 in all samples
EC = assay(criEC)
bas = which(apply(EC,1,min)>1)
ECl = EC[bas,]
dim(ECl)

# compute principal components and look at some potential batch effects
prc = prcomp(t(log(ECl)))
table(criEC$submitted_by)
pairs(prc$x[,1:3], pch=19, col=factor(criEC$submitted_by))

head(criEC$date_created)
datelen = 10
table(substr(criEC$date_created, 1, datelen))

# new boxplot function stratifying by submitter ID and day of submission, only w/ controls
cribox2 = function(target, SE, datelen=10) {
  sind = which(rowRanges(SE)$gene_name == target)
  contrids = grep("control", SE$target)
  limSE = SE[sind, contrids]
  fac = paste0(substr(limSE$date_created, 1, datelen),
               substr(limSE$submitted_by, 1, 3))
  par(mar = c(7,4,2,2))
  boxplot(split(as.numeric(assay(limSE)), fac),
          ylab = paste(target, "expected count"),
          las = 2)
}

# look at some genes
cribox2("STAT1", criEC)
cribox2("STAT5A", criEC)
cribox2("TP53", criEC)

# finer stratification by adjusting datelen
cribox2("TP53", criEC, datelen=13)

# Adjusting for batch effects with a linear model
# function to fit linear model controlling for batch effects
crilm = function(target, SE, dateprec=10) {
  sind = which(rowRanges(SE)$gene_name == target)
  contrids = grep("control", SE$target)
  limSE = SE[sind, c(contrids, which(SE$target == target))]
  fac = paste0(substr(limSE$date_created,1,dateprec), substr(limSE$submitted_by,1,3))
  istarg = rep(0,ncol(limSE))
  istarg[which(limSE$target == target)] = 1
  ndf = data.frame(EC = as.numeric(assay(limSE)), date_conc_submitter = fac, istarg = istarg)
  lm(log(EC) ~ date_conc_submitter + istarg, data = ndf)
}

# look at two genes
m1 = crilm("STAT1", criEC)
summary(m1)
plot(m1,2)

m2 = crilm("STAT5A", criEC)
summary(m2)
plot(m2, 2)

# use non-parametric assessment of adjusted effect with permutation null distribution
crilm_shuff = function(target, SE, dateprec=10) {
  sind = which(rowRanges(SE)$gene_name == target)
  contrids = grep("ontrol", SE$target)
  limSE = SE[sind, c(contrids, which(SE$target==target)) ]
  fac = paste0(substr(limSE$date_created,1,dateprec), 
               substr(limSE$submitted_by,1,3))
  istarg = rep(0,ncol(limSE))
  istarg[which(limSE$target==target)]=1
  
  nel = length(istarg)
  while(1) {
    shuff = sample(seq_len(nel), size = nel, replace = FALSE)
    ndf = data.frame(EC = as.numeric(assay(limSE)),
                     date_conc_submitter = fac,
                     istarg = istarg[shuff])
    ans = lm(log(EC) ~ date_conc_submitter + istarg, data = ndf)
    if (!is.na(coef(ans)["istarg"]))
      break
  }
  ans
}

# compute 1000 permutations for STAT1 
set.seed(1234)
stat1null = vapply(1:1000, function(x) {
  crilm_shuff("STAT1", criEC)$coef["istarg"]
}, FUN.VALUE = numeric(1))

# plot STAT1 effect vs. null
obsstat1 = coef(m1)["istarg"]
hist(stat1null,
     xlim=c(min(c(obsstat1, min(stat1null)))-.1,
            max(stat1null)+.1),
            main=paste("Permutation distribution of effect of",
                       "CRISPR interference on STAT1 expression"))
abline(v = obsstat1)

# p-value is then:
(1 + sum(stat1null < obsstat1))/length(stat1null)

# percentiles
quantile(stat1null, probs = seq(0.01, 0.99, 0.01))

# for STAT5A
obsstat5a = coef(m2)["istarg"]
set.seed(1234)
stat5anull = vapply(1:1000, function(x) {
  crilm_shuff("STAT5A", criEC)$coef["istarg"]
}, FUN.VALUE = numeric(1))
(1 + sum(stat5anull < obsstat5a))/length(stat5anull)

# Transcriptome-wide effects of CRISPRi
# use regularized linear models to understand off-target effects
library(limma)

# function to fit linear models to all genes in SE object
crilimmav = function(target, SE) {
  contrids = grep("control", SE$target)
  limSE = SE[, c(contrids, which(SE$target==target)) ]
  fac = paste0(substr(limSE$date_created,1,10), substr(limSE$submitted_by,1,3))
  istarg = rep(0,ncol(limSE))
  istarg[which(limSE$target==target)]=1
  ndf = data.frame(date_conc_submitter=fac, istarg=istarg)
  
  mm = model.matrix(~date_conc_submitter + istarg, data = ndf)
  vv = voom(assay(limSE), design = mm, plot=FALSE)
  f1 = lmFit(vv, mm)
  eBayes(f1)
}

# perform on STAT1 knockdown
li1 = crilimmav("STAT1", criEC[bas,])

# reporter function to make topTable with top 15 most affected genes
report_l = function(clout, SE = criEC, n=15, coef="istarg") {
  options(digits = 3)
  tt = topTable(clout, coef=coef, n=n, sort.by="p")
  tts = rowData(SE[rownames(tt),])$gene_name
  rownames(tt) = tts
  tt
}
report_l(li1)

# repeat for TEAD4 (much larger number of off-target effects)
li2 = crilimmav("TEAD4", criEC[bas,])
report_l(li2)

# repeat for all genes
alltarg = unique(criEC$target)
suppressMessages({
  suppressWarnings({
    allruns = lapply(alltarg, function(x) report_l(crilimmav(x, criEC[bas,])))
  })
})
names(allruns) = alltarg

# order experiments by adjusted p-value
allfd = sapply(allruns, function(x) x[1, "adj.P.Val"])
sort(allfd)[1:10]

# targets with significant findings of interference on at least 15 genes
sigs = sapply(allruns, function(x) ifelse(sum(x[,"adj.P.Val"] < 0.05) == 15, 1, 0))
sum(sigs, na.rm = TRUE)

# looking at a given experiment
allruns$SYNCRIP

# compare distributions of gene impacted by interference on a target
cribox3 = function(target, rider, SE) {
  sind = which(rowRanges(SE)$gene_name == rider)
  contrids = grep("control", SE$target)
  limSE = SE[sind, c(contrids, which(SE$target == target))]
  boxplot(split(as.numeric(assay(limSE)), limSE$target),
          ylab = paste(rider, "expected count"),
          xlab = "interference target")
}
cribox3("SYNCRIP", "GSDMB", criEC[bas,])
