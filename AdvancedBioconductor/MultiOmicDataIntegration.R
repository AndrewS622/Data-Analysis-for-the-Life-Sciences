## Integrative Analysis Examples
# Comparing transcription regulators in yeast
library(harbChIP)
data("harbChIP")

# looking at Mcm1 binding
sv <- qqnorm(exprs(harbChIP)[,"MCM1"], main="Mcm1 binding scores")
topb <- names(sort(exprs(harbChIP)[,"MCM1"], decreasing = TRUE)[1:5])
points(sv$x[topb], sv$y[topb], col="red", pch=19)

# gene with strongest Mcm1 binding vs. time and comparing to Mbp1
library(yeastCC)
data("spYCCES")
alp <- spYCCES[, spYCCES$syncmeth=="alpha"]
nm <- names(which.max(exprs(harbChIP)[,"MCM1"]))
nm2 <- names(which.max(exprs(harbChIP)[,"MBP1"]))
plot(exprs(alp)[nm,]~alp$time, ylab=paste0(nm," expression"), type="l", ylim=c(-1,1))
lines(exprs(alp)[nm2,]~alp$time, ylab=paste0(nm2, " expression"),
      col="purple")

# DNA variants under ESRRA binding peaks
# coincidences of GWAS hits and ESRRA binding sites
library(ERBS)
data("GM12878")
library(gwascat)
data("gwrngs19")
fo <- findOverlaps(GM12878, reduce(gwrngs19))
length(fo)

# distinct peaks including GWAS hits
length(unique(queryHits(fo)))

# is this significant?
# reposition GRanges to randomly selected start sites, count number of hits covered
library(ph525x)
library(gwascat)
rg <- reduce(gwrngs19)

# perform repositioning 100 times
set.seed(1234)
rsc <- sapply(1:100, function(x) {
  length(findOverlaps(reposition(GM12878), rg))
})

# estimate p-value
mean(rsc > length(fo))

# since none exceeded observed, p-value is less than:
1/length(rsc)

# siRNA knockdown in pancreatic cancer line
library(GEOquery)
pc1 <- getGEO("GSE35463")[[1]]
pc1$source_name_ch1
colnames(exprs(pc1))[14]
pc1$data_processing
pData(pc1)

# Affy gene 1.0 ST platform
BiocManager::install("hugene10sttranscriptcluster.db")
library(hugene10sttranscriptcluster.db)
select(hugene10sttranscriptcluster.db, 
       keytype = "SYMBOL", 
       keys = "NUPR1", 
       columns = "PROBEID")

# differences between siCtrl and siNupr1 at 4 time points
exp <- exprs(pc1)["8000574",11:14]
ctrl <- exprs(pc1)["8000574",7:10]

# using shaprio test of non-normality
shapiro.test(exp-ctrl)

# use paired t-test to assess shift in mean
t.test(ctrl, exp, paired=TRUE)

## Exploring TCGA data
library(curatedTCGAData)
library(TCGAutils)

# disease code = READ for rectal adenocarcinoma (from diseaseCodes)
readData = curatedTCGAData("READ", 
                           c("RNASeq2GeneNorm", "Mutation", "Methylation_methyl450"), 
                           dry.run = FALSE)
readData

# show map of sample/patient ID (primary), assay (experiment type), and column name (identifier)
sampleMap(readData)

table(sampleMap(readData)$primary)           # number of datasets per sample
table(table(sampleMap(readData)$primary))    # distribution of datasets/sample

# column data contains patient info
clin = colData(readData)
dim(clin)
head(colnames(clin), 10)    # inspect some variable names

# distribution of cancer phenotypes, vitality (1 = deceased)
table(clin$pathology_T_stage)
clin$t_stage = factor(substr(clin$pathology_T_stage,1,2))   # keep only first 2 characters
table(clin$t_stage)
table(clin$vital_status)
table(clin$t_stage, clin$vital_status)

# mutation data
mut = readData[[1]]
mut

# sample IDs from mutation and from overall dataset
mut_samp_ids = colnames(mut)
head(mut_samp_ids)
head(rownames(clin))

# shorten sample IDs to match clinical data
mut_samp_ids = substr(mut_samp_ids,1,12)    # keep only first 12 characters
all(mut_samp_ids %in% rownames(clin))

# check mutation data
rownames(assay(mut))                     # loci for mutations
assay(mut)[1:4, 1:4]
table(assay(mut)[1,], useNA = "ifany")   # note almost all NAs

# access attributes of S4 objects with @ instead of $
mut_assay = mut@assays
class(mut_assay)
length(mut_assay)

# access first assay data
mut_assay[[1]]    # or mut_assay[[mut_samp_ids[1]]]
mut_assay[[1]]$Hugo_Symbol
table(mut_assay[[1]]$Mutation_Status)
table(mut_assay[[1]]$Variant_Classification)

# function to extract patient ID, gene symbol, and mutation type for each patient
mut_df = mapply(function(id, a) {
  d = as.data.frame(mcols(a)[c("Hugo_Symbol", "Variant_Classification")])
  names(d) = c("symbol", "variant_class")
  d$patientID = id
  d
}, id = mut_samp_ids, a = mut_assay, SIMPLIFY = FALSE, USE.NAMES = FALSE)
mut_df = do.call(rbind, mut_df)

# look at mutations
head(mut_df)

mut_tab = table(mut_df$symbol, mut_df$variant_class)
head(mut_tab)

# sum number of mutations per gene of all mutation classes and order
mut_num = apply(mut_tab[, c("Missense_Mutation", "Nonsense_Mutation", 
                            "Frame_Shift_Del", "Frame_Shift_Ins")], 1, sum)

mut_order = order(mut_num, decreasing = TRUE)

# look at top 10 most mutated genes
mut_tab[mut_order[1:10], c("Missense_Mutation", "Nonsense_Mutation", 
                           "Frame_Shift_Del", "Frame_Shift_Ins")]

# linking mutations and tumor stage
# find number of mutations in each patient
nmut = sapply(split(mut_df$patientID, mut_df$patientID),length)
head(nmut)

# not all patients had mutation studies
length(nmut)
nrow(clin)

# which patients did have mutation studies
clinwmut = clin[names(nmut),]

# plot number of mutations vs. stage of cancer
with(clinwmut, boxplot(split(nmut, t_stage), log="y"))

# look at TP53 mutation
tp53_mut_pts = mut_df[mut_df$symbol == "TP53","patientID"]
clinwmut$tp53_mut = clinwmut$patientID %in% tp53_mut_pts

# look at TP53 mutation vs. cancer stage
table(clinwmut$tp53_mut, clinwmut$t_stage)

# linking expression and tumor stage
rnaseq = readData[[2]]
rnaseq
assay(rnaseq)[1:3, 1:3]

# convert to log scale
assay(rnaseq) = log2(assay(rnaseq) + 1)
assay(rnaseq)[1:3, 1:3]

# subset names to clinical patient name
colnames(rnaseq) = substr(colnames(rnaseq),1,12)    # keep only first 12 characters

colData(rnaseq) = clin[colnames(rnaseq),]

# fit model of RNA expression vs. t-stage
library(limma)
mm = model.matrix(~t_stage, data=colData(rnaseq))
f1 = lmFit(assay(rnaseq), mm)
ef1 = eBayes(f1)
topTable(ef1)

# look at example down- and up-regulated genes
par(mfrow=c(1,2))
boxplot(split(assay(rnaseq)["CNDP2",], rnaseq$t_stage), main="CNDP2")    # higher expression in lower t_stage
boxplot(split(assay(rnaseq)["TAC1",], rnaseq$t_stage), main="TAC1")    # higher expression in higher t_stage

# linking methylation and expression
methyl = readData[[3]]
methyl
assay(methyl)

# take only primary tumor samples ("01A")
isprimary = sapply(strsplit(colnames(methyl), split = "-"), `[[`, 4) == "01A"
methyl = methyl[, isprimary]

# substring and put in data into colData from clin
colnames(methyl) = substr(colnames(methyl),1,12)

colData(methyl) = clin[colnames(methyl),]

# samples in both methylation and RNA-seq datasets
length(intersect(colnames(methyl), colnames(rnaseq)))

# find overlapping patients
methyl_subset = methyl[,which(colnames(methyl) %in% colnames(rnaseq))]
rnaseq_subset = rnaseq[,which(colnames(rnaseq) %in% colnames(methyl))]

# genes in both datasets
methyl_genes = rowData(methyl_subset)$Gene_Symbol
head(methyl_genes)

# function to find correlations between methylation and gene expression for a given gene
me_rna_cor = function(sym, mpick = 3){
  # subset methylation data to first mpick methylation sites for given gene symbol
  methyl_ind = which(methyl_genes == sym)
  if (length(methyl_ind) > mpick){    
    methyl_ind = methyl_ind[1:mpick]
  }
  methyl_dat = assay(methyl_subset)[methyl_ind,]    # subset to selected methylation sites
  
  # subset expression data to selected gene symbol
  expr_ind = which(rownames(rnaseq_subset) == sym)    
  expr_dat = assay(rnaseq_subset)[expr_ind,]
  
  # combine methylation and expression data as data frame
  combined_dat = as(t(methyl_dat), "DataFrame")
  combined_dat$expr = expr_dat
  
  # plot pairs and calculate correlation coefficients between methylation and expression
  pairs(combined_dat)
  sapply(1:mpick, function(i){
    cor(as.numeric(combined_dat[,i]), combined_dat[,mpick+1])
  })
}

me_rna_cor("TAC1", mpick=2)

# Assessment
library(ph525x)
data(package = "ph525x")
# readES is ExpressionSet from rectal adenocarcinoma
# readMuts is df with mutation info from rectal adenocarcinoma

options(digits = 3)
data(readES)
pData(readES)
rownames(pData(readES))

# gene-wise FDRs for constant expression across tumor stages
library(limma)
mm <- model.matrix(~t_stage, data=pData(readES))
f1 <- lmFit(readES, mm)
ef1 <- eBayes(f1)
topTable(ef1, 2:4, n=20)

# significant genes
sum(topTable(ef1, 2:4, n=Inf)$adj.P.Val < 0.105)

# make t_stage a numerical score
readES$numts <- as.numeric(factor(readES$t_stage))
mm2 <- model.matrix(~numts, data=pData(readES))
f2 <- lmFit(readES, mm2)
ef2 <- eBayes(f2)
toptable(ef2, 2, n=50)

# testing slope of regression line significance
sum(topTable(ef2, 2, n=Inf)$adj.P.Val < 0.105)

# interpreting a significant gene
boxplot(exprs(readES)["COMP",]~readES$t_stage)

# looking at mutation dataset
data(readMuts)
length(unique(readMuts$Tumor_Sample_Barcode))

# different mutation types
sort(table(readMuts$Variant_Classification), decreasing = TRUE)

# filter to KRAS gene
idx <- which(readMuts$Hugo_Symbol == "KRAS")
mutsKRAS <- readMuts[idx,]

# what amino acid mutation is most common
sort(table(mutsKRAS$AAChange), decreasing = TRUE)

# convert patient IDs from readMuts to same as readES
mut_id <- readMuts$Tumor_Sample_Barcode
mut_id <- gsub("-", ".", mut_id)
mut_id <- unique(tolower(mut_id))
mut_id <- substr(mut_id, 1, 12)
idx <- unique(mut_id) %in% rownames(pData(readES))
sum(idx)
