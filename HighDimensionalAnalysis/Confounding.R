## Confounding
library(dagdata)
data(admissions)
print(admissions)

# proportions of men and women accepted
index <- which(admissions$Gender == 1)
appliedm <- sum(admissions$Number[index])
acceptedm <- round(sum(admissions$Number[index] * admissions$Percent[index]/100))
acceptedm/appliedm
rejectedm <- appliedm - acceptedm

index <- which(admissions$Gender == 0)
appliedf <- sum(admissions$Number[index])
acceptedf <- round(sum(admissions$Number[index] * admissions$Percent[index]/100))
acceptedf/appliedf
rejectedf <- appliedf - acceptedf

# test for significance using chi-square
tab <- as.table(matrix(c(acceptedm, rejectedm, acceptedf, rejectedf), nrow=2))
chisq.test(tab)

# quantify difficulty of a matrix by proportion of accepted students
H <- sapply(c("A", "B", "C", "D", "E", "F"), function(majr) {
  idx <- which(admissions$Major == majr)
  sum(admissions$Percent[idx] * admissions$Number[idx])/(100*sum(admissions$Number[idx]))
})
H

# correlation between difficulty and genders
index <- which(admissions$Gender == 1)
cor(admissions$Number[index], H)

index <- which(admissions$Gender == 0)
cor(admissions$Number[index], H)
# females more likely to apply to harder majors

## Genomics
library(Biobase)
library(GSE5859)
data(GSE5859)

geneExpression = exprs(e)
sampleInfo = pData(e)

# extract years with multiple ethnicities represented
head(sampleInfo)
year <- format(sampleInfo$date, "%y")
u <- unique(year)
eth <- sampleInfo$ethnicity
yeareth <- matrix(c(year, eth), ncol = 2)
mult_eth <- sapply(1:length(u), function(i) {
  ui <- u[i]
  idx <- which(yeareth[,1] == ui)
  length(unique(yeareth[idx,2]))
})
sum(mult_eth > 1)
which(mult_eth > 1)

# repeat but for month-year combination
month.year <- format(sampleInfo$date, "%m%y")
u <- unique(month.year)
moyreth <- matrix(c(month.year, eth), ncol = 2)
uniques <- sapply(1:length(u), function(i) {
  ui <- u[i]
  idx <- which(moyreth[,1] == ui)
  length(unique(moyreth[idx,2]))
})
mean(uniques > 1)

sample_comp <- function(idx1, idx2) {
  df <- data.frame(geneExpression[,idx1], geneExpression[,idx2])
  g <- c(rep(0, length(idx1)), rep(1, length(idx2)))
  tt <- rowttests(as.matrix(df), fac = as.factor(g))
  library(qvalue)
  qs <- qvalue(tt$p.value)
  paste(sum(qs$qvalues < 0.05), qs$pi0)
}

# compare CEU in 2002 and 2003
library(genefilter)
idx02 <- which(year == "02" & eth == "CEU")
idx03 <- which(year == "03" & eth == "CEU")
sample_comp(idx02, idx03)

# compare CEU in 2003 and 2004
idx03 <- which(year == "03" & eth == "CEU")
idx04 <- which(year == "04" & eth == "CEU")
sample_comp(idx03, idx04)

# compare ASN to CEU
idxASN <- which(eth == "ASN")
idxCEU <- which(eth == "CEU")
sample_comp(idxASN, idxCEU)

# to account for confounding with date, choose one year with both represented
idxASN05 <- which(eth == "ASN" & year == "05")
idxCEU05 <- which(eth == "CEU" & year == "05")
sample_comp(idxASN05, idxCEU05)

# take three random CEU samples from 2002
set.seed(3)
idxASN05 <- which(eth == "ASN" & year == "05")
idxCEU02 <- which(eth == "CEU" & year == "02")
idxCEUsub <- sample(idxCEU02, 3)
sample_comp(idxASN05, idxCEUsub)