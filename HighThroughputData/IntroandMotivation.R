## R Refresher
library(devtools)
install_github("genomicsclass/GSE5859Subset")
library(GSE5859Subset)
data(GSE5859Subset)

sum(na.omit(match(sampleInfo$date, as.Date("2005-06-27"))))

sum(na.omit(match(geneAnnotation$CHR, "chrY")))

subj <- which(sampleInfo$date == as.Date("2005-06-10"))
gene <- which(geneAnnotation$SYMBOL == "ARPC1A")
geneExpression[gene,subj]

median(apply(geneExpression, 2, median))

g <- factor(sampleInfo$group)
pval <- function(e, group = g) {
  t.test(e[group == 1], e[group == 0])$p.value
}

pvals <- apply(geneExpression, 1, pval)
min(pvals)

## Inference in Practice
set.seed(1)
library(downloader)
url = "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/femaleControlsPopulation.csv"
filename = "femaleControlsPopulation.csv"
if (!file.exists(filename)) download(url,destfile=filename)
population <- read.csv(filename)
pvals <- replicate(1000, {
  control <- sample(population[,1],12)
  treatment <- sample(population[,1],12)
  t.test(treatment,control)$p.value
})
head(pvals)
hist(pvals)
mean(pvals < 0.01)
# p-values are uniformly distributed random variables

MCmouse <- function() {
  cases <- rnorm(10,30,2)
  controls <- rnorm(10,30,2)
  t.test(cases,controls)$p.value
}
set.seed(100)
mousepvals <- replicate(20, MCmouse())
sum(mousepvals < 0.05)

set.seed(100)
MCmouseDiets <- function() {
  mousepvals <- replicate(20, MCmouse())
  sum(mousepvals < 0.05)
}
mouseDietpvals <- replicate(1000, MCmouseDiets())
mean(mouseDietpvals)
mean(mouseDietpvals>0)
# 1.007/20 p-values significant at 0.05 level on average
