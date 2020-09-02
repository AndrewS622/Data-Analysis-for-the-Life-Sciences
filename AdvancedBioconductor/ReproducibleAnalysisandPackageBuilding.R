## Package Design and Creation
# Creating a function
juxta = function (chrname="chr22", ...) 
{
  require(ERBS)
  data(HepG2)
  data(GM12878)
  require(ggbio)
  require(GenomicRanges)
  ap1 <- autoplot(HepG2[seqnames(HepG2) == chrname])
  ap2 <- autoplot(GM12878[seqnames(GM12878) == chrname])
  tracks(HepG2 = ap1, Bcell = ap2, ...)
}

# creating a package
package.skeleton("erbsViz", list="juxta")
# after updating man pages
install.packages("erbsViz", repos=NULL, type="source")
library(erbsViz)

# look at example
jdemo = juxta()
class(jdemo)
getSlots(getClass(class(jdemo)))
jdemo

# OrganismDb for C. elegans
library(OrganismDbi)
gd <- list(join1 = c(GO.db="GOID", org.Ce.eg.db="GO"),
           join2 = c(org.Ce.eg.db="ENTREZID",
           TxDb.Celegans.UCSC.ce6.ensGene="GENEID"))
# ensure all C. elegans packages have been installed
makeOrganismPackage("Cen.ele6", gd, "C. elegans", "1.0.0", "me <me@abc.com>", "me <me@abc.com>", ".")
install.packages("Cen.ele6", repos=NULL, type="source")
library(Cen.ele6)
sum(seqlengths(Cen.ele6))

# OrganismDb for Yeast
gd = list( join1 = c(GO.db="GOID", org.Sc.sgd.db="GO"),
           join2 = c(org.Sc.sgd.db="ENTREZID",
                     TxDb.Scerevisiae.UCSC.sacCer3.sgdGene="GENEID"))
makeOrganismPackage(pkgname="Sac.cer3",  # simplify typing!
                    graphData=gd, organism="Saccharomyces cerevisiae",
                    version="1.0.0", maintainer="Student <ph525x@harvardx.edu>",
                    author="Student <ph525x@harvardx.edu>",
                    destDir=".",
                    license="Artistic-2.0")
install.packages("Sac.cer3", repos=NULL, type="source")
library(Sac.cer3)