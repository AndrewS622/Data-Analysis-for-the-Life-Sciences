juxta <-
function (chrname="chr22", ...) 
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
