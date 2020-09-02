\name{Sac.cer3}
\docType{data}
\alias{Sac.cer3}
\title{Annotation package that collates several annotation resources.}
\description{This data object was automatically created by Student <ph525x@harvardx.edu>. It
  represents a collection of annotation packages that can be used as a
  single object named Sac.cer3.  This object can be used with the
  standard four accessor method for all AnnotationDbi objects.  Namely:
  \code{cols}, \code{keytype}, \code{keys} and \code{select}.  Users are
  encouraged to read the vignette from the \code{OrganismDbi} package for
  more details.}
\usage{library(Sac.cer3)}
\keyword{datasets}
\examples{
  Sac.cer3
  cls <- columns(Sac.cer3)
  cls
  cls <- cls[c(1,19,45)]
  kts <- keytypes(Sac.cer3)
  kt <- kts[2]
  kts
  ks <- head(keys(Sac.cer3, keytype=kts[2]))
  ks
  res <- select(Sac.cer3, keys=ks, columns=cls, keytype=kt)
  head(res)
}
