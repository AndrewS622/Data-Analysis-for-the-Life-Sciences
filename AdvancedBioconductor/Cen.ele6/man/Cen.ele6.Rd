\name{Cen.ele6}
\docType{data}
\alias{Cen.ele6}
\title{Annotation package that collates several annotation resources.}
\description{This data object was automatically created by me <me@abc.com>. It
  represents a collection of annotation packages that can be used as a
  single object named Cen.ele6.  This object can be used with the
  standard four accessor method for all AnnotationDbi objects.  Namely:
  \code{cols}, \code{keytype}, \code{keys} and \code{select}.  Users are
  encouraged to read the vignette from the \code{OrganismDbi} package for
  more details.}
\usage{library(Cen.ele6)}
\keyword{datasets}
\examples{
  Cen.ele6
  cls <- columns(Cen.ele6)
  cls
  cls <- cls[c(1,19,45)]
  kts <- keytypes(Cen.ele6)
  kt <- kts[2]
  kts
  ks <- head(keys(Cen.ele6, keytype=kts[2]))
  ks
  res <- select(Cen.ele6, keys=ks, columns=cls, keytype=kt)
  head(res)
}
