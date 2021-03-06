% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getHighConfSub.R
\name{getHighConfSub}
\alias{getHighConfSub}
\title{Classify substitutions based on identified RSF interval and return high
confidence transitions}
\usage{
getHighConfSub(countTable, support, supportStart = NA, supportEnd =
NA, substitution = "TC")
}
\arguments{
\item{countTable}{A GRanges object, corresponding to a count table as
returned by the \link{getAllSub} function}

\item{support}{List, consisting of two numeric slots defining the left and
right boundaries (start and end values, respectively) of the RSF interval,
as returned by the \link{getExpInterval} function.}

\item{supportStart}{Numeric, if \code{support} not provided, the RSF value
determining the left boundary (start) of the RSF interval. Use this argument
to specify a user-defined RSF interval.}

\item{supportEnd}{Numeric, if \code{support} not provided, the RSF value
determining the right boundary (end) of the RSF interval. Use this argument
to specify a user-defined RSF interval.}

\item{substitution}{A character indicating which substitution is induced by
the experimental procedure (e.g. 4-SU treatment - a standard in PAR-CLIP
experiments - induces T to C transitions and hence substitution = 'TC' in
this case.)}
}
\value{
a GRanges object containing high confidence substitutions, with
strand-specific coverage, counts and RSF values as metadata.
}
\description{
Classify genomic positions exhibiting a substitution based on the relative
substitution frequency (RSF) interval. The latter is returned by the
\code{getExpInterval} function, but can be user-specified through visual
inspection of the posterior class probability returned by the same function.
}
\note{
In the example below, left and right boundaries were arbitrarily
chosen as showcase.
}
\examples{

filename <- system.file( "extdata", "example.bam", package = "wavClusteR" )
example <- readSortedBam( filename = filename )
countTable <- getAllSub( example, minCov = 10, cores = 1 )
highConfSub <- getHighConfSub( countTable, supportStart = 0.2, supportEnd = 0.7, substitution = "TC" )
highConfSub

}
\seealso{
\code{\link{getAllSub}}, \code{\link{getExpInterval}}
}
\author{
Federico Comoglio and Cem Sievers
}
\keyword{core}
