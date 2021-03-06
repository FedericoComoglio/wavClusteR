% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getClusters.R
\name{getClusters}
\alias{getClusters}
\title{Identify clusters containing high-confidence substitutions and resolve
boundaries at high resolution}
\usage{
getClusters(highConfSub, coverage, sortedBam, cores =
1, threshold)
}
\arguments{
\item{highConfSub}{GRanges object containing high-confidence substitution
sites as returned by the \link{getHighConfSub} function}

\item{coverage}{An Rle object containing the coverage at each genomic
position as returned by a call to \link{coverage}}

\item{sortedBam}{a GRanges object containing all aligned reads, including
read sequence (qseq) and MD tag (MD), as returned by the
\link{readSortedBam} function}

\item{cores}{integer, the number of cores to be used for parallel
evaluation. Default is 1.}

\item{threshold}{numeric, the difference in
coverage to be considered noise. If not specified, a Gaussian mixture model
is used to learn a threshold from the data. Empirically, 10\% of the minimum
coverage required at substitutions (see argument \code{minCov} in the
\link{getHighConfSub} function) might suffice to provide highly resolved
clusters. However, if \code{minCov} is much lower than the median
strand-specific coverage at substitutions \eqn{m}, which can be computed
using \code{summary(elementMetadata(highConfSub)[, 'coverage'])['Median']}),
10\% of \eqn{m} might represent an optimal choice.}
}
\value{
GRanges object containing the identified cluster boundaries.
}
\description{
Identifies clusters using the mini-rank norm (MRN) algorithm,
which employs thresholding of
background coverage differences and finds the optimal cluster boundaries by
exhaustively evaluating all putative clusters using a rank-based approach.
This method has higher sensitivity and an approximately 10-fold faster
running time than the CWT-based cluster identification algorithm.
}
\note{
Clusters returned by this function need to be further merged by the
function \code{filterClusters}, which also computes all relevant cluster
statistics.
}
\examples{

filename <- system.file( "extdata", "example.bam", package = "wavClusteR" )
example <- readSortedBam( filename = filename )
countTable <- getAllSub( example, minCov = 10, cores = 1 )
highConfSub <- getHighConfSub( countTable, supportStart = 0.2, supportEnd = 0.7, substitution = "TC" )
coverage <- coverage( example )
clusters <- getClusters( highConfSub = highConfSub, 
                         coverage = coverage, 
                         sortedBam = example, 
	                 cores = 1, 
	                 threshold = 2 ) 

}
\references{
Sievers C, Schlumpf T, Sawarkar R, Comoglio F and Paro R. (2012) Mixture
models and wavelet transforms reveal high confidence RNA-protein interaction
sites in MOV10 PAR-CLIP data, Nucleic Acids Res. 40(20):e160. doi:
10.1093/nar/gks697

Comoglio F, Sievers C and Paro R (2015) Sensitive and highly resolved identification
of RNA-protein interaction sites in PAR-CLIP data, BMC Bioinformatics 16, 32.
}
\seealso{
\code{\link{getHighConfSub}}, \code{\link{filterClusters}}
}
\author{
Federico Comoglio and Cem Sievers
}
\keyword{core}
