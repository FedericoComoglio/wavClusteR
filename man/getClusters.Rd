\name{getClusters}
\alias{getClusters}

\title{Identify clusters containing high-confidence substitutions and resolve boundaries at high resolution}

\description{Identifies clusters using either the mini-rank norm (MRN) algorithm (default and recommended to achieve highest sensitivity) or via a continuous wavelet transform (CWT) based approach. The former employs thresholding of background coverage differences and finds the optimal cluster boundaries by exhaustively evaluating all putative clusters using a rank-based approach. This method has higher sensitivity and an approximately 10-fold faster running time than the CWT-based cluster identification algorithm. The latter, maintained for compatibility with \code{wavClusteR}, computes the CWT on a 1 kb window of the coverage function centered at a high-confidence substitution site, and identifies cluster boundaries by extending away from peak positions.}

\usage{getClusters(highConfSub, coverage, sortedBam, method = 'mrn', 
       cores = 1, threshold, step = 1, snr = 3)}

\arguments{
  \item{highConfSub}{GRanges object containing high-confidence substitution sites as returned by the \link{getHighConfSub} function}
  \item{coverage}{An Rle object containing the coverage at each genomic position as returned by a call to \link{coverage}}
  \item{sortedBam}{a GRanges object containing all aligned reads, including read sequence (qseq) and MD tag (MD), as returned by the \link{readSortedBam} function}
  \item{method}{a character, either set to "mrn" or to "cwt" to compute clusters using the mini-rank norm or the wavelet transform-based algorithm, respectively. Default is "mrn" (recommended). }
  \item{cores}{integer, the number of cores to be used for parallel evaluation. Default is 1.}
  \item{threshold}{numeric, if \code{method = "mrn"}, the difference in coverage to be considered noise. If not specified, a Gaussian mixture model is used to learn a threshold from the data. Empirically, 10\% of the minimum coverage required at substitutions (see argument \code{minCov} in the \link{getHighConfSub} function) might suffice to provide highly resolved clusters. However, if \code{minCov} is much lower than the median strand-specific coverage at substitutions \eqn{m}, which can be computed using \code{summary(elementMetadata(highConfSub)[, 'coverage'])['Median']}), 10\% of \eqn{m} might represent an optimal choice.}
  \item{step}{numeric, if \code{method = "cwt"}, step size of window shift. If two high-confidence substitution sites are located within a distance less than \code{step}, the wavelet transform is computed only once. Default: 1, i.e. each high-confidence substitution site is considered independently.}
  \item{snr}{numeric, if \code{method = "cwt"}, signal-to-noise ratio controlling the peak calling as performed by \code{wavCWTPeaks} implemented in the wmtsa package. Default: 3.}
}

\value{GRanges object containing the identified cluster boundaries.}

\references{William Constantine and Donald Percival (2011), wmtsa: Wavelet Methods for Time Series Analysis,
\url{http://CRAN.R-project.org/package=wmtsa}}

\author{Federico Comoglio and Cem Sievers}

\note{Clusters returned by this function need to be further merged by the function \code{filterClusters}, which also computes all relevant cluster statistics.}

\seealso{\code{\link{getHighConfSub}}, \code{\link{filterClusters}}}

\examples{
filename <- system.file( "extdata", "example.bam", package = "wavClusteR" )
example <- readSortedBam( filename = filename )
countTable <- getAllSub( example, minCov = 10, cores = 1 )
highConfSub <- getHighConfSub( countTable, supportStart = 0.2, supportEnd = 0.7, substitution = "TC" )
coverage <- coverage( example )
clusters <- getClusters( highConfSub = highConfSub, 
                         coverage = coverage, 
                         sortedBam = example, 
	                 method = 'mrn', 
	                 cores = 1, 
	                 threshold = 2 ) 
}

\keyword{core}
