\name{exportCoverage}
\alias{exportCoverage}

\title{Export coverage as BigWig track}

\description{Export coverage as BigWig track, compatible with the UCSC genome browser}

\usage{exportCoverage(coverage, filename = 'wavClusters.BigWig')}

\arguments{
  \item{coverage}{An Rle object containing the coverage at each genomic position as returned by a call to \code{coverage}}
  \item{filename}{A character defining the BED file name. Default to "wavClusters.BigWig"}
}

\value{A BigWig file of the exported Rle object}

\author{Federico Comoglio}

\keyword{postprocessing}
