\name{exportClusters}
\alias{exportClusters}

\title{Export clusters as BED track}

\description{Export clusters as BED track, compatible with the UCSC genome browser}

\usage{exportClusters(clusters, filename = 'wavClusters.bed', 
	       trackname = 'wavClusters', description = 'wavClusters')}

\arguments{
  \item{clusters}{GRanges object containing individual clusters as identified by the \link{filterClusters} function}
  \item{filename}{A character defining the BED file name. Default to "wavClusters.bed"}
  \item{trackname}{A character defining the track.name of the BED file. Default to "wavClusters"}
  \item{description}{A character defining the description of the BED file. Default to "wavClusters"}
}

\value{A BED file of the exported GRanges object}

\note{Clusters are color coded according to their strand information (red for the plus strand, blue for the minus strand).}

\author{Federico Comoglio}

\seealso{\code{\link{filterClusters}}}

\keyword{postprocessing}
