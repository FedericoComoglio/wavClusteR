\name{exportSequences}
\alias{exportSequences}

\title{Export cluster sequences for motif search analysis}

\description{Export cluster sequences for motif search analysis (FASTA format), e.g. using MEME-ChIP}

\usage{exportSequences(clusters, filename = 'wavClusters.fasta')}

\arguments{
    \item{clusters}{GRanges object containing individual clusters as identified by the \link{filterClusters} function}
    \item{filename}{A character defining the BED file name. Default to "wavClusters.fasta"}
}

\value{A FASTA file containing the cluster sequences}

\author{Federico Comoglio}

\seealso{\code{\link{filterClusters}}}

\keyword{postprocessing}
