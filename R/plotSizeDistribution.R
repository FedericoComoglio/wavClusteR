#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#2013 - Federico Comoglio & Cem Sievers, D-BSSE, ETH Zurich
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



#' Plot the distribution of cluster sizes
#' 
#' Produce an histogram of cluster sizes
#' 
#' 
#' @usage plotSizeDistribution( clusters, showCov = FALSE, ... )
#' @param clusters GRanges object containing individual clusters as identified
#' by the \link{getClusters} function
#' @param showCov logical, if TRUE a scatter plot of average cluster coverage
#' vs. cluster size is shown along with a loess fit. Default is FALSE.
#' @param ... Additional parameters to be passed to the \code{hist} function
#' @return Called for its effect, returns a histogram.
#' @author Federico Comoglio
#' @seealso \code{\link{getClusters}}
#' @keywords postprocessing graphics
#' @examples
#' 
#' require(BSgenome.Hsapiens.UCSC.hg19)
#' 
#' data( model, package = "wavClusteR" ) 
#' 
#' filename <- system.file( "extdata", "example.bam", package = "wavClusteR" )
#' example <- readSortedBam( filename = filename )
#' countTable <- getAllSub( example, minCov = 10, cores = 1 )
#' highConfSub <- getHighConfSub( countTable, supportStart = 0.2, supportEnd = 0.7, substitution = "TC" )
#' coverage <- coverage( example )
#' clusters <- getClusters( highConfSub = highConfSub, 
#'                          coverage = coverage, 
#'                          sortedBam = example, 
#' 	                 method = 'mrn', 
#' 	                 cores = 1, 
#' 	                 threshold = 2 ) 
#' 
#' fclusters <- filterClusters( clusters = clusters, 
#' 		             highConfSub = highConfSub, 
#'         		     coverage = coverage,
#' 			     model = model, 
#' 			     genome = Hsapiens, 
#' 		             refBase = 'T', 
#' 		             minWidth = 12 )
#' plotSizeDistribution(fclusters, breaks = 30, col = 'skyblue2')
#' 
#' @export plotSizeDistribution
plotSizeDistribution <- function( clusters, showCov = FALSE,... ) {
# Plot the size distribution of clusters
#
# Args:
#   clusters: GRanges object containing valid clusters as returned by filterClusters
#   showCov: logical, if true a scatter plot of cluster sizes vs. average cluster coverage is returned
# Returns:
#   called for its effect, returns either an histogram or a scatter plot
#
# Error handling
#   ... 
	sizeV <- width( clusters )
	if( !showCov) {
		hist( sizeV, ..., main = 'Size distribution', 
				  xlab = 'Length (bases)', 
				  ylab = 'Number of wavClusters' )
	} else {
		meanCov <- log10( elementMetadata( clusters )[, 'MeanCov'] )
		df <- data.frame( sizeV, meanCov )
		p <- with( df, 
			   ggplot( df, aes( sizeV, meanCov)) + 
		           geom_point() +
			   geom_smooth() +
			   theme_bw() +
		           xlab( 'Length (bases)' ) +
  			   ylab( 'Mean coverage (log10)' )
			 )
		print( p )
	}
}
