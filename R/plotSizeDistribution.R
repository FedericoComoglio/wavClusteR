#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#2013 - Federico Comoglio & Cem Sievers, D-BSSE, ETH Zurich
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
