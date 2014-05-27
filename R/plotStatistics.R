#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#2013 - Federico Comoglio & Cem Sievers, D-BSSE, ETH Zurich
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

plotStatistics <- function( clusters, corMethod = 'spearman', lower = panel.smooth, ... ) {
# Pairs plot the metadata of clusters
#
# Args:
#   clusters: GRanges object containing individual cluster statistics as returned by filterClusters
#   corMethod: character, one out of 'pearson', 'spearman' or 'kendall'. Default is ranked-based correlation ('spearman')
#   lower: function to generate lower panel in the pairs plot. Default is panel.smooth
#   ...: additional parameters to be passed to the pairs function
#
# Returns:
#   Called for its effect, returns a pairs plot of cluster statistics
#
# Error handling:
#   verify that the clusters object contains statistics, otherwise raise an error
	emd <- as.data.frame( elementMetadata( clusters ) )
	stopifnot( 'Ntransitions' %in% colnames( emd ) )

	emd <- emd[, -5 ] #remove Sequence

	panelCor <- function( x, y, digits = 2, prefix = '', ... ) {
		usr <- par( 'usr' )
		on.exit( par( usr ) )
    		par( usr = c(0, 1, 0, 1) )
    		r <- cor( x, y, method = 'spearman' )
    		txt <- format( c( r, 0.123456789 ), digits = digits )[ 1 ]
    		txt <- paste0( prefix, txt )
		cexCor <- 0.8 / strwidth( txt )
	    	text( 0.5, 0.5, txt, cex = cexCor * abs( r ) )
	}
	
	pairs( emd, 
	       lower.panel = lower, 
               upper.panel = panelCor, ... )

}

