#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#2013 - Federico Comoglio & Cem Sievers, D-BSSE, ETH Zurich
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



#' Pairs plot visualization of clusters statistics
#' 
#' Graphical representation of cluster statistics, featuring pairwise
#' correlations in the upper panel.
#' 
#' 
#' @usage plotStatistics(clusters, corMethod = 'spearman', lower =
#' panel.smooth, ...)
#' @param clusters GRanges object containing individual clusters as identified
#' by the \link{getClusters} function
#' @param corMethod A character defining the correlation coefficient to be
#' computed. See the help page of the \code{cor} function for possible options.
#' Default is "spearman". Hence, rank-based Spearman's correlation coefficients
#' are computed
#' @param lower A function compatible with the \code{lower} panel argument of
#' the \code{pairs} function
#' @param ... Additional parameters to be passed to the \code{pairs} function
#' @return called for its effect
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
#' plotStatistics( clusters = fclusters )
#' 
#' @export plotStatistics
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

