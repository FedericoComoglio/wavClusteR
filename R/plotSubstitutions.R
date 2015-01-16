#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#2013 - Federico Comoglio & Cem Sievers, D-BSSE, ETH Zurich
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



#' Barplot visualization of the number of genomic positions exhibiting a given
#' substitution and, if model provided, additional diagnostic plots.
#' 
#' Graphical representation of the total number of genomic positions exhibiting
#' one or more substitutions of a given type. This information is used to
#' estimate the mixing coefficients of the non-parametric mixture model. If the
#' mixture model fit is provided, returns additional diagnostic plots such as
#' the total number of reads exhibiting a given substitution and relative
#' substitution frequency-dependent representations of the total number of
#' genomic positions with substitutions of a given type.
#' 
#' 
#' @usage plotSubstitutions(countTable, highlight = "TC", model)
#' @param countTable A GRanges object, corresponding to a count table as
#' returned by the \link{getAllSub} function
#' @param highlight A character indicating which substitution should be
#' highlighted in the barplot. A standard PAR-CLIP experiment employing 4-SU
#' treatment induces T to C transitions, encoded as "TC". Default is "TC".
#' @param model A list containing the model as returned by the function
#' \code{fitMixtureModel}
#' @return called for its effect
#' @author Federico Comoglio and Cem Sievers
#' @seealso \code{\link{getAllSub}}
#' @keywords graphics
#' @examples
#' 
#' filename <- system.file( "extdata", "example.bam", package = "wavClusteR" )
#' example <- readSortedBam( filename = filename )
#' countTable <- getAllSub( example, minCov = 10, cores = 1 )
#' plotSubstitutions(countTable = countTable, highlight = "TC")
#' 
#' @export plotSubstitutions
plotSubstitutions <- function( countTable, highlight = 'TC', model ) {
# produce barplot of observed substitutions. If model is not supplied, returns the simplest diagnostic plot. Produces four diagnostic plots otherwise.
#
# Args:
#   countTable:  a GRanges object, corresponding to a count table where each substitution has a corresponding strand-specific coverage and a count value, as returned by the getFilteredSub function
#   highlight: character, the substitution to be highlighted in the plot
#   model: list, mixture model as returned by fitMixtureModel
#
# Returns:
#   called for its effect, returns a barplot
#
# Error handling
#   ...

	#1-extract substitutions and compute summary table
	countTable <- countTable[ !elementMetadata( countTable )[, 'substitutions'] %in% c( 'AN', 'CN', 'GN', 'TN' ) ]
	emd <- elementMetadata( countTable )
	subst <- emd[, 'substitutions']
	count <- emd[, 'count']
	countPos <- table( subst )
	n <- length( countPos )

	#2-prepare plot and highlight transition of interest
	substNames <- names( countPos )
	posHL <- which( substNames == highlight ) #position of the substitution to highlight
	percentagePos <- round( countPos[ posHL ] / sum( countPos ) * 100, 2 )
	col <- rep( 'gray60', n )
	col[ posHL ] <- 'skyblue2'

	#3-compute extra diagnostics if model supplied
	if( !missing( model ) ) {
		#total number of reads carrying a transition
		countReads <- sapply( split( count, subst ), sum )
		percentageReads <- round( countReads[ posHL ] / sum( countReads ) * 100, 2 )
		#distribution of genomic positions within or outside hc support
		support <- getExpInterval( model, plot = FALSE )
		rsf <- count / emd[, 'coverage']
		#within
		within <- ( rsf >= support$supportStart ) & ( rsf <= support$supportEnd )
print(table(within))
		substIn <- subst[ within ]
		countPosIn <- table( substIn )
		substNamesIn <- names( countPosIn )
		#outside
		substOut <- subst[ !within ]
		countPosOut <- table( substOut )
		substNamesOut <- names( countPosOut )
				
		par( mfrow = c( 2, 2 ) )
		
		barplot( countPos,
                 	 cex.names = 0.8, 
	   	 	 names.arg = substNames,
		  	 main      = paste0( 'Substitutions (', highlight, ' = ', percentagePos, ' %)' ), 
			 ylab      = 'Number of genomic positions',
			 col       = col )

		barplot( countReads,
                 	 cex.names = 0.8, 
	   	 	 names.arg = substNames,
		  	 main      = paste0( 'Substitutions (', highlight, ' = ', percentageReads, '%)', sep = '' ), 
			 ylab      = 'Number of reads with substitution',
			 col       = col )

		barplot( countPosIn,
                 	 cex.names = 0.8, 
	   	 	 names.arg = substNamesIn,
		  	 main      = paste0( 'RSF in [', support$supportStart, ',', support$supportEnd, ']' ), 
			 ylab      = 'Number of genomic positions',
			 col       = col )
		
		barplot( countPosOut,
                 	 cex.names = 0.8, 
	   	 	 names.arg = substNamesOut,
		  	 main      = paste0( 'RSF not in [', support$supportStart, ',', support$supportEnd, ']' ), 
			 ylab      = 'Number of genomic positions',
			 col       = col )	
	} else {
		barplot( countPos,
                 	 cex.names = 0.8, 
	   	 	 names.arg = substNames,
		  	 main      = paste0( 'Substitutions (', highlight, ' = ', percentagePos, ' %)' ), 
			 ylab      = 'Number of genomic positions',
			 col       = col )
	}
}	
