#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#2013 - Federico Comoglio & Cem Sievers, D-BSSE, ETH Zurich
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



#' Classify substitutions based on identified RSF interval and return high
#' confidence transitions
#' 
#' Classify genomic positions exhibiting a substitution based on the relative
#' substitution frequency (RSF) interval. The latter is returned by the
#' \code{getExpInterval} function, but can be user-specified through visual
#' inspection of the posterior class probability returned by the same function.
#' 
#' 
#' @usage getHighConfSub(countTable, support, supportStart = NA, supportEnd =
#' NA, substitution = "TC")
#' @param countTable A GRanges object, corresponding to a count table as
#' returned by the \link{getAllSub} function
#' @param support List, consisting of two numeric slots defining the left and
#' right boundaries (start and end values, respectively) of the RSF interval,
#' as returned by the \link{getExpInterval} function.
#' @param supportStart Numeric, if \code{support} not provided, the RSF value
#' determining the left boundary (start) of the RSF interval. Use this argument
#' to specify a user-defined RSF interval.
#' @param supportEnd Numeric, if \code{support} not provided, the RSF value
#' determining the right boundary (end) of the RSF interval. Use this argument
#' to specify a user-defined RSF interval.
#' @param substitution A character indicating which substitution is induced by
#' the experimental procedure (e.g. 4-SU treatment - a standard in PAR-CLIP
#' experiments - induces T to C transitions and hence substitution = 'TC' in
#' this case.)
#' @return a GRanges object containing high confidence substitutions, with
#' strand-specific coverage, counts and RSF values as metadata.
#' @note In the example below, left and right boundaries were arbitrarily
#' chosen as showcase.
#' @author Federico Comoglio and Cem Sievers
#' @seealso \code{\link{getAllSub}}, \code{\link{getExpInterval}}
#' @keywords core
#' @examples
#' 
#' filename <- system.file( "extdata", "example.bam", package = "wavClusteR" )
#' example <- readSortedBam( filename = filename )
#' countTable <- getAllSub( example, minCov = 10, cores = 1 )
#' highConfSub <- getHighConfSub( countTable, supportStart = 0.2, supportEnd = 0.7, substitution = "TC" )
#' highConfSub
#' 
#' @export getHighConfSub
getHighConfSub <- function( countTable, support, supportStart = NA, supportEnd = NA, substitution = 'TC') { 
# filter substitutions to return high confidence ones, based on the identified support of the posterior density that is dominated by experimental induction
#
# Args:
#   countTable:  a GRanges object, corresponding to a count table where each substitution has a corresponding strand-specific coverage and a count value, as returned by the getAllSub function
#   supportStart: numeric, start of the support
#   supportEnd: numeric, end of the support
#   substitution: character, the substitution/transition induced by the experimental procedure
#
# Returns:
#   a GRanges object containing high confidence substitutions and corresponding coverage, counts and rsf values
#
# Error handling
#   if missing support and user provided values are NA, raise an error
	if( missing( support ) & ( is.na( supportStart ) | is.na( supportEnd ) ) )
		stop( 'Missing support and user defined support is NA. Please either use the values returned by the getSupport function or provide values for the supportStart and supportEnd parameters' )

	if( missing( support ) ) {
		supportStart <- supportStart
		supportEnd <- supportEnd
	}
	else {
		supportStart <- support$supportStart
		supportEnd <- support$supportEnd
	}
	
	#1-extract all substitutions of the same type as input substitution
	countTableSub <- countTable[ elementMetadata( countTable )[, 'substitutions'] == substitution ]
	emd <- elementMetadata( countTableSub )
	subst <- emd[, 'substitutions']
	cov <- emd[, 'coverage']
  	count <- emd[, 'count']

	#2-compute relative substitution frequencies
    	rsf <- count / cov 
	
	#3-determine high confidence substitutions
	highConf <- ( rsf >= supportStart & rsf <= supportEnd )

	#4-prepare output GRanges object
	highConfSub <- 	countTableSub[ highConf, -1 ] #copies ranges, coverage and counts
	elementMetadata( highConfSub )[, 'rsf'] <- rsf[ highConf ]
	
	#5-fill metadata slot with the support used for computation
    	metadata( highConfSub ) <- list( supportStart = supportStart, supportEnd = supportEnd )

	return( highConfSub ) 
}
