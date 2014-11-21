#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#2013 - Federico Comoglio & Cem Sievers, D-BSSE, ETH Zurich
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
