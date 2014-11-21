#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#2013 - Federico Comoglio & Cem Sievers, D-BSSE, ETH Zurich
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

computePosterior <- function( tableCovCount ) {
# compute posterior densities
#
# Args:
#   tableCovCount: a table object, as produced within the estimateP function
#
# Returns:
#   a vector, with values of the estimated posterior density
#
# Error handling
#   ...

	matrixCovCount <- as.matrix( tableCovCount )

    	covVal <- as.numeric( rownames( matrixCovCount ) )
    	substVal <- as.numeric( colnames( matrixCovCount ) )

	mu <- seq( 0.001, 0.999, by = 0.001 )
	logmu <- log( mu )
	log1mu <- log( 1 - mu )
	a <- 1
	b <- 1

	nonZero <- which( matrixCovCount != 0, arr.ind = TRUE)

	posterior <- apply( nonZero, 1, function( pair ) {
		i <- pair[ 1 ]
		j <- pair[ 2 ]
		thisCov <- covVal[ i ]
		thisSubst <- substVal[ j ]
		l <- thisCov - thisSubst

		logFact <- lgamma (a + thisSubst + b + l ) - lgamma( thisSubst + a ) - lgamma( l + b )
		logPost <- logFact + ( a + thisSubst - 1 ) * logmu + ( b + l - 1 ) * log1mu
		term <- exp( logPost )
		factor <- matrixCovCount[ i, j ]
		factor * term	
	} )

	posterior <- rowSums( posterior )
	posterior <- posterior / ( sum(posterior) * 0.001 )
	return( posterior )
}


estimateP <- function( countTableSplit ) { 
# wrapper to estimate posterior densities
#
# Args:
#   countTableSplit:  a GRanges object, corresponding to a count table where each substitution has a corresponding strand-specific coverage and a count value, as returned by the getFilteredSub function and as further split by the fitMixtureModel function
#
# Returns:
#   a vector, with values of the estimated posterior density
#
# Error handling
#   ...

	emd <- elementMetadata( countTableSplit )
	cov <- emd[, 'coverage']
	count <- emd[, 'count']
	tableCovCount <- table( cov, count )
  	posterior <- computePosterior( tableCovCount )
  
	return( posterior )
}


fitMixtureModel <- function( countTable, substitution = 'TC' ) {
# fit non-parametric mixture model and returns full density, components and mixing coefficients
#
# Args:
#   countTable:  a GRanges object, corresponding to a count table where each substitution has a corresponding strand-specific coverage and a count value, as returned by the getFilteredSub function
#   substitution: character, the specific transition induced by the experimental workflow. Default to 'TC' for 4-SU.
#
# Returns:
#   a list, where slots correspond to mixing coefficients, full density and two individual densities
#
# Error handling
#   ...

	#1-extract all substitutions and construct a summary table
	subst <- elementMetadata( countTable )[ ,'substitutions' ]
	substTable <- table( subst )

	#2-identifies the entry in the summary table corresponding to the substitution to be considered
	pos <- which( names( substTable ) == substitution )
	
	#3-get counts for experimentally induced substitutions and estimate non-experimentally induced ones
	expSub <- substTable[ pos ]
	nonexpSub <- max( substTable[ -pos ] )

	#4-estimate mixing coefficients
	l2 <- ( expSub - nonexpSub ) / expSub
	l1 <- 1 - l2

	#5-split countTable by substitution. [[1]] is FALSE, [[2]] is true
	countTableSplit <- split( countTable, subst == substitution)

	#6-estimate model densities	
	message( 'Estimating model densities...' )

	p1 <- estimateP( countTableSplit[[1]] )
	p <- estimateP( countTableSplit[[2]] )
	p2 <- (p - l1 * p1) / l2
	p2[p2 < 0] <- 0 

	return( list(l1 = l1, l2 = l2, p = p, p1 = p1, p2 = p2) )
}
