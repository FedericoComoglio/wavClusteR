#2013 - Federico Comoglio & Cem Sievers, D-BSSE, ETH Zurich

computePosterior <- function( tableCovCount ) {
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




#' Fit a non-parametric mixture model from all identified substitutions
#' 
#' Estimates the two-component mixture model consisting of the mixing
#' coefficients and the density functions.
#' 
#' 
#' @usage fitMixtureModel(countTable, substitution = "TC")
#' @param countTable A GRanges object, corresponding to a count table as
#' returned by the \link{getAllSub} function
#' @param substitution A character indicating which substitution is induced by
#' the experimental procedure (e.g. 4-SU treatment - a standard in PAR-CLIP
#' experiments - induces T to C transitions and hence substitution = 'TC' in
#' this case.)
#' @return A list containing: \item{l1}{The first mixing coefficient}
#' \item{l2}{The second mixing coefficient} \item{p}{The mixture model}
#' \item{p1}{The first component of the mixture} \item{p2}{The second component
#' of the mixture}
#' @author Federico Comoglio and Cem Sievers
#' @seealso \code{\link{getAllSub}} \code{\link{getExpInterval}}
#' @keywords core model
#' @examples
#' 
#' \dontrun{
#' filename <- system.file( "extdata", "example.bam", package = "wavClusteR" )
#' example <- readSortedBam(filename = filename)
#' countTable <- getAllSub( example, minCov = 10, cores = 1 )
#' 
#' fitMixtureModel( countTable, substitution = "TC" )
#' }
#' 
#' #load and inspect the model
#' data( model )
#' str( model )
#' 
#' #plot densities and estimate the relative substitution frequency support dominated by PAR-CLIP induction
#' getExpInterval( model, bayes = TRUE, plot = TRUE )
#' 
#' @export fitMixtureModel
fitMixtureModel <- function( countTable, substitution = 'TC' ) {
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
