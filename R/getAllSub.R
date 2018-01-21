#2013 - Federico Comoglio & Cem Sievers, D-BSSE, ETH Zurich

getComplSubst <- function( processedMD ) {
# Return complementary bases for substitutions identified on the minus strand
#
# Args:
#   a matrix, containing the complete set of substitutions identified via MD tag by the processMD function
#
# Returns:
#   a matrix, same dimensions as the input
#
# Error handling
#   ...
	minusStrand <- which( processedMD[, 'strand'] == '-' )
	mds <- processedMD[ minusStrand, 'substitutions' ]
	mds <- DNAStringSet( mds )
	cpl <-  as.character( Biostrings::complement( mds ) )
	processedMD[ minusStrand, 'substitutions' ] <- cpl
	return( processedMD )
}


getCoverageAtSubst <- function( allSubstSplit, cov ) {
# Compute strand-specific coverage at substitutions
#
# Args:
#   allSubstSplit: GRanges object containing all substitutions identified on a given strand
#   cov: Rle object containing the coverage on the same strand
#
# Returns:
#   a GRanges object, identical as the input object, with elementMetadata 'coverage' filled with computed coverage
#
# Error handling
#   ...

	#1-sort substitutions
	allSubstSplit <- sort( allSubstSplit )

	#2-split them by chromosome
	chr <- seqnames( allSubstSplit )
	allSubstSplitByChr <- split( allSubstSplit, chr )
	n <- length( allSubstSplitByChr ) #formerly length( seqnames( allSubstSplitByChr ) )

	message( 'Computing local coverage at substitutions...' )

	#3-compute local coverage (done for each set of substitutions localizing on a given chromosome)
	chrCov <- names( cov ) #names of chr from coverage
	chrSubst <- names( allSubstSplitByChr ) #names of chr from substitutions
	for( i in seq_len(n) ) {
		whichChr <- which( chrCov == chrSubst[ i ] )
		start <- start( allSubstSplitByChr[[ i ]] )
		localCov <- as.numeric( cov[[ whichChr ]][ start ] )
		elementMetadata( allSubstSplitByChr[[ i ]] )[, 'coverage'] <- localCov
	}

	return( allSubstSplitByChr )
}


getCountTable <- function( substWithCov, minCov ) {
# summarize substitutions to produce a count table
#
# Args:
#   substWithCov: GRanges object containing all identified substitutions and their local coverage as produced by the getFilteredSub function
#   minCov: numeric, the minimum coverage required at a genomic position exhibiting a substitution in order to return it in the output.
#
# Returns:
#   a GRanges object, where substitutions have been summarized
#
# Error handling
#   ...

	#1-retain substitutions with coverage >= minCov
	cov <- elementMetadata( substWithCov )[, 'coverage']
	filtered <- substWithCov[ cov >= minCov ]

	#2-remove substitutions involving N's
	remove <- elementMetadata( filtered )[, 'substitutions'] %in% c( 'AN', 'CN', 'GN', 'TN' )
	filtered <- filtered[ !remove ]

	#3-convert to dataframe
	filteredDf <- as.data.frame( filtered )

	#4-remove duplicates. Note: requires to use duplicated as unique over a GRanges does not accounts for uniqueness of elementMetadata
	remove <- duplicated( filteredDf )
	filteredDf <- filteredDf[ !remove, ]

	#5-initialize count table (GRanges)
	countTable <- GRanges( seqnames      = filteredDf[, 'seqnames'],
			       ranges        = IRanges( filteredDf[, 'start'], filteredDf[, 'end'] ),
 	       		       strand        = filteredDf[, 'strand'],
 		               substitutions = filteredDf[, 'substitutions'],
			       coverage      = filteredDf[, 'coverage'] )
	#6-find overlaps between filtered substitutions and count table template (unique ones)
	olaps <- findOverlaps( filtered, countTable )
	subst <- elementMetadata( filtered )[ queryHits( olaps ), 'substitutions' ]
	sh <- subjectHits( olaps )
	splits <- split( subst, sh )
	uniqueSub <- elementMetadata( countTable )[ unique( sh ), 'substitutions' ]
	n <- length( uniqueSub )
	counts <- sapply( seq_len( n ), function(x) length( which( splits[[ x ]] == uniqueSub[ x ] ) ) )
	elementMetadata( countTable )[, 'count'] <- counts

	return( countTable )
}




#' Identify all substitutions observed across genomic positions exhibiting a
#' specified minimum coverage
#'
#' All substitutions observed across genomic positions exhibiting user-defined
#' minimum coverage are extracted and a count table is returned. This function
#' supports parallel computing.
#'
#' The choice of the minimum coverage influences the variance of the relative
#' substitution frequency estimates, which in turn affect the mixture model
#' fit. A conservative value depending on the library size is recommended for a
#' first analysis. Values smaller than 10 have not been tested and are
#' therefore not recommended.
#'
#' @usage getAllSub(sortedBam, minCov = 20, cores = 1)
#' @param sortedBam GRanges object containing aligned reads as returned by
#' \link{readSortedBam}
#' @param minCov An integer defining the minimum coverage required at a genomic
#' position exhibiting a substitution. Genomic positions of coverage less than
#' \code{minCov} are discarded. Default is 20 (see Details).
#' @param cores An integer defining the number of cores to be used for parallel
#' processing, if available. Default is 1.
#' @return A GRanges object containing a count table, where each range
#' correspond to a substitution. The metadata correspond to the following
#' information: \item{substitutions}{observed substitution, e.g. AT, i.e. A in
#' the reference sequence and T in the mapped read.}
#' \item{coverage}{strand-specific coverage.} \item{count}{number of
#' strand-specific substitutions.}
#' @author Federico Comoglio and Cem Sievers, with contributions from Martin Morgan
#' @seealso \code{\link{readSortedBam}}
#' @keywords core
#' @examples
#'
#' filename <- system.file( "extdata", "example.bam", package = "wavClusteR" )
#' example <- readSortedBam(filename = filename)
#' countTable <- getAllSub( example, minCov = 10, cores = 1 )
#' countTable
#'
#' @export getAllSub
getAllSub <- function( sortedBam, minCov = 20, cores = 1 ) {
# Error handling
#   if qseq (read sequences) and MD (MD tag) are present as metadata, raise an error
	emd <- elementMetadata( sortedBam )	#dataframe object
	if( !all( c('qseq', 'tag.MD') %in% colnames( emd ) ) )
		stop( 'read sequences or MD tags missing. Please read in the BAM file using the readSortedBam function' )

#   if minimum coverage < 10, inform the user that this might lead to unstable MLE of the RSF
	if( minCov < 10 )
		message( 'Note: chosing a minimum coverage < 10 might lead to unstable estimate of the relative substitution frequencies values.' )

	#reshape function for later use
#	rsh <- function( x )
#		  	matrix( x, ncol = 4 )

	#1-split sortedBam by strand to compute strand-specific coverage
	sortedBamByStrand <- split( sortedBam, strand( sortedBam ) )[1 : 2]
	covPlus	 <- coverage( sortedBamByStrand[ 1 ] )
	covMinus <- coverage( sortedBamByStrand[ 2 ] )

	#2-retain only those positions exhibiting a substitution and prepare input to extract all info from MD tag
	mds <- emd$tag.MD
	flag <- grepl( pattern = '[A-Za-z]', mds ) #grep any line with at least one letter (numbers only indicate perfect matches). Returns logical

	hasSubst <- sortedBam[ flag ] #select only those positions with a substitution
	substMtx <- split( hasSubst, hasSubst$tag.MD )	#suggested by MM, generates a GRangesList


#	substMtx <- as.matrix( as.data.frame( hasSubst )[, c( 'seqnames', 'start', 'strand', 'qseq', 'MD' )] ) #transforming to matrix makes split >> faster
#	substMtx <- split( substMtx[, c( 'seqnames', 'start', 'strand', 'qseq' )], substMtx[, 'MD'] ) #use of $ not allowed here (atomic vectors)
#	substMtx <- lapply( substMtx, rsh )

	#3-extract information from MD field
	processedMD <- processMD( substMtx, cores = cores ) #parallelized
	processedMD <- getComplSubst( processedMD )
	n <- nrow( processedMD )

	#4-create a GRanges object containing all identified substitutions
	pos <- as.numeric( processedMD[, 'posMismatchesGenome'] )
	allSubst <- GRanges( seqnames      = processedMD[, 'seqnames'],
			     ranges        = IRanges( pos , pos ),
	     		     strand        = processedMD[, 'strand'],
			     substitutions = processedMD[, 'substitutions'],
                             coverage      = rep( Inf, n ) )

	#5-split the object to extract substitutions on +/- strand, to be processed separately
	allSubstSplit <- split( allSubst, strand( allSubst ) )[1 : 2]
	message('   considering the + strand' )
	substPlus <- getCoverageAtSubst( allSubstSplit[[ 1 ]], covPlus )
	message('   considering the - strand' )
	substMinus <- getCoverageAtSubst( allSubstSplit[[ 2 ]], covMinus )

	#6-combine results in a single GRanges object
	substWithCov <- c( substPlus, substMinus )
	substWithCov <- unlist( substWithCov, use.names = FALSE )

	#7-summarize results to produce a count table
	countTable <- getCountTable( substWithCov, minCov )

	return( countTable )
}
