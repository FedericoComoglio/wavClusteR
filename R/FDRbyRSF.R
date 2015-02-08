#2015 - Federico Comoglio & Cem Sievers, D-BSSE, ETH Zurich

#' Estimate False Discovery Rate within the relative substitution frequency
#' support by integrating PAR-CLIP data and RNA-Seq data
#' 
#' Estimate upper and lower bounds for the False Discovery Rate within the
#' relative substitution frequency (RSF) support by integrating PAR-CLIP data
#' and RNA-Seq data (current version makes use of unstranded RNA-Seq)
#' 
#' For details on the FDR computation, please see Comoglio, Sievers and Paro.
#' 
#' @usage estimateFDR(countTable, RNASeq, substitution = 'TC', minCov = 20,
#' span = 0.1, cores = 1, plot = TRUE, verbose = TRUE, ...)
#' @param countTable A GRanges object, corresponding to a count table as
#' returned by the \link{getAllSub} function
#' @param RNASeq GRanges object containing aligned RNA-Seq reads as returned by
#' \link{readSortedBam}
#' @param substitution A character indicating which substitution is induced by
#' the experimental procedure (e.g. 4-SU treatment - a standard in PAR-CLIP
#' experiments - induces T to C transitions and hence substitution = 'TC' in
#' this case.)
#' @param minCov An integer defining the minimum coverage required at a genomic
#' position exhibiting a substitution. Genomic positions of coverage less than
#' \code{minCov} are discarded. Default is 20 (see Details).
#' @param span A numeric indicating the width of RSF intervals to be considered
#' for FDR computation. Defauls is 0.1 (i.e. 10 intervals are considered
#' spanning the RSF support (0,1]
#' @param cores An integer defining the number of cores to be used for parallel
#' processing, if available. Default is 1.
#' @param plot Logical, if TRUE a dotchart with cluster annotations is produced
#' @param verbose Logical, if TRUE processing steps are printed
#' @param ... Additional parameters to be passed to the \code{plot} function
#' @return A list with three slots, containing upper and lower FDR bounds, and
#' the total number of positive instances each RSF interval. If \code{plot},
#' these three vectors are depicted as a line plot.
#' @note The approach used to compute the upper bound for the FDR is very
#' conservative. See supplementary information in Comoglio et al. for details.
#' @author Federico Comoglio and Cem Sievers
#' @seealso \code{\link{readSortedBam}}, \code{\link{getAllSub}}
#' Comoglio F, Sievers C and Paro R (2015) Sensitive and highly resolved identification
#' of RNA-protein interaction sites in PAR-CLIP data, BMC Bioinformatics 16, 32.
#' @keywords core graphics
#' @export estimateFDR

estimateFDR <- function( countTable, RNASeq, substitution = 'TC', minCov = 20, span = 0.1, cores = 1, plot = TRUE, verbose = TRUE, ... ) {
# Error handling
#   ...

	#internal functions that are specifically written to cope with unstranded RNA-Seq data
	getAllSubstNoStrand <- function( sortedBam, cores ) {
		emd <- elementMetadata( sortedBam )		
		cov <- coverage( sortedBam )
		mds <- emd$MD
		flag <- grepl( pattern = '[A-Za-z]', mds ) 
		hasSubst <- sortedBam[ flag ]
		substMtx <- split( hasSubst, hasSubst$MD )	#suggested by MM, generates a GRangesList

		processedMD <- processMD( substMtx, cores )
		processedMD <- getComplSubst( processedMD )
		n <- nrow( processedMD )
		pos <- as.numeric( processedMD[, 'posMismatchesGenome'] )
		allSubst <- GRanges( seqnames      = processedMD[, 'seqnames'], 
				     ranges        = IRanges( pos , pos ), 
	     			     strand        = processedMD[, 'strand'], 
			     	     substitutions = processedMD[, 'substitutions'], 
                                     coverage      = rep( Inf, n ) )

		return( allSubst )
	}

	getCountTableRNASeq <- function( subst ) {
		filteredDf <- as.data.frame( subst )
		tmp <- filteredDf[, -6]
		remove <- duplicated( tmp )
		filteredDf <- filteredDf[ !remove, ]
		countTable <- GRanges( seqnames      = filteredDf[, 'seqnames'], 
				       ranges        = IRanges( filteredDf[, 'start'], filteredDf[, 'end'] ),
 	       			       strand        = filteredDf[, 'strand'],
 		      	       	       substitutions = filteredDf[, 'substitutions'], 
				       coverage      = filteredDf[, 'coverage'] ) 
		counts <- countOverlaps( countTable, subst )
		elementMetadata( countTable )[, 'count'] <- counts
		return( countTable )		
	}

#	getAllSubstNoStrand <- function( sortedBam, cores ) {
#		emd <- elementMetadata( sortedBam )		
#		rsh <- function( x ) 
#			matrix( x, ncol = 4 )
#		cov <- coverage( sortedBam )
#		mds <- emd$MD
#		flag <- grepl( pattern = '[A-Za-z]', mds ) 
#		hasSubst <- sortedBam[ flag ]
#		substMtx <- as.matrix( as.data.frame( hasSubst )[, c( 'seqnames', 'start', 'strand', 'qseq', 'MD' )] )
#		substMtx <- split( substMtx[, c( 'seqnames', 'start', 'strand', 'qseq' )], substMtx[, 'MD'] ) 
#		substMtx <- lapply( substMtx, rsh )
#		processedMD <- processMD( substMtx, cores )
#		processedMD <- getComplSubst( processedMD )
#		n <- nrow( processedMD )
#		pos <- as.numeric( processedMD[, 'posMismatchesGenome'] )
#		allSubst <- GRanges( seqnames      = processedMD[, 'chromosomes'], 
#				     ranges        = IRanges( pos , pos ), 
#	     			     strand        = processedMD[, 'strands'], 
#			     	     substitutions = processedMD[, 'substitutions'], 
#                                    coverage      = rep( Inf, n ) )
#		return( allSubst )
#	}
#	getCountTableRNASeq <- function( subst ) {
#		filteredDf <- as.data.frame( subst )
#		tmp <- filteredDf[, -6]
#		remove <- duplicated( tmp )
#		filteredDf <- filteredDf[ !remove, ]
#		countTable <- GRanges( seqnames      = filteredDf[, 'seqnames'], 
#				       ranges        = IRanges( filteredDf[, 'start'], filteredDf[, 'end'] ),
#	       			       strand        = filteredDf[, 'strand'],
#		      	       	       substitutions = filteredDf[, 'substitutions'], 
#				       coverage      = filteredDf[, 'coverage'] ) 
#		counts <- countOverlaps( countTable, subst )
#		elementMetadata( countTable )[, 'count'] <- counts
#		return( countTable )		
#	}

	getFDRBounds <- function( GR, span = 0.1 ) {
		intervals <- seq( 0, 1, by = span )
		emd <- elementMetadata( GR )
		n <- length( intervals ) - 1
		FDR <- list( upper = rep( Inf, n ), lower = rep( Inf, n ) )
		ps <- rep( Inf, n )
		for(i in seq_len( n ) ) {
			subset <- (intervals[ i ] < emd[, 'rsf'] & emd[, 'rsf'] <= intervals[ i + 1 ])
			emdSubset <- emd[ subset, 'rsfRNA']
			P <- length( emdSubset )
			U <- sum( emdSubset > 0 ) #no matter what RSF value
			L <- sum( intervals[ i ] < emdSubset & emdSubset <= intervals[ i + 1 ])
			FDR$upper[ i ] <- U / P
			FDR$lower[ i ] <- L / P
			ps[ i ] <- P
		}
		return( list(upper = FDR$upper, lower = FDR$lower, ps = ps) )
	}
	
	#1-select all PAR-CLIP substitutions of type substitution
	allSub <- countTable[ elementMetadata( countTable )[, 'substitutions'] == substitution ]

	#2-retain RNA-Seq reads overlapping (unstranded query) to PAR-CLIP substitutions and compute RNA-Seq coverage
	if( verbose )
		message( 'Extracting RNA-Seq reads overlapping with PAR-CLIP transitions...', '\n' )
	RNASeqSubset <- subsetByOverlaps( RNASeq, allSub, ignore.strand = TRUE )
	covRNASeq <- coverage( RNASeqSubset )

	#3-compute all substitutions observed in RNA-Seq data
	if( verbose )
		message( 'Identifying substitutions in RNA-Seq reads...' )
	RNASeqSubsetNoStrand <- RNASeqSubset
	strand(RNASeqSubsetNoStrand) <- '*'
	AllSubstRNASeq <- getAllSubstNoStrand( RNASeqSubsetNoStrand, cores ) 

	#4-retain only substitutions of same type as PAR-CLIP substitution (and its complement)
	substitutionCompl <- as.character( complement( DNAString( substitution ) ) )
	substitutions <- elementMetadata(AllSubstRNASeq)[, 'substitutions']
	AllSubstRNASeq <- AllSubstRNASeq[ substitutions  == substitution | substitutions  == substitutionCompl ]

	#5-get count table for these RNA-Seq substitutions
	if( verbose )
		message( 'Getting count table for RNA-Seq substitutions and their coverage...' )
	countTableRNASeq <- getCountTableRNASeq( AllSubstRNASeq )

	#6-compute RNA-Seq coverage at PAR-CLIP substitutions
	counts <- countOverlaps( allSub, RNASeqSubset, ignore.strand = TRUE )
	elementMetadata( allSub )[, 'coverageRNA'] <- counts

	#7-retain only those PAR-CLIP substitutions with RNA-Seq coverage >= minCov
	allSub <- allSub[ counts >= minCov ]

	#8-compute the number of transitions per position observed in RNASeq data
	if( verbose )
		message( 'Summarizing RNA-Seq substitutions...' )
	olaps <- findOverlaps( allSub, countTableRNASeq, ignore.strand = TRUE )
	elementMetadata( allSub )[, 'countRNA'] <- 0
	elementMetadata( allSub )[queryHits( olaps ), 'countRNA'] <- elementMetadata( countTableRNASeq[ subjectHits( olaps ) ] )[, 'count']
	
	#9-compute RSFs
	if( verbose )
		message( 'Computing RSFs and FDR...' )
	elementMetadata( allSub )[, 'rsf'] <- Inf
	elementMetadata( allSub )[, 'rsfRNA'] <- Inf
	emd <- elementMetadata( allSub )
	elementMetadata( allSub )[, 'rsf'] <- emd[, 'count'] / emd[, 'coverage']
	elementMetadata( allSub )[, 'rsfRNA'] <- emd[, 'countRNA'] / emd[, 'coverageRNA']

	#10-compute FDR
	FDR <- getFDRBounds( allSub, span )
	FDR$ps <- FDR$ps / sum( FDR$ps )
	
	if( plot ) {
		x <- 0 : ( length( FDR$upper ) - 1)
		lab <- paste0( '(', span * x, ',', span * ( x + 1 ), ']' ) 
		plot( FDR$upper, type = 'l', 
	       		lwd  =  2, 
	        	col  = 'skyblue3', 
	        	xaxt = 'n', 
	        	xlab = '', 
	        	ylab = 'FDR',	
			ylim = range( FDR ),
			...)
		lines( FDR$lower, lwd = 2, lty = 2, col = 'skyblue3' )
		lines( FDR$ps, lwd = 2, col = 'red3' )

		legend( 'topright', 
		 	legend = c( 'Upper bound', 'Lower bound', 'Positives' ), 
		 	col = c( 'skyblue3', 'skyblue3', 'red3' ), 
		 	lwd = 2, 
			lty = c( 1, 2, 1 ),
		 	bty = 'n' )
		axis( 1, at = x + 1, labels = lab, las = 2 )
	}
	
	return( FDR )
}
