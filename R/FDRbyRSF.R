#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#2013 - Federico Comoglio & Cem Sievers, D-BSSE, ETH Zurich
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

estimateFDR <- function( countTable, RNASeq, substitution = 'TC', minCov = 20, span = 0.1, cores = 1, plot = TRUE, verbose = TRUE, ... ) {
# estimate False Discovery Rate within the RSF support using RNA-Seq data (current version requires unstranded RNA-Seq)
#
# Args:
#   countTable:  a GRanges object, corresponding to a count table where each substitution has a corresponding strand-specific coverage and a count value, as returned by the getAllSub function
#   RNASeq: a GRanges object containing all RNA-Seq reads, as obtained by loading a BAM file using the readSortedBam function
#   substitution: character, the substitution/transition induced by the experimental procedure
#   minCov: numeric, the minimum required RNA-Seq coverage for a substitution identified by PAR-CLIP to be considered. Default is 20.
#   span: numeric, the width of RSF intervals to be considered
#   cores: numeric, the number of cores to be used for parallel computing. Default is 1
#   plot: logical, if TRUE plots line plot with additional features. Default is TRUE
#   verbose: logical, if TRUE, prints steps. Default is TRUE
#   ...: additional parameters to be passed to the plot function
#
# Returns:
#   A list with two slots, containing the total number of positives and the FDR for each RSF interval, respectively. If plot, then also plots positives and FDR as a line plot.
#
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

		processedMD <- wavClusteR:::processMD( substMtx, cores )
		processedMD <- wavClusteR:::getComplSubst( processedMD )
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
