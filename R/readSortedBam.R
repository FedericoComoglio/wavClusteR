#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#2013 - Federico Comoglio & Cem Sievers, D-BSSE, ETH Zurich
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

readSortedBam <- function( filename, which ) {
# read sorted BAM file. Optionally, only reads mapping to a specific set of genomics coordinates are read.
#
# Args:
#   filename: character, the filename and/or path to it
#   which: a GRanges, RangesList or RangedData specifying the regions on the reference sequence for which matches are desired.
#
# Returns:
#   a GRanges object containing all raw data for wavClusteR2 analysis
#
# Error handling
#   ...
	pos <- NULL; rm( pos ) #pass check

	tag <- c('MD') #MD: mismatch field
	if( missing(which) ) {	
		param <- ScanBamParam( flag        = scanBamFlag(isUnmappedQuery=FALSE),
				       tag         = tag,
			       	       simpleCigar = FALSE,
			               what        = c('rname', 'pos', 'qwidth', 'strand', 'seq') )
	}
	else {
		param <- ScanBamParam( flag        = scanBamFlag(isUnmappedQuery=FALSE),
				       tag         = tag,
				       simpleCigar = FALSE,
		                       what        = c('rname', 'pos', 'qwidth', 'strand', 'seq'),
		   	               which       = which )

	}
	rawData <- scanBam( filename, param = param )
	#coerce to GRanges and return object.

	sortedBam <- lapply( rawData, function(rawDataChr) with( rawDataChr, GRanges( seqnames = rname, 
			                                                              ranges   = IRanges( start = pos, width = qwidth ),
		                                                                      strand   = strand,
			                                                              qseq     = seq,
			                                                              tag      = tag ) ) )
	sortedBam <- GRangesList( sortedBam )
	sortedBam <- unlist( sortedBam, use.names = FALSE ) #this routine ensures that all relevant genomic regions are considered

	return( sortedBam )
}
