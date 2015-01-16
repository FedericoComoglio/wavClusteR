#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#2013 - Federico Comoglio & Cem Sievers, D-BSSE, ETH Zurich
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



#' Load a sorted BAM file
#' 
#' Load a sorted BAM file. Optionally, only reads mapping to a specific set of
#' genomics coordinates are loaded. Only fields strictly necessary to run a
#' wavClusteR analysis are loaded.
#' 
#' 
#' @usage readSortedBam(filename, which)
#' @param filename Name of the sorted BAM file, including full path to file if
#' it is located outside the current working directory.
#' @param which a GRanges, RangesList or RangedData specifying the regions on
#' the reference sequence for which matches are desired. See the documentation
#' of the \code{Rsamtools} package for details.
#' @return a GRanges object containing aligned reads, including read sequence
#' (qseq) and MD tag (MD)
#' @note The input BAM file must be sorted and indexed. Alignment with bowtie
#' or bowtie2, conversion from SAM to BAM output, sorting and indexing using
#' SAMtools is recommended.
#' @author Federico Comoglio
#' @references Martin Morgan and Herve Pages, Rsamtools: Binary alignment
#' (BAM), variant call (BCF), or tabix file import,
#' \url{http://bioconductor.org/packages/release/bioc/html/Rsamtools.html}
#' @keywords Import
#' @examples
#' 
#' library(Rsamtools)
#' filename <- system.file( "extdata", "example.bam", package = "wavClusteR" )
#' sortedBam <- readSortedBam( filename = filename )
#' 
#' @export readSortedBam
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
