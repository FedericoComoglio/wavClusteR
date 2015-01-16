#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#2013 - Federico Comoglio & Cem Sievers, D-BSSE, ETH Zurich
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

exportGR <- function( GR, filename, trackname, description) {
# Create UCSC compatible BED tracks of input object
#
# Args:
#   GR: GRanges object as returned by the getClusters or the filterClusters function, or as returned by getHighConfSub
#   filename: character, name of the file to be written.
#   trackname: character, track name. 
#   description: character, track description.
#
# Returns:
#   writes a BED file to the specified path
#
# Error handling
#   ...

	n <- length( GR )
	chrs <- as.character( seqnames( GR ) )
	starts <- start( GR )
	ends <- end( GR )
	strands <- as.character( strand( GR ) )	
	
	#1-prepare header and cluster names
	header <- paste( 'track name = \'', filename, '\' description = \'', filename, '\' visibility=2 itemRgb=\'On\'', sep = '' )
	names <- paste( chrs, starts, ends, sep = '-' )

	#2-write the header to output file
	write.table( header, 
		     file      = filename, 
		     row.names = FALSE, 
		     col.names = FALSE, 
		     quote     = FALSE )

	#3-strand-specific color coding	
	reps <- as.numeric( table( strands ) )
	colSet <- rep( c( '255,0,0', '0,0,255' ), times = reps )	

	#4-prepare and write output
	bed <- cbind( chrs, starts, ends, names, rep( 0, n ), strands, starts, ends, colSet )
	write.table( bed, 
		     file      = filename, 
		     sep       = '\t', 
		     row.names = FALSE, 
		     col.names = FALSE, 
		     quote     = FALSE, 
		     append    = TRUE ) 
}




#' Export coverage as BigWig track
#' 
#' Export coverage as BigWig track, compatible with the UCSC genome browser
#' 
#' 
#' @usage exportCoverage(coverage, filename = 'wavClusters.BigWig')
#' @param coverage An Rle object containing the coverage at each genomic
#' position as returned by a call to \code{coverage}
#' @param filename A character defining the BED file name. Default to
#' "wavClusters.BigWig"
#' @return A BigWig file of the exported Rle object
#' @author Federico Comoglio
#' @keywords postprocessing
#' @export exportCoverage
exportCoverage <- function( coverage, filename = 'wavClusters.BigWig' ) {
# Create UCSC compatible BigWig tracks of the coverage
#
# Args:
#   coverage: Rle object, containing the coverage of the sortedBam file as returned by the coverage function
#   filename: character, name of the file to be written. Default is 'wavClusters.bigWig'
#
# Returns:
#   writes a BigWig file to the specified path
#
# Error handling
#   ...
	export( coverage, con = filename, format = 'BigWig' )
}




#' Export clusters as BED track
#' 
#' Export clusters as BED track, compatible with the UCSC genome browser
#' 
#' 
#' @usage exportClusters(clusters, filename = 'wavClusters.bed', trackname =
#' 'wavClusters', description = 'wavClusters')
#' @param clusters GRanges object containing individual clusters as identified
#' by the \link{filterClusters} function
#' @param filename A character defining the BED file name. Default to
#' "wavClusters.bed"
#' @param trackname A character defining the track.name of the BED file.
#' Default to "wavClusters"
#' @param description A character defining the description of the BED file.
#' Default to "wavClusters"
#' @return A BED file of the exported GRanges object
#' @note Clusters are color coded according to their strand information (red
#' for the plus strand, blue for the minus strand).
#' @author Federico Comoglio
#' @seealso \code{\link{filterClusters}}
#' @keywords postprocessing
#' @export exportClusters
exportClusters <- function( clusters, filename = 'wavClusters.bed', trackname = 'wavClusters', description = 'wavClusters') {
# Create UCSC compatible BED tracks of clusters
#
# Args:
#   clusters: GRanges object as returned by the getClusters or the filterClusters function
#   filename: character, name of the file to be written. Default is 'wavClusters.bed'
#   trackname: character, track name. Default is 'wavClusters'
#   description: character, track description. Default is 'wavClusters'
#
# Returns:
#   writes a BED file to the specified path
#
# Error handling
#   ...
	exportGR( clusters, filename = filename, trackname = trackname, description = description )
}




#' Export high-confidence substitutions as BED track
#' 
#' Export high-confidence substitutions as BED track, compatible with the UCSC
#' genome browser
#' 
#' 
#' @usage exportHighConfSub(highConfSub, filename = 'highConfSub.bed',
#' trackname = 'highConfSub', description = 'highConfSub')
#' @param highConfSub GRanges object containing high-confidence substitution
#' sites as returned by the \link{getHighConfSub} function
#' @param filename A character defining the BED file name. Default to
#' "wavClusters.bed"
#' @param trackname A character defining the track.name of the BED file.
#' Default to "wavClusters"
#' @param description A character defining the description of the BED file.
#' Default to "wavClusters"
#' @return A BED file of the exported GRanges object
#' @note Substitutions are color coded according to their strand information
#' (red for the plus strand, blue for the minus strand).
#' @author Federico Comoglio
#' @seealso \code{\link{getHighConfSub}}
#' @keywords postprocessing
#' @export exportHighConfSub
exportHighConfSub <- function( highConfSub, filename = 'highConfSub.bed', trackname = 'highConfSub', description = 'highConfSub') {
# Create UCSC compatible BED tracks of high confidence substitutions
#
# Args:
#   highConfSub: GRanges object as returned by the getHighConfSub function
#   filename: character, name of the file to be written. Default is 'highConfSub.bed'
#   trackname: character, track name. Default is 'highConfSub'
#   description: character, track description. Default is 'highConfSub'
#
# Returns:
#   writes a BED file to the specified path
#
# Error handling
#   ...
	exportGR( highConfSub, filename = filename, trackname = trackname, description = description )
}




#' Export cluster sequences for motif search analysis
#' 
#' Export cluster sequences for motif search analysis (FASTA format), e.g.
#' using MEME-ChIP
#' 
#' 
#' @usage exportSequences(clusters, filename = 'wavClusters.fasta')
#' @param clusters GRanges object containing individual clusters as identified
#' by the \link{filterClusters} function
#' @param filename A character defining the BED file name. Default to
#' "wavClusters.fasta"
#' @return A FASTA file containing the cluster sequences
#' @author Federico Comoglio
#' @seealso \code{\link{filterClusters}}
#' @keywords postprocessing
#' @export exportSequences
exportSequences <- function( clusters, filename = 'wavClusters.fasta' ) {
# Export DNA sequences underlying clusters for motif search
#
# Args:
#   clusters: GRanges object as returned by the filterClusters function
#   filename: character, name of the file to be written. Default is 'wavClusters.fasta'
#
# Returns:
#   writes a fasta file to the specified path, to be used for motif search (e.g. using the web GUI of MEME-ChIP)
#
# Error handling
#   ...
	sequences <- elementMetadata( clusters )[, 'Sequence']
	sequences <- as.character( sequences )
	sequences <- as.list( sequences )
	names <- 1 : length( sequences )	

	write.fasta( names     = names, 
		     sequences = sequences, 
		     file.out  = filename )
}
