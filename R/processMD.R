#2013 - Federico Comoglio & Cem Sievers, D-BSSE, ETH Zurich

processMD <- function( SubstMtx, cores ) {
# higher-order wrapper to extract information from MD tag
#
# Args:
#   SubstMtx: a list of matrices, where each list entry is named by a specific MD tag value (as results from splitting by MD tag) and corresponds to a N x 4 matrix, where N is the number of reads showing that MD tag and 4 corresponds to seqnames, start, strand and qseq of the corresponding reads.
#   cores: numeric, the number of cores to be used for parallel computing
#
# Returns:
#   a matrix, containing all identified substitutions from the MD tag
#
# Error handling
#   ...
 	#backend for doMC (windows)
	if ( suppressWarnings( require( "doMC" ) ) ) {
		registerDoMC( cores = cores )
     	} else {
     	}
	
	n <- length( SubstMtx )
	set <- c( seq( 1, n, 1e3 ), n + 1 ) #construct chunks of size 1e3
	m <- length( set )

	message( 'Considering substitutions, n = ', n, ', processing in ', m - 1, ' chunks' )
	k = NULL; rm( k ) #pass Rcheck
	substitutions <- foreach( k = seq_len( m - 1 ), .combine = rbind) %dopar% {
		message('   chunk #: ', k)
		chunk <- SubstMtx[ set[ k ] : ( set[ k + 1 ] - 1 ) ]
		processChunk( chunk )
	}

	return( substitutions )
}


#processChunk <- function( chunk ) {
# process a single chunk of MD tags (size 1e3)
#
# Args:
#   chunk: a subset of SubstMtx as produced and passed by the processMD function 
#
# Returns:
#   a matrix, containing identified substitutions according to MD tag for the processed chunk of data
#
# Error handling
#   ...
#	n <- length( chunk )
#	mds <- names( chunk ) #contains all MD tags of the chunk being processed
#
#	chunkSubstitutions <- c() 
#
#	for( i in seq_len( n ) ) {
#		md <- mds[ i ]
#		chunkSubstitutions <- rbind( chunkSubstitutions, extractSingleMD( md, chunk[[i]] ) )
#	}
#
#	return( chunkSubstitutions )
#}


processChunk <- function( chunk ) { #contributed by Martin Morgan
# process a single chunk of MD tags (size 1e3)
#
# Args:
#   chunk: a subset of SubstMtx as produced and passed by the processMD function 
#
# Returns:
#   a data frame, containing identified substitutions according to MD tag for the processed chunk of data
#
# Error handling
#   ...

     mds <- names( chunk )
	strands <- unname( strand( chunk ) )	
	strands <- lapply( strands, as.character )

     ## Match a nucleotide preceed by a digit (returns list)
     mdsSplit <- regmatches( mds, gregexpr( '\\d*|[ACGTN]{1}', mds) )
     len <- elementLengths( mdsSplit )
     grp <- rep( seq_along( len ), len )

     ## number and geometry of substitutions
     u <- unlist(mdsSplit)
     isChar <- grepl("^[^[:digit:]]+$", u)
     isCharGrp <- splitAsList( isChar, grp )

     mismatchedBases <- splitAsList( u[isChar], grp[isChar] )

     mismatchedPos <- integer(length(u))
     mismatchedPos[isChar] <- 1L
     mismatchedPos[!isChar] <- as.integer(u[!isChar])
     mismatchedPos <- cumsum(splitAsList(mismatchedPos, grp))[isCharGrp]

     chunk_n <- unname(elementLengths(chunk))
     mismatched_n <- unname(elementLengths(mismatchedPos))

     ## nucleotide substitutions
     at0 <- mismatchedPos[rep(seq_along(chunk), chunk_n)]
     qseq <- unlist(chunk)$qseq
     offset <- c(0, cumsum(width(qseq)[-length(qseq)])) + at0
     at <- IRanges( unlist( offset ), width = 1 )

     nucleotides <- extractAt(unlist(qseq), at)
     mismatchedBases <- unlist(rep(unname(mismatchedBases), chunk_n))
#     subst <- DNAStringSet(paste0(nucleotides, mismatchedBases))
     subst <- paste0( mismatchedBases, nucleotides )

     ##4-Identify genomic positions of mismatches
     posMismatchesGenome <- unlist(start(chunk)) + at0 - 1

	##5-Get strand information right
	strands <- unlist( lapply( seq_len( length( mismatched_n ) ), 
		function( i ) rep( strands[[ i ]], each = mismatched_n[ i ] ) ) )

     ##6-Prepare output
     DataFrame(
         chunkIdx = rep(seq_along(chunk), mismatched_n * chunk_n),
         readIdx  = rep(seq_along(qseq), rep(mismatched_n, chunk_n)),
         seqnames = unlist(unname(rep(seqnames(chunk), mismatched_n))),
         posMismatchesGenome = unlist(unname(posMismatchesGenome)),
         strand   = strands,
         substitutions = subst)
  }


##backup version (wrong strand assignment for complex MD)
#processChunk <- function( chunk ) { #contributed by Martin Morgan
## process a single chunk of MD tags (size 1e3)
##
## Args:
##   chunk: a subset of SubstMtx as produced and passed by the processMD function 
##
## Returns:
##   a data frame, containing identified substitutions according to MD tag for the processed chunk of data
##
## Error handling
##   ...
#
#     mds <- names( chunk )
#
#     ## Match a nucleotide preceed by a digit (returns list)
#     mdsSplit <- regmatches( mds, gregexpr( '\\d*|[ACGTN]{1}', mds) )
#     len <- elementLengths( mdsSplit )
#     grp <- rep( seq_along( len ), len )
#
#     ## number and geometry of substitutions
#     u <- unlist(mdsSplit)
#     isChar <- grepl("^[^[:digit:]]+$", u)
#     isCharGrp <- splitAsList( isChar, grp )
#
#     mismatchedBases <- splitAsList( u[isChar], grp[isChar] )
#
#     mismatchedPos <- integer(length(u))
#     mismatchedPos[isChar] <- 1L
#     mismatchedPos[!isChar] <- as.integer(u[!isChar])
#     mismatchedPos <- cumsum(splitAsList(mismatchedPos, grp))[isCharGrp]
#
#     chunk_n <- unname(elementLengths(chunk))
#     mismatched_n <- unname(elementLengths(mismatchedPos))
#
#     ## nucleotide substitutions
#     at0 <- mismatchedPos[rep(seq_along(chunk), chunk_n)]
#     qseq <- unlist(chunk)$qseq
#     offset <- c(0, cumsum(width(qseq)[-length(qseq)])) + at0
#     at <- IRanges( unlist( offset ), width = 1 )
#
#     nucleotides <- extractAt(unlist(qseq), at)
#     mismatchedBases <- unlist(rep(unname(mismatchedBases), chunk_n))
##     subst <- DNAStringSet(paste0(nucleotides, mismatchedBases))
#     subst <- paste0( mismatchedBases, nucleotides )
#
#     ##4-Identify genomic positions of mismatches
#     posMismatchesGenome <- unlist(start(chunk)) + at0 - 1
#
#     ##5-Prepare output
#     DataFrame(
#         chunkIdx = rep(seq_along(chunk), mismatched_n * chunk_n),
#         readIdx  = rep(seq_along(qseq), rep(mismatched_n, chunk_n)),
#         seqnames = unlist(unname(rep(seqnames(chunk), mismatched_n))),
#         posMismatchesGenome = unlist(unname(posMismatchesGenome)),
#         strand   = unlist(unname(rep(strand(chunk), mismatched_n))),
#         substitutions = subst)
#  }


#extractSingleMD <- function( MD, SubstMtxInstance ) {  
# Extract information from MD tag of a single instance of the SubstMtx matrix. 
#
# Args:
#   MD: character, the MD tag being considered
#   SubstMtxInstance: a N x 4 matrix, where N is the number of reads showing MD and 4 corresponds to seqnames, start, strand and qseq of the corresponding reads.
#
# Returns:
#   a matrix, where each row corresponds to the genomic position (chromosome, position, strand) of a given identified substitution
#
# Error handling
#   ...

	#1-Extract vectors 
#	chr <- SubstMtxInstance[, 1]
#	start <- as.numeric( SubstMtxInstance[, 2] )
#	strand <- SubstMtxInstance[, 3]
#	qseq <- SubstMtxInstance[, 4]

	#2-Match a nucleotide preceed by a digit (returns list)
#	MDSplit <- regmatches( MD, gregexpr( '\\d*|[ACGTN]{1}', MD) )[[ 1 ]]
#	MDSplitNum <-  suppressWarnings( as.numeric( MDSplit ) )
#	chars <- is.na( MDSplitNum )	#vector, position of characters
#	mismatchedBases <- MDSplit[ chars ]
#	MDSplitNum[ which( is.na( MDSplitNum ) ) ] <- 1 #positions with characters just contribute to 1 to the cumulative sum
#	posMismatches <- cumsum( MDSplitNum )[ chars ]
#	n <- length( posMismatches )

	#3-Extract characters from sequence
#	readSplit <- strsplit( qseq, '' )
#	substitutions <- vapply( readSplit, FUN.VALUE = rep( '', n ), function(x) x[ posMismatches ] )
#	substitutions <- paste( mismatchedBases, substitutions, sep = '' ) 

	#4-Identify genomic positions of mismatches
#	posMismatchesGenome <- as.vector( sapply( start, function(x) x + posMismatches - 1 ) )

	#5-Prepare output
#	chromosomes <- rep( chr, each = n )
#	strands <- rep( strand, each = n )
#	extractedMDInfo <- cbind( chromosomes, posMismatchesGenome, strands, substitutions )

#	return( extractedMDInfo )
#}
