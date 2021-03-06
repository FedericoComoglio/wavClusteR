#2013 - Federico Comoglio & Cem Sievers, D-BSSE, ETH Zurich

#' Compute and plot read densities in genomic regions around transcription
#' start sites
#' 
#' Aligned reads are used generate a metaTSS profile across genomic regions
#' containing transcription start sites (TSSs).
#' 
#' 
#' @usage getMetaTSS(sortedBam, txDB = NULL, upstream = 1e3, downstream = 1e3,
#' nBins = 40, genome = 'hg19', tablename = 'ensGene', unique = FALSE, plot =
#' TRUE, verbose = TRUE, ...)
#' @param sortedBam GRanges object containing aligned reads as returned by
#' \link{readSortedBam}
#' @param txDB TranscriptDb object obtained through a call to the
#' \code{makeTxDbFromUCSC} function in the \code{GenomicFeatures}
#' package. Default is NULL, namely the object will be fetched internally
#' @param upstream An integer corresponding to the number of bases to be
#' considered upstream the TSS. Default is 1000
#' @param downstream An integer corresponding to the number of bases to be
#' considered downstream the TSS. Default is 1000
#' @param nBins An integer corresponding to the number of bins to be used to
#' partition the genes. Default is 40, i.e. bin size 50 bases
#' @param genome A character specifying the genome abbreviation used by UCSC.
#' Available abbreviations are returned by a call to \code{ucscGenomes()[ ,
#' "db"]}. Default is "hg19" (human genome)
#' @param tablename A character specifying the name of the UCSC table
#' containing the transcript annotations to retrieve. Available table names are
#' returned by a call to \code{supportedUCSCtables()}. Default is "ensGene",
#' namely ensembl gene annotations
#' @param unique Logical, if TRUE only non-overlapping TSSs extended by
#' upstream/downstream are considered. Default is FALSE, i.e. all TSSs are
#' considered
#' @param plot Logical, if TRUE a dotchart with cluster annotations is produced
#' @param verbose Logical, if TRUE processing steps are printed
#' @param ... Additional parameters to be passed to the \code{plot} function
#' @return A numeric vector of the same length as \code{nBins} containing
#' normalized counts. If \code{plot}, the metaTSS profile is also depicted as a
#' line plot.
#' @author Federico Comoglio
#' @seealso \code{\link{readSortedBam}}
#' @references Comoglio F*, Sievers C* and Paro R, wavClusteR: an R package for
#' PAR-CLIP data analysis, submitted
#' @keywords postprocessing graphics
#' @examples
#' 
#' require(BSgenome.Hsapiens.UCSC.hg19)
#' 
#' filename <- system.file( "extdata", "example.bam", package = "wavClusteR" )
#' example <- readSortedBam( filename = filename )
#' \dontrun{tss <- getMetaTSS( sortedBam = example )}
#' 
#' @export getMetaTSS

getMetaTSS <- function( sortedBam, txDB = NULL, upstream = 1e3, downstream = 1e3, nBins = 40, genome = 'hg19', tablename = 'ensGene', unique = FALSE, plot = TRUE, verbose = TRUE, ... ) {
# Returns:
#   a numeric vector of the same length as nBins with normalized counts. If plot, then also plots the distribution
#
# Error handling
#   if txDB not provided, check that genome and tablename are within those available from UCSC, otherwise raise an error
	availGenomes <- ucscGenomes()[ , 'db']
	availTables <- rownames( supportedUCSCtables() )
	if( is.null( txDB ) & (!(genome %in% availGenomes) | (!(tablename %in% availTables))) ) {
		stop('transcriptDB object not provided and genome or tablename not within 
		      those available from UCSC. Please use ucscGenomes()[ , \'db\'] and supportedUCSCtables() 
                      for a list of supported ones.')
	}

	#1-if not provided, create the TranscriptDb object
	if( is.null( txDB ) ) {
		if(verbose)
			message( 'Creating TranscriptDb object...' )
		txDB <- makeTxDbFromUCSC(genome = genome, tablename = tablename)
	}

	#2-extract genes from annotations
	gn <- genes( txDB )
	
	#3-compute regions around TSS by extending up/downstream
	if(verbose)
		message( 'Extracting promoter regions...' )
	suppressWarnings( tss <- promoters( gn, upstream = upstream, downstream = downstream ) )

	#4-if unique, filter tss for non-overlapping instances
	if( unique ) {
		if(verbose)
			message( 'Retaining non-overlapping promoters...' )
		olaps <- findOverlaps( tss, drop.self = TRUE )
		overlapping <- unique( c( queryHits( olaps ), subjectHits( olaps ) ) )
		tss <- tss[ -overlapping ]
	}
	
	#5-construct bins
	chrs <- seqnames( tss )
	starts <- start( tss )
	ends <- end( tss )
	strands <- strand( tss )
	minus <- which( strands == '-' )
	n <- length( starts )
	
	if(verbose)
		message( 'Computing MetaTSS distribution...' )

	posList <- lapply( 1 : n, function( i ) round( seq( starts[ i ], ends[ i ], length = nBins + 1 ) ) )
	posList <- lapply( posList, function( x ) cbind( x[ -length( x ) ], x[ -1 ] ) )
	x <- lapply(posList, nrow )
	posMatrix <- do.call( rbind, posList )

	bins <- GRanges( seqnames = rep( chrs, each = nBins ), 
			 ranges    = IRanges( posMatrix[, 1], posMatrix[, 2] ), 
			 strand    = rep( strands, each = nBins ) )

	#6-compute overlaps (strand-specific) and reshape output
	counts <- countOverlaps( bins, sortedBam  )
	counts <- matrix( counts, ncol = nBins, byrow = TRUE )
	
	#7-reverse values for TSS on the - strand
	counts[minus, ] <- t( apply( counts[minus, ], 1, rev ) )
	
	#8-get proper distribution
	metaProfile <- colSums( counts ) / sum( counts )

	if( plot ) { 
		plot( metaProfile, 
		      type = 'l', 
 		      xaxt = 'n', 
		      xlab = 'Position w.r.t TSS (bp)',
                      ylab = 'Normalized read count', 
		      lwd  = 2, 
		      col  = 'skyblue3',
		      ... )
		grid()
		axis( 1, at = c( 1, nBins / 2, nBins ), labels = c( -upstream, 'TSS', downstream ) )
	}

	return( metaProfile )
}




#' Compute and plot metagene profile using identified clusters
#' 
#' Transcriptome-wide identified clusters are used to generate a metagene
#' profile by binning gene bodies, upstream and downstream regions.
#' 
#' 
#' @usage getMetaGene(clusters, txDB = NULL, upstream = 1e3, downstream = 1e3,
#' nBins = 40, nBinsUD = 10, minLength = 1, genome = 'hg19', tablename =
#' 'ensGene', plot = TRUE, verbose = TRUE, ...)
#' @param clusters GRanges object containing individual clusters as identified
#' by the \link{getClusters} function
#' @param txDB TranscriptDb object obtained through a call to the
#' \code{makeTxDbFromUCSC} function in the \code{GenomicFeatures}
#' package. Default is NULL, namely the object will be fetched internally
#' @param upstream An integer corresponding to the number of bases to be
#' considered upstream the gene. Default is 1000
#' @param downstream An integer corresponding to the number of bases to be
#' considered downstream the gene. Default is 1000
#' @param nBins An integer corresponding to the number of bins to be used to
#' partition the genes. Default is 40
#' @param nBinsUD An integer corresponding to the number of bins to be used to
#' partion upstream and downstream regions. Defauls is 10, i.e. the bin size is
#' 100 bases for the default extension
#' @param minLength An integer indicating the the minimum required length of a
#' gene in order for it to be considered. Default is 1, i.e. all genes are
#' considered
#' @param genome A character specifying the genome abbreviation used by UCSC.
#' Available abbreviations are returned by a call to \code{ucscGenomes()[ ,
#' "db"]}. Default is "hg19" (human genome)
#' @param tablename A character specifying the name of the UCSC table
#' containing the transcript annotations to retrieve. Available table names are
#' returned by a call to \code{supportedUCSCtables()}. Default is "ensGene",
#' namely ensembl gene annotations
#' @param plot Logical, if TRUE a dotchart with cluster annotations is produced
#' @param verbose Logical, if TRUE processing steps are printed
#' @param ... Additional parameters to be passed to the \code{plot} function
#' @return A numeric vector of the same length as \code{nBins} + 2 *
#' \code{nBinsUD} containing normalized counts. If \code{plot}, the metagene
#' profile is also depicted as a line plot.
#' @author Federico Comoglio
#'
#' Comoglio F, Sievers C and Paro R (2015) Sensitive and highly resolved identification
#' of RNA-protein interaction sites in PAR-CLIP data, BMC Bioinformatics 16, 32.
#'
#' @seealso \code{\link{getClusters}}
#' @references Comoglio F*, Sievers C* and Paro R, wavClusteR: an R package for
#' PAR-CLIP data analysis, submitted
#' @keywords postprocessing graphics
#' @examples
#' 
#' require(BSgenome.Hsapiens.UCSC.hg19)
#' 
#' data( model, package = "wavClusteR" ) 
#' 
#' filename <- system.file( "extdata", "example.bam", package = "wavClusteR" )
#' example <- readSortedBam( filename = filename )
#' countTable <- getAllSub( example, minCov = 10, cores = 1 )
#' highConfSub <- getHighConfSub( countTable, supportStart = 0.2, supportEnd = 0.7, substitution = "TC" )
#' coverage <- coverage( example )
#' clusters <- getClusters( highConfSub = highConfSub, 
#'                          coverage = coverage, 
#'                          sortedBam = example, 
#' 	                 threshold = 2 ) 
#' 
#' fclusters <- filterClusters( clusters = clusters, 
#' 		             highConfSub = highConfSub, 
#'         		     coverage = coverage,
#' 			     model = model, 
#' 			     genome = Hsapiens, 
#' 		             refBase = 'T', 
#' 		             minWidth = 12 )
#' \dontrun{meta <- getMetaGene( clusters = fclusters )}
#' 
#' @export getMetaGene

getMetaGene <- function( clusters, txDB = NULL, upstream = 1e3, downstream = 1e3, nBins = 40, nBinsUD = 10, minLength = 1, genome = 'hg19', tablename = 'ensGene', plot = TRUE, verbose = TRUE, ... ) {
# Error handling
#   if txDB not provided, check that genome and tablename are within those available from UCSC, otherwise raise an error
	availGenomes <- ucscGenomes()[ , 'db']
	availTables <- rownames( supportedUCSCtables() )
	if( is.null( txDB ) & (!(genome %in% availGenomes) | (!(tablename %in% availTables))) ) {
		stop('transcriptDB object not provided and genome or tablename not within 
		      those available from UCSC. Please use ucscGenomes()[ , \'db\'] and supportedUCSCtables() 
                      for a list of supported ones.')
	}

	#1-if not provided, create the TranscriptDb object
	if( is.null( txDB ) ) {
		if(verbose)
			message( 'Creating TranscriptDb object...' )
		txDB <- makeTxDbFromUCSC(genome = genome, tablename = tablename)
	}

	#2-extract genes from annotations and retain only those passing the minLength cutoff
	if(verbose)
		message( 'Extracting genes and flanking regions...' )
	gn <- genes( txDB )
	gn <- gn[ width(gn) >= minLength ]

	#3-extract upstream/downstream region
	suppressWarnings( upstreamFlank <- flank( gn, upstream, start = TRUE ) )
	suppressWarnings( downstreamFlank <- flank( gn, downstream, start = FALSE ) )

	#4-construct bins
	chrs <- seqnames( gn )
	starts <- start( gn )
	ends <- end( gn )
	strands <- strand( gn )
	minus <- which( strands == '-' )
	n <- length( starts )
	
	if(verbose)
		message( 'Computing Metagene distribution...' )

	posList <- lapply( 1 : n, function( i ) round( seq( starts[ i ], ends[ i ], length = nBins + 1 ) ) )
	posList <- lapply( posList, function( x ) cbind( x[ -length( x ) ], x[ -1 ] ) )
	x <- lapply(posList, nrow )
	posMatrix <- do.call( rbind, posList )

	bins <- GRanges( seqnames = rep( chrs, each = nBins ), 
			 ranges    = IRanges( posMatrix[, 1], posMatrix[, 2] ), 
			 strand    = rep( strands, each = nBins ) )

	#5-compute overlaps (strand-specific),reshape output, and reverse values for genes on the - strand
	countsG <- countOverlaps( bins, clusters  )
	countsG <- matrix( countsG, ncol = nBins, byrow = TRUE )
	countsG[minus, ] <- t( apply( countsG[minus, ], 1, rev ) )

	#6-same procedure as above with upstream region
	chrs <- seqnames( upstreamFlank )
	starts <- start( upstreamFlank )
	ends <- end( upstreamFlank )
	
	posList <- lapply( 1 : n, function( i ) round( seq( starts[ i ], ends[ i ], length = nBinsUD + 1 ) ) )
	posList <- lapply( posList, function( x ) cbind( x[ -length( x ) ], x[ -1 ] ) )
	x <- lapply(posList, nrow )
	posMatrix <- do.call( rbind, posList )

	bins <- GRanges( seqnames = rep( chrs, each = nBinsUD ), 
			 ranges    = IRanges( posMatrix[, 1], posMatrix[, 2] ), 
			 strand    = rep( strands, each = nBinsUD ) )

	countsU <- countOverlaps( bins, clusters  )
	countsU <- matrix( countsU, ncol = nBinsUD, byrow = TRUE )
	countsU[minus, ] <- t( apply( countsU[minus, ], 1, rev ) )

	#8-same procedure with downstream region
	chrs <- seqnames( downstreamFlank )
	starts <- start( downstreamFlank )
	ends <- end( downstreamFlank )
	
	posList <- lapply( 1 : n, function( i ) round( seq( starts[ i ], ends[ i ], length = nBinsUD + 1 ) ) )
	posList <- lapply( posList, function( x ) cbind( x[ -length( x ) ], x[ -1 ] ) )
	x <- lapply(posList, nrow )
	posMatrix <- do.call( rbind, posList )

	bins <- GRanges( seqnames = rep( chrs, each = nBinsUD ), 
			 ranges    = IRanges( posMatrix[, 1], posMatrix[, 2] ), 
			 strand    = rep( strands, each = nBinsUD ) )

	countsD <- countOverlaps( bins, clusters  )
	countsD <- matrix( countsD, ncol = nBinsUD, byrow = TRUE )
	countsD[minus, ] <- t( apply( countsD[minus, ], 1, rev ) )

	#9-merge upstream-genes-downstream counts
	counts <- cbind( countsU, countsG, countsD )

	#10-get proper distribution
	metaProfile <- colSums( counts ) / sum( counts )

	if( plot ) { 
		plot( metaProfile, 
		      type = 'l', 
 		      xaxt = 'n',
		      xlab = 'Position',
                      ylab = 'Normalized number of wavClusters', 
		      lwd  = 2, 
		      col  = 'skyblue3',
		      ... )
		grid()
		axis( 1, at = c( 1, nBinsUD, nBins + nBinsUD, nBins + 2 * nBinsUD ), labels = c( -upstream, 'TSS', 'TES', paste0( '+', downstream ) ) )
		abline( v = c( nBinsUD, nBins + nBinsUD ), lwd = 1.5, lty = 2, col = 'gray50' )
	}

	return( metaProfile )
}



#' Compute and plot distribution of average coverage or relative log-odds as
#' metagene profile using identified clusters
#' 
#' Transcriptome-wide identified clusters are used to generate a metagene
#' profile by binning gene bodies. Within each bin, the distribution of the
#' average cluster coverage or of the relative log-odds is computed.
#' 
#' 
#' @usage getMetaCoverage(clusters, txDB = NULL, upstream = 1e3, downstream =
#' 1e3, nBins = 40, nBinsUD = 10, minLength = 1, genome = 'hg19', tablename =
#' 'ensGene', odds = FALSE, plot = TRUE, verbose = TRUE, ...)
#' @param clusters GRanges object containing individual clusters as identified
#' by the \link{getClusters} function
#' @param txDB TranscriptDb object obtained through a call to the
#' \code{makeTxDbFromUCSC} function in the \code{GenomicFeatures}
#' package. Default is NULL, namely the object will be fetched internally
#' @param upstream An integer corresponding to the number of bases to be
#' considered upstream the gene. Default is 1000
#' @param downstream An integer corresponding to the number of bases to be
#' considered downstream the gene. Default is 1000
#' @param nBins An integer corresponding to the number of bins to be used to
#' partition the genes. Default is 40
#' @param nBinsUD An integer corresponding to the number of bins to be used to
#' partion upstream and downstream regions. Defauls is 10, i.e. the bin size is
#' 100 bases for the default extension
#' @param minLength An integer indicating the the minimum required length of a
#' gene in order for it to be considered. Default is 1, i.e. all genes are
#' considered
#' @param genome A character specifying the genome abbreviation used by UCSC.
#' Available abbreviations are returned by a call to \code{ucscGenomes()[ ,
#' "db"]}. Default is "hg19" (human genome)
#' @param tablename A character specifying the name of the UCSC table
#' containing the transcript annotations to retrieve. Available table names are
#' returned by a call to \code{supportedUCSCtables()}. Default is "ensGene",
#' namely ensembl gene annotations
#' @param odds Logical, if TRUE relative log-odds distributions are shown
#' instead of mean coverage
#' @param plot Logical, if TRUE a dotchart with cluster annotations is produced
#' @param verbose Logical, if TRUE processing steps are printed
#' @param ... Additional parameters to be passed to the \code{plot} function
#' @return Called for its effects.
#' @author Federico Comoglio
#' @seealso \code{\link{getClusters}}
#' @references Comoglio F*, Sievers C* and Paro R, wavClusteR: an R package for
#' PAR-CLIP data analysis, submitted
#' @keywords postprocessing graphics
#' @examples
#' 
#' require(BSgenome.Hsapiens.UCSC.hg19)
#' 
#' data( model, package = "wavClusteR" ) 
#' 
#' filename <- system.file( "extdata", "example.bam", package = "wavClusteR" )
#' example <- readSortedBam( filename = filename )
#' countTable <- getAllSub( example, minCov = 10, cores = 1 )
#' highConfSub <- getHighConfSub( countTable, supportStart = 0.2, supportEnd = 0.7, substitution = "TC" )
#' coverage <- coverage( example )
#' clusters <- getClusters( highConfSub = highConfSub, 
#'                          coverage = coverage, 
#'                          sortedBam = example, 
#' 	                 threshold = 2 ) 
#' 
#' fclusters <- filterClusters( clusters = clusters, 
#' 		             highConfSub = highConfSub, 
#'         		     coverage = coverage,
#' 			     model = model, 
#' 			     genome = Hsapiens, 
#' 		             refBase = 'T', 
#' 		             minWidth = 12 )
#' \dontrun{getMetaCoverage( clusters = fclusters, odds = FALSE )}
#'
 
getMetaCoverage <- function( clusters, txDB = NULL, upstream = 1e3, downstream = 1e3, nBins = 40, nBinsUD = 10, minLength = 1, genome = 'hg19', tablename = 'ensGene', odds = FALSE, plot = TRUE, verbose = TRUE, ... ) {
# Error handling
#   if txDB not provided, check that genome and tablename are within those available from UCSC, otherwise raise an error
	availGenomes <- ucscGenomes()[ , 'db']
	availTables <- rownames( supportedUCSCtables() )
	if( is.null( txDB ) & (!(genome %in% availGenomes) | (!(tablename %in% availTables))) ) {
		stop('transcriptDB object not provided and genome or tablename not within 
		      those available from UCSC. Please use ucscGenomes()[ , \'db\'] and supportedUCSCtables() 
                      for a list of supported ones.')
	}

	#1-if not provided, create the TranscriptDb object
	if( is.null( txDB ) ) {
		if(verbose)
			message( 'Creating TranscriptDb object...' )
		txDB <- makeTxDbFromUCSC(genome = genome, tablename = tablename)
	}

	#2-extract genes from annotations and retain only those passing the minLength cutoff
	if(verbose)
		message( 'Extracting genes and flanking regions...' )
	gn <- genes( txDB )
	gn <- gn[ width(gn) >= minLength ]

	#3-extract upstream/downstream region
	suppressWarnings( upstreamFlank <- flank( gn, upstream, start = TRUE ) )
	suppressWarnings( downstreamFlank <- flank( gn, downstream, start = TRUE ) )

	#4-construct bins
	chrs <- seqnames( gn )
	starts <- start( gn )
	ends <- end( gn )
	strands <- strand( gn )
	minus <- which( strands == '-' )
	n <- length( starts )
	
	if(verbose)
		message( 'Creating bins...' )

	posList <- lapply( 1 : n, function( i ) round( seq( starts[ i ], ends[ i ], length = nBins + 1 ) ) )
	posList <- lapply( posList, function( x ) cbind( x[ -length( x ) ], x[ -1 ] ) )
	x <- lapply(posList, nrow )
	posMatrix <- do.call( rbind, posList )

	bins <- GRanges( seqnames = rep( chrs, each = nBins ), 
			 ranges    = IRanges( posMatrix[, 1], posMatrix[, 2] ), 
			 strand    = rep( strands, each = nBins ) )

	bins <- split( bins, strand( bins ) )
	
	#5-extract mean coverage or relative log odds from clusters
	if( odds ) {
		meanFeature <- elementMetadata( clusters )[, 'RelLogOdds' ]		
	} else {	
		meanFeature <- elementMetadata( clusters )[, 'MeanCov' ]
	}	

	#6-compute overlaps (strand-specific) for genes annotated on the +
	olaps <- findOverlaps( clusters,  bins$'+' )
	sameBin <- subjectHits( olaps ) %% nBins
	meanCovPlus <- split( meanFeature[  queryHits( olaps ) ], sameBin ) # $0 corresponds to last bin

	#6b-for genes annotated on the -
	olaps <- findOverlaps( clusters,  bins$'-' )
	sameBin <- nBins - subjectHits( olaps ) %% nBins + 1 #nBins+1 corresponds to first bin
	meanCovMinus <- split( meanFeature[  queryHits( olaps ) ], sameBin )
	names( meanCovMinus )[ which( names( meanCovMinus ) == nBins + 1 ) ] <- 0
	meanCovStar <- mapply( c, meanCovPlus, meanCovMinus, SIMPLIFY = FALSE )

	if( plot ) { 
		boxplot( meanCovStar, 
		      outline = FALSE,
		      notch = TRUE, 
 		      xaxt = 'n',
		      xlab = 'Position',
                      ylab = ifelse( odds, 'Relative log-odds', 'Mean coverage' ), 
		      lwd  = 1.25, 
		      col  = 'skyblue3',
		      ... )
		axis( 1, at = c( 1, 0.5 * nBins, nBins ), labels = c( 'TSS', '50%', 'TES' ) )
	}
}

