#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#2013 - Federico Comoglio & Cem Sievers, D-BSSE, ETH Zurich
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

getMetaTSS <- function( sortedBam, txDB = NULL, upstream = 1e3, downstream = 1e3, nBins = 40, genome = 'hg19', tablename = 'ensGene', unique = FALSE, plot = TRUE, verbose = TRUE, ... ) {
# return (and plot) metaTSS profile
#
# Args:
#   sortedBam: GRanges object, containing all raw data for wavClusteR2 analysis as returned by readSortedBam
#   txDB: TxDb object obtained through the makeTranscriptDbFromUCSC function. Default is NULL, namely the object will be fetched internally
#   upstream: numeric, the number of bases to be considered upstream the TSS. Default is 1000
#   downstream: numeric, the number of bases to be considered downstream the TSS. Default is 1000
#   nBins: numeric, the number of bins to be used. Defauls is 40, i.e. bin size 50 bases
#   genome: character, genome abbreviation used by UCSC and obtained by ucscGenomes()[ , "db"]. Default is human genome (hg19)
#   tablename: character, name of the UCSC table containing the transcript annotations to retrieve, according to supportedUCSCtables(). Default is ensembl gene annotation (ensGene)
#   unique: logical, if TRUE only non-overlapping TSS extended by upstream/downstream are considered. Default is FALSE 
#   plot: logical, if TRUE plots line plot with additional features. Default is TRUE
#   verbose: logical, if TRUE, prints steps. Default is TRUE
#   ...: additional parameters to be passed to the plot function
#
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

	#1-if not provided, create the TxDb object
	if( is.null( txDB ) ) {
		if(verbose)
			message( 'Creating TxDb object...' )
		txDB <- makeTranscriptDbFromUCSC(genome = genome, tablename = tablename)
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
		olaps <- findOverlaps( tss, ignoreSelf = TRUE )
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


getMetaGene <- function( clusters, txDB = NULL, upstream = 1e3, downstream = 1e3, nBins = 40, nBinsUD = 10, minLength = 1, genome = 'hg19', tablename = 'ensGene', plot = TRUE, verbose = TRUE, ... ) {
# return (and plot) metagene profile
#
# Args:
#   clusters: GRanges object as returned by the getClusters or the filterClusters function
#   txDB: TxDb object obtained through the makeTranscriptDbFromUCSC function. Default is NULL, namely the object will be fetched internally
#   upstream: numeric, the number of bases to be considered upstream the gene. Default is 1000
#   downstream: numeric, the number of bases to be considered downstream the gene. Default is 1000
#   nBins: numeric, the number of bins to be used to partition the genes. Default is 40.
#   nBinsUD: numeric, the number of bins to be used to partion upstream and downstream regions. Defauls is 10, i.e. bin size 100 bases for the default extension
#   minLength: numeric, the minimum length of a gene for being retained in the analysis. Defauls is 1, i.e. all genes are considered
#   genome: character, genome abbreviation used by UCSC and obtained by ucscGenomes()[ , "db"]. Default is human genome (hg19)
#   tablename: character, name of the UCSC table containing the transcript annotations to retrieve, according to supportedUCSCtables(). Default is ensembl gene annotation (ensGene)
#   plot: logical, if TRUE plots line plot with additional features. Default is TRUE
#   verbose: logical, if TRUE, prints steps. Default is TRUE
#   ...: additional parameters to be passed to the plot function
#
# Returns:
#   a numeric vector of the same length as nBins + 2 * nBinsUD with normalized counts. If plot, then also plots the distribution
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

	#1-if not provided, create the TxDb object
	if( is.null( txDB ) ) {
		if(verbose)
			message( 'Creating TxDb object...' )
		txDB <- makeTranscriptDbFromUCSC(genome = genome, tablename = tablename)
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

getMetaCoverage <- function( clusters, txDB = NULL, upstream = 1e3, downstream = 1e3, nBins = 40, nBinsUD = 10, minLength = 1, genome = 'hg19', tablename = 'ensGene', odds = FALSE, plot = TRUE, verbose = TRUE, ... ) {
# plot distribution of cluster mean coverage as a function of a metagene. It returns a boxplot per bin
#
# Args:
#   clusters: GRanges object as returned by the getClusters or the filterClusters function
#   txDB: TxDb object obtained through the makeTranscriptDbFromUCSC function. Default is NULL, namely the object will be fetched internally
#   upstream: numeric, the number of bases to be considered upstream the gene. Default is 1000
#   downstream: numeric, the number of bases to be considered downstream the gene. Default is 1000
#   nBins: numeric, the number of bins to be used to partition the genes. Default is 40.
#   nBinsUD: numeric, the number of bins to be used to partion upstream and downstream regions. Defauls is 10, i.e. bin size 100 bases for the default extension
#   minLength: numeric, the minimum length of a gene for being retained in the analysis. Defauls is 1, i.e. all genes are considered
#   genome: character, genome abbreviation used by UCSC and obtained by ucscGenomes()[ , "db"]. Default is human genome (hg19)
#   tablename: character, name of the UCSC table containing the transcript annotations to retrieve, according to supportedUCSCtables(). Default is ensembl gene annotation (ensGene)
#   odds: logical, if TRUE log odds and not average cluster is considered
#   plot: logical, if TRUE plots line plot with additional features. Default is TRUE
#   verbose: logical, if TRUE, prints steps. Default is TRUE
#   ...: additional parameters to be passed to the plot function
#
# Returns:
#   a numeric vector of the same length as nBins + 2 * nBinsUD with normalized counts. If plot, then also plots the distribution
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

	#1-if not provided, create the TxDb object
	if( is.null( txDB ) ) {
		if(verbose)
			message( 'Creating TxDb object...' )
		txDB <- makeTranscriptDbFromUCSC(genome = genome, tablename = tablename)
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

