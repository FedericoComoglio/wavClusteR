#2013 - Federico Comoglio & Cem Sievers, D-BSSE, ETH Zurich

getLogOdd <- function( model, k, n ) {
# Compute log odds for a given number of reads k exhibiting a transition over a total coverage n at that genomic position
#
# Args:
#   model: list containing the parameters of the mixture model
#   k: numeric, the number of reads exhibiting a transition at a given genomic position
#   n: numeric, the strand-specific coverage at that position
#
# Returns:
#   ...
#
# Error handling
#   ...

	p1 <- model$p1
	p2 <- model$p2
	l1 <- model$l1
	l2 <- model$l2
	nVals <- length( p1 )
	dx <- 1 / nVals

	xVals <- seq( 0, 1, length = nVals )
   
	#1-compute the joint P(Theta, Y) from P(Theta, Y, X)
	#  p(Y|X) ~ Bin(n,k)
	#  p(x|Theta) is given by mixture model components
	#  p(Theta) is given by mixing coefficients
	bin <- dbinom( k, n, xVals )
	joint1 <- sum( bin * p1 ) * l1 * dx 
  	joint2 <- sum( bin * p2 ) * l2 * dx 
  
	#2-compute P(Y) by marginalizing P(Theta, Y) over Theta
  	py <- joint1 + joint2
  	
	#3-compute posteriors p(Theta|Y)
	post1 <- joint1 / py
  	post2 <- joint2 / py
	
	#4-compute log odds
	logOdds <- log( post2 / post1 )
	
	return( logOdds )
}


computeLogOdds <- function( highConfSub, model ) {
# Compute log odds for a given number of reads k exhibiting a transition over a total coverage n at that genomic position
#
# Args:
#   highConfSub: GRanges object containing high-confidence transitions
#   model: list containing the parameters of the mixture model
#
# Returns:
#   a GRanges object, identical to highConfSub with an additional metadata containing the computed logOdds for every transition
#
# Error handling
#   ...

	#1-extract metadata
	which <- c( 'coverage', 'count' )
	emd <- as.data.frame( elementMetadata( highConfSub )[, which] )

	#2-add a column to be used to sort the results
	emd[, 'inOrder'] <- 1 : nrow( emd )

	#3-create a hash table of unique pairs (coverage, counts) and use it to compute log odds
	hash <- as.data.frame( emd[ !duplicated( emd [, which] ), which ] )
	hashOdds <- apply( hash, 1, function(x) getLogOdd( model, x[ 2 ], x[ 1 ] ) )
	hash <- cbind( hash, hashOdds )

	#4-merge results to create an object of the same length as emd
	logOdds <- merge( emd, hash )

	#5-sort the result and append to input
	logOdds <- logOdds[ order( logOdds[, 'inOrder'] ),  ]	
	elementMetadata( highConfSub )[ 'logOdds' ] <- logOdds[, 'hashOdds']		

	return( highConfSub )
}


getOdd <- function(xvec, yvec, x) {
# Identifies the values of yvec which is closest to mapped values of xvec
#
# Args:
#  xvec: a numeric vector
#  yvec: a numeric vector
#  x: numeric
#
# Returns:
#  the value of yvec at the xvec value closest to x
#
	d <- abs( xvec - x )
   	pos <- which.min( d )
    	return( yvec[ pos ] )
}

computelogOdds <- function( model ) {
# Compute log odds using components of mixture model
#
# Args:
#   model: list containing the parameters of the mixture model
#
# Returns:
#   a numeric vector, of the same length as model$p1 (because computed on the same xval interval)
#
# Error handling
#   ...
	logOdds <- log( ( model$l2 * model$p2 ) / ( model$l1 * model$p1 ) )
	return( logOdds )
}




#' Merge clusters and compute all relevant cluster statistics
#' 
#' If clusters have been identified using the mini-rank norm algorithm, cluster
#' statistics are computed. In contrast, if the CWT-based cluster
#' identification algorithm was used, clusters are first filtered to retain
#' only those instances containing a wavelet peak and a high-confidence
#' substitution site within their cluster boundaries.
#' 
#' 
#' @usage filterClusters(clusters, highConfSub, coverage, model, genome,
#' refBase = 'T', minWidth = 12, verbose = TRUE)
#' @param clusters GRanges object containing individual clusters as identified
#' by the \link{getClusters} function
#' @param highConfSub GRanges object containing high-confidence substitution
#' sites as returned by the \link{getHighConfSub} function
#' @param coverage An Rle object containing the coverage at each genomic
#' position as returned by a call to \link{coverage}
#' @param model List of 5 items containing the estimated mixture model as
#' returned by the \link{fitMixtureModel} function
#' @param genome BSgenome object of the relevant reference genome (e.g.
#' \code{Hsapiens} for the human genome hg19)
#' @param refBase A character specifying the base in the reference genome for
#' which transitions are experimentally induced (e.g. 4-SU treatment - a
#' standard in PAR-CLIP experiments - induces T to C transitions and hence
#' \code{refBase = "T"} in this case). Default is "T"
#' @param minWidth An integer corresponding to the minimum width of reported
#' clusters. Shorter clusters are extended to \code{minWidth} starting from the
#' cluster center
#' @param verbose Logical, if TRUE processing steps are printed
#' @return GRanges object containing the transcriptome-wide identified
#' clusters, having metadata: \item{Ntransitions}{The number of high-confidence
#' transitions within the cluster} \item{MeanCov}{The mean coverage within the
#' cluster} \item{NbasesInRef}{The number of genomic positions within the
#' cluster corresponding to \code{refBase}} \item{CrossLinkEff}{The
#' crosslinking efficiency within the cluster, estimated as the ratio between
#' the number of high-confidence transitions within the cluster and the total
#' number of genomic positions therein corresponding to \code{refBase}}
#' \item{Sequence}{The genomic sequence undelying the cluster (plus strand)}
#' \item{SumLogOdds}{The sum of the log-odd values within the cluster}
#' \item{RelLogOdds}{The sum of the log-odds divided by the number of
#' high-confidence transitions within the cluster. This variable can be
#' regarded as a proxy for statistical significance and can be therefore used
#' to rank clusters. See Comoglio, Sievers and Paro for details.}
#' @note 1) This function calls the appropriate processing function according
#' to the method used to compute clusters. This information is stored in the
#' \code{metadata(ranges(clusters))} slot as an object of type list.
#' 
#' 2) Notice that \code{genome} corresponds to the according reference genome
#' matching the organism in which experiments have been carried out. For
#' example \code{genome = Hsapiens} is used for the human reference genome
#' (assembly 19), where \code{Hsapiens} is provided by
#' \code{BSgenome.Hsapiens.UCSC.hg19}.
#' @author Federico Comoglio and Cem Sievers
#' @seealso \code{\link{getClusters}}, \code{\link{getHighConfSub}},
#' \code{\link{fitMixtureModel}}
#' @references Herve Pages, BSgenome: Infrastructure for Biostrings-based
#' genome data packages
#' 
#' Sievers C, Schlumpf T, Sawarkar R, Comoglio F and Paro R. (2012) Mixture
#' models and wavelet transforms reveal high confidence RNA-protein interaction
#' sites in MOV10 PAR-CLIP data, Nucleic Acids Res. 40(20):e160. doi:
#' 10.1093/nar/gks697
#' 
#' Comoglio F, Sievers C and Paro R.
#' 
#' @keywords core
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
#' 	                 method = 'mrn', 
#' 	                 cores = 1, 
#' 	                 threshold = 2 ) 
#' 
#' fclusters <- filterClusters( clusters = clusters, 
#' 		             highConfSub = highConfSub, 
#'         		     coverage = coverage,
#' 			     model = model, 
#' 			     genome = Hsapiens, 
#' 		             refBase = 'T', 
#' 		             minWidth = 12 )
#' fclusters
#' 
#' @export filterClusters
filterClusters <- function( clusters, highConfSub, coverage, model, genome, refBase = 'T', minWidth = 12, verbose = TRUE ) {
# Error handling
#   if method is not within 'coverage' or 'cwt', raise an error
	method <- metadata( ranges( clusters ) )[[1]] #slot is a list
	stopifnot( method %in% c( 'mrn', 'cwt' ) )

	if( method == 'mrn' ) 
		clusters <- filterClustersMRN( clusters     = clusters, 
					      highConfSub  = highConfSub, 
					      coverage     = coverage, 
					      model        = model, 
					      genome       = genome, 
					      refBase      = refBase, 
					      minWidth     = minWidth, 
					      verbose      = verbose )
	else
		clusters <- filterClustersCWT( clusters     = clusters, 
					      highConfSub  = highConfSub, 
					      coverage     = coverage, 
					      model        = model, 
					      genome       = genome, 
					      refBase      = refBase, 
					      minWidth     = minWidth )

	return( clusters )
}


filterClustersMRN <- function( clusters, highConfSub, coverage, model, genome, refBase = 'T', minWidth = 12, verbose = TRUE ) {
# Identifies clusters with high-confidence transitions and returns all relevant descriptors and statistics
#
# Args:
#   clusters: GRanges object containing individual clusters as identified by getClusters
#   highConfSub: GRanges object containing high-confidence transitions
#   coverage: Rle object containing sequencing coverage
#   model: list containing the parameters of the mixture model
#   genome: BSgenome object of the relevant genome (e.g. Hsapiens)
#   refBase: character, the base in the reference genome for which transitions are induced. Default is T (for T-->C)
#   minWidth: numeric, minimum cluster width. Clusters narrower than minWidth will be extended to that size. Default is 12
#   verbose: logical, if TRUE, prints steps. Default is true
#
# Returns:
#   A GRanges object containing valid clusters for post-processing (e.g. annotation) or export
#
# Error handling
#   ...
	
	#1-identify base complementary to reference base
	bases <- c('A', 'T', 'C', 'G')
	basesCpl <- c('T', 'A', 'G', 'C')
	refBaseCpl <- basesCpl[ which( bases == refBase ) ]
 
	#2-compute log odds
	if( verbose )
		message( 'Computing log odds...' )
	highConfSub <- computeLogOdds( highConfSub, model )

	#3-check cluster sizes and increase size to 12 if w < 12
	if(verbose)
		message( 'Refining cluster sizes...' )
	w <- width( clusters )
	whichFixed <- which( w < minWidth )
	fixedClusters <- resize( clusters[whichFixed], minWidth, fix = 'center' )
	clusters[ whichFixed ] <- fixedClusters

	#3b-check that resized clusters do not exceed coverage boundaries
	chr <- unique( seqnames( clusters ) )	
	sn <- seqnames( clusters )
	for( lev in chr ) {
		endCov <- length( coverage[[ lev ]] )
		outside <- end( clusters[ sn == lev ] ) > endCov
		end( clusters[ sn == lev ] )[ outside ] <- endCov 
	}

	#4-assign strand to clusters. If done before reduce, clusters on different strands will not be merged
#	if( verbose )
#		message( 'Assigning strand information to clusters...' )
#	highConfSubNoStrand <- highConfSub
#	strand( highConfSubNoStrand ) <- '*'
#	highConfSubNoStrand <- sort( highConfSubNoStrand ) #generates 1:1 map substitution:cluster (order preserved)
#	olaps <- findOverlaps( highConfSubNoStrand, highConfSub )	
#	strand( clusters ) <- strand( highConfSub[ subjectHits( olaps ) ] )
	
	#5-reduce clusters
	if( verbose )
		message( 'Combining clusters...' )
	clusters <- reduce( clusters )
	minusStrand <- which( strand( clusters ) == '-' )
	n <- length( clusters )
	
	#6-count the number of transitions per cluster
	if( verbose )
		message( 'Quantifying transitions within clusters...' )
	nTransitions <- countOverlaps( clusters, highConfSub )

	#7-compute statistics for each cluster
	if( verbose )
		message( 'Computing statistics...' )
	chrClusters <- as.character( seqnames( clusters ) )
	startClusters <- start( clusters )
	endClusters <- end( clusters )
	sequences <- as.vector( getSeq( genome, chrClusters, startClusters, endClusters, as.character = TRUE ) )	#to compute #bases of transition type
	sumOddsVec <- meanCovVec <- c()	#vectors containing future elementMetadata

	pb <- txtProgressBar( min = 0, max = n, initial = 1, char = "=", style = 3 )

	olaps <- findOverlaps( highConfSub, clusters )
	mapping <- split( queryHits( olaps ), subjectHits( olaps ) )

	logOdds <- elementMetadata( highConfSub )[, 'logOdds']
#prevLO	rsfs <- elementMetadata( highConfSub )[, 'rsf']
	
	for( i in seq_len( n ) ) {
#prevLO		rsf <- rsfs[ mapping[[ i ]] ]
		odds <- logOdds[ mapping[[ i ]] ]
		chr <- chrClusters[ i ]
#prevLO		sumOdds <- sum( sapply( rsf, getOdd, xvec, logOdds ) )	#get sum log odds	
		sumOdds <- sum( odds )
          	meanCov <- mean( coverage[[ chr ]][ startClusters[ i ] : endClusters[ i ] ] ) 
		sumOddsVec <- c( sumOddsVec, sumOdds )
		meanCovVec <- c( meanCovVec, meanCov )
		setTxtProgressBar( pb, i )
	}
	close( pb )

	toCount <- rep( refBase, n )
  	toCount[ minusStrand ] <- refBaseCpl
  	refBaseCountVec <- str_count( sequences, toCount )
    	relOddsVec <- sumOddsVec / refBaseCountVec
	crosslinkEff <- round( nTransitions / refBaseCountVec, 2 )

	#8-consolidate results and prepare output
	if( verbose )
		message( 'Consolidating results...' )

	emdDf <- data.frame( Ntransitions = nTransitions,
                             MeanCov      = meanCovVec,
                             NbasesInRef  = refBaseCountVec,
			     CrossLinkEff = crosslinkEff,
			     Sequence     = sequences,
                             SumLogOdds   = sumOddsVec,
			     RelLogOdds   = relOddsVec ) #data.frame with elementMetadata
	
	elementMetadata( clusters ) <- emdDf
	
	return( clusters )
}


filterClustersCWT <- function( clusters, highConfSub, coverage, model, genome, refBase = 'T', minWidth = 12) { 
#maintained pipeline (reproducibility of wavClusteR results), see Sievers et al. (2012), Nucl Acids Res 40(2):e160
  	ref.base.compl <- as.character( complement( DNAString(refBase) ) )
  	chr.cov <- names( coverage )
  	xvec <- seq( .001, .999, .001 )
  	logOdds <- computelogOdds( model )
  	w <- width( clusters )

	clust.small <- clusters[w < minWidth]
	clust.large <- clusters[w >= minWidth]
    
	if(length(clust.small) > 0) 
		clust.small <- resize(clust.small, minWidth, fix = 'center')
    	
	clusters <- c(clust.small, clust.large)
    	clusters <- reduce(clusters)
    	over <- findOverlaps(highConfSub, clusters)
    	over.mm <- as.matrix( over )
    	rel.freq <- elementMetadata(highConfSub)[, 'rsf']     
    	chr.v <- start.v <- end.v <- strand.v <- sum.odds.v <- mean.cov.v <- c()
    	index <- c()
	for(i in seq_len(nrow(over.mm))) {
      		clust.ind <- over.mm[i, 2]
      		sub.ind <- over.mm[i, 1]
      		tmp.sub <- highConfSub[ sub.ind ]
     		tmp.clust <- clusters[ clust.ind ]
      		tmp.odds <- getOdd( xvec, logOdds, rel.freq[sub.ind] )   
      	
		if(clust.ind %in% index) {
      			n <- length(sum.odds.v)
      			sum.odds.v[ n ] <- sum.odds.v[ n ] + tmp.odds
      		}
      		else {
       			index <- c(index, clust.ind)
       			chr <- as.character( seqnames(tmp.clust) )
       			s <- start(tmp.clust)
       			e <- end(tmp.clust)
		        idx <- which(chr.cov == chr)
          		mean.cov <- mean( coverage[[ idx ]][s : e] ) 
          		chr.v <- c(chr.v, chr)
    			start.v <- c(start.v, s)
    			end.v <- c(end.v, e)
   			strand.v <- c( strand.v, as.character(strand(tmp.sub)) )
			sum.odds.v <- c(sum.odds.v, tmp.odds)
			mean.cov.v <- c(mean.cov.v, mean.cov)
		}
	}  

	sequences <- as.vector( getSeq(genome, chr.v, start.v, end.v, as.character = TRUE) )
  	to.count <- rep( refBase, length(chr.v) )
  	to.count[ strand.v == '-' ] <- ref.base.compl
  	ref.base.count.v <- str_count(sequences, to.count)
    	rel.odds.v <- sum.odds.v / ref.base.count.v
    	res.gr <- GRanges(chr.v, IRanges(start.v, end.v ), strand.v,
    	mean.coverage = mean.cov.v, base.count = ref.base.count.v, sum.odds = sum.odds.v, relative.odds = rel.odds.v)     
 
	return(res.gr)
}
