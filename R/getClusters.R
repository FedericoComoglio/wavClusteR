#2013 - Federico Comoglio & Cem Sievers, D-BSSE, ETH Zurich


#' Identify clusters containing high-confidence substitutions and resolve
#' boundaries at high resolution
#' 
#' Identifies clusters using either the mini-rank norm (MRN) algorithm (default
#' and recommended to achieve highest sensitivity) or via a continuous wavelet
#' transform (CWT) based approach. The former employs thresholding of
#' background coverage differences and finds the optimal cluster boundaries by
#' exhaustively evaluating all putative clusters using a rank-based approach.
#' This method has higher sensitivity and an approximately 10-fold faster
#' running time than the CWT-based cluster identification algorithm. The
#' latter, maintained for compatibility with \code{wavClusteR}, computes the
#' CWT on a 1 kb window of the coverage function centered at a high-confidence
#' substitution site, and identifies cluster boundaries by extending away from
#' peak positions.
#' 
#' 
#' @usage getClusters(highConfSub, coverage, sortedBam, method = 'mrn', cores =
#' 1, threshold, step = 1, snr = 3)
#' @param highConfSub GRanges object containing high-confidence substitution
#' sites as returned by the \link{getHighConfSub} function
#' @param coverage An Rle object containing the coverage at each genomic
#' position as returned by a call to \link{coverage}
#' @param sortedBam a GRanges object containing all aligned reads, including
#' read sequence (qseq) and MD tag (MD), as returned by the
#' \link{readSortedBam} function
#' @param method a character, either set to "mrn" or to "cwt" to compute
#' clusters using the mini-rank norm or the wavelet transform-based algorithm,
#' respectively. Default is "mrn" (recommended).
#' @param cores integer, the number of cores to be used for parallel
#' evaluation. Default is 1.
#' @param threshold numeric, if \code{method = "mrn"}, the difference in
#' coverage to be considered noise. If not specified, a Gaussian mixture model
#' is used to learn a threshold from the data. Empirically, 10\% of the minimum
#' coverage required at substitutions (see argument \code{minCov} in the
#' \link{getHighConfSub} function) might suffice to provide highly resolved
#' clusters. However, if \code{minCov} is much lower than the median
#' strand-specific coverage at substitutions \eqn{m}, which can be computed
#' using \code{summary(elementMetadata(highConfSub)[, 'coverage'])['Median']}),
#' 10\% of \eqn{m} might represent an optimal choice.
#' @param step numeric, if \code{method = "cwt"}, step size of window shift. If
#' two high-confidence substitution sites are located within a distance less
#' than \code{step}, the wavelet transform is computed only once. Default: 1,
#' i.e. each high-confidence substitution site is considered independently.
#' @param snr numeric, if \code{method = "cwt"}, signal-to-noise ratio
#' controlling the peak calling as performed by \code{wavCWTPeaks} implemented
#' in the wmtsa package. Default: 3.
#' @return GRanges object containing the identified cluster boundaries.
#' @note Clusters returned by this function need to be further merged by the
#' function \code{filterClusters}, which also computes all relevant cluster
#' statistics.
#' @author Federico Comoglio and Cem Sievers
#' @seealso \code{\link{getHighConfSub}}, \code{\link{filterClusters}}
#' @references William Constantine and Donald Percival (2011), wmtsa: Wavelet
#' Methods for Time Series Analysis,
#' \url{http://CRAN.R-project.org/package=wmtsa}
#' @keywords core
#' @examples
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
#' @export getClusters
getClusters <- function(highConfSub, coverage, sortedBam, method = 'mrn', cores = 1, threshold, step = 1, snr = 3) {
# Error handling
#   if method is not within 'coverage' or 'cwt', raise an error

	stopifnot( method %in% c('mrn', 'cwt') )

	if(method == 'mrn') { 
		clusters <- getClustersMRN( highConfSub = highConfSub, 
					    coverage     = coverage, 
					    sortedBam    = sortedBam,
				            threshold    = threshold,
					    cores        = cores )
		metadata( ranges( clusters ) ) <- list( 'mrn' ) #fills slot with method info
	}
	else {
		clusters <- getClustersCWT( highConfSub = highConfSub, 
					    coverage    = coverage, 
					    step        = step, 
					    snr         = snr,
					    cores       = cores )
		metadata( ranges( clusters ) ) <- list( 'cwt' ) #fills slot with method info
	}
	clusters <- sort( clusters ) #sorting facilitates strand attribution

	return( clusters )
}


getWindow <- function( posTransition, nonZeroPos, diff, jumps ) {
# Given the genomic position of a transition, identify the genomic positions of continuous non-zero coverage that contain the transition
#
# Args:
#   pos.transition: numeric, the genomic position of the transition being considered
#   nonzero.pos: numeric vector, the positions in the coverage function that exhibit a non-zero coverage
#   diff: numeric vector, the difference between nonzero.pos
#   jumps: numeric vector, which of the above is not a consecutive index
#
# Returns:
#   A numeric vector with the sequence of genomic positions within the window
#
# Error handling
#   ...
	#1-localize the transition within positions of non-zero coverage (guaranteed to be non-empty)
	hit <- which( nonZeroPos == posTransition )
	
	#2-localize the left and right boundaries of the window
	leftIdx <- tail( jumps[ jumps < hit ], 1 ) + 1
	rightIdx <- jumps[ jumps > hit ][1]

	if( length( leftIdx ) == 0 || is.na( leftIdx ) | length(rightIdx) == 0 || is.na( rightIdx ) ) { #rare instances where from chr start to first transition coverage always > 0 OR from last transition to chr end coverage always > 0
		leftIdx <- hit
		rightIdx <- hit
	}

	#3-transform relative to absolute indices
	left <- nonZeroPos[ leftIdx ]
	right <- nonZeroPos[ rightIdx ]

	#4-construct the window
	window <- left : right

	return(window)
}

getSEcoverage <- function( sortedBam ) {
# Computes the start and the end coverage functions, obtained from considering read start and end only, respectively
#
# Args:
#   sortedBam: a GRanges object as read in by the readSortedBam function
#
# Returns:
#   a list with two slots. Each slot is an Rle object containing the start and end coverage, respectively.
#
# Error handling
#   if no sortedBam is provided, raise an error
	if( missing( sortedBam ) )
		stop('No sorted BAM file provided. Please read BAM files using the readSortedBam function and provide the returned GRanges object as input.')	
	
	#1-resize GRanges object to contain start position of reads only and compute coverage
	strand( sortedBam ) <- '*' #in order to use the resize function, s.t. start of reads on the - strand is not the read end
	sortedBamResized <- resize( sortedBam, width = 1, fix = 'start' )
	coverStart <- coverage( sortedBamResized )

	#1-resize GRanges object to contain end position of reads only and compute coverage
	sortedBamResized <- resize( sortedBam, width = 1, fix = 'end' )
	coverEnd <- coverage( sortedBamResized )

	return( list( coverStart = coverStart, coverEnd = coverEnd ) )	
}

learnThreshold <- function( highConfSubset, coverage ) {
# Fit a two-components Gaussian Mixture Model to learn a global threshold used to estimate background/noise jumps in the coverage function during the MRN cluster boundaries identification
#
# Args:
#   highConfSubset: a GRanges object containing a subset of highConfSub to be used to estimate the threshold
#   coverage: Rle object containing sequencing coverage
#
# Returns:
#   a numeric value of the threshold, in [0,1]
#
# Error handling
#   ...

	subset <- resize( highConfSubset, 50, fix = 'center' )
	subset <- reduce( subset )
		subset <- as.data.frame( subset )[, 2 : 3]
		chr <- as.character( seqnames( highConfSubset[ 1 ] ) )
		diffChr <- diff( coverage[[ chr ]] )

		X <- apply(subset, 1, function( interval ) {
			xd <- diffChr[ interval[ 1 ] : interval[ 2 ] ]
			xd <- xd[ xd > 0 ]
			if( length( xd ) > 0 ) {
				max <- max( xd )
				as.numeric( xd / max )
			}
		})

		X <- unlist( X )
		fit <- Mclust( X, G = 2, modelNames = 'V' )
		threshold <- min( X[ fit$classification == 2 ] )

	return( threshold )	
}

getClustersMRN <- function( highConfSub, coverage, sortedBam, threshold, cores ) {
# Identify clusters using the 'coverage' method
#
# Args:
#   highConfSub: GRanges object containing high-confidence transitions
#   coverage: Rle object containing sequencing coverage
#   sortedBam: list, contains BAM file slots as read in by the readSortedBam function
#   threshold: numeric, the difference in coverage to be considered noise. Default is no value specified, hence learning it from data
#   cores: numeric, the number of cores to be used for parallel evaluation
#
# Returns:
#   A GRanges object containing clusters
#
# Error handling
#   ...

	message( 'Computing start/end read positions' )
	SEcoverage <- getSEcoverage( sortedBam )

   	highConfSubSplit <- split( highConfSub, seqnames( highConfSub ) )
	
	#0-if the background threshold should be learnt, consider transitions on the first chromosome for this purpose
	treshold.provided <- TRUE
	if( missing( threshold ) ) {
		treshold.provided <- FALSE
		message( 'Learning background threshold by fitting a GMM' )
		subset <- highConfSubSplit[[ 1 ]]
		sampleSize <- length( subset )
		sampleSize <- ifelse( sampleSize > 1e3, 1e3, sampleSize ) 
		subset <- subset[ 1 : sampleSize ]
		threshold <- learnThreshold( subset, coverage )
		message( '   Estimated threshold (% of maximum local coverage differences) from ', sampleSize, ' sampled transitions: ', signif( threshold * 100, 3 ) )			
	}

	chrom <- as.character( unique( seqnames( highConfSub ) ) )
	message( 'Number of chromosomes exhibiting high confidence transitions: ', length( chrom ) )

 	#backend for doMC (windows)
	if ( suppressWarnings( require( "doMC" ) ) ) {
		registerDoMC( cores = cores )
     	} else {
     	}

	lev = NULL; rm( lev ) #pass Rcheck
	clusters <- foreach( lev = chrom, .combine = rbind ) %dopar% {
      		message( '...Processing = ', lev )
      		tempGR <- highConfSubSplit[[ lev ]]
      		chrPileup <- coverage[[ lev ]]
		startPileup <- SEcoverage$coverStart[[ lev ]] 
		endPileup <- SEcoverage$coverEnd[[ lev ]] 
		nonZeroPos <- which( chrPileup != 0 )
		nonZeroPosDiff <- diff( nonZeroPos )
		jumps <- which( nonZeroPosDiff != 1 )
		n <- length( tempGR )	

		#0-preallocate output
		startVector <- rep( Inf, n )
		endVector <- startVector	
		starts <- start( tempGR )
		strands <- strand( tempGR )		

		for( i in seq_len( n ) ) {
			posTransition <- starts[ i ]	#w.r.t. genomic coordinates
			window <- getWindow( posTransition, nonZeroPos, nonZeroPosDiff, jumps )

			if( length( window ) > 1 ) { #if = 1, then hcTC observed at start/end chromosome (p << 1e-4)
				pw <- which(window == posTransition)	#position of the transition in the window
				
				#1-get start and end coverage in window
				ws <- as.numeric( startPileup[ window ] )
				we <- as.numeric( endPileup[ window ] )
				
				#2-filter positions that exhibit differences higher than background threshold (yields candidate boundaries)
				if( !treshold.provided ) {
					deltaS <- floor( threshold * max( ws[ 1 : pw ] ) )
					deltaE <- floor( threshold * max( we[ pw : length( we ) ] ) )

					putativeStartIdx <- which( ws >= deltaS ) 
					putativeEndIdx <- which( we >= deltaE )	

					putativeStartIdx <- putativeStartIdx[ putativeStartIdx <= pw ]
					putativeEndIdx <- putativeEndIdx[ putativeEndIdx >= pw ]
		
				}
				else { #threshold provided
					putativeStartIdx <- which( ws >= threshold ) 
					putativeEndIdx <- which( we >= threshold )				
				
					#1b-if no suitable position found, lower threshold to 1 (implies always feasible solution)	
					if( length( putativeStartIdx[ putativeStartIdx <= pw ] ) == 0 | length( putativeEndIdx[ putativeEndIdx >= pw ] ) == 0 ) {
						putativeStartIdx <- which( ws >= 1 )
						putativeEndIdx <- which( we >= 1 )
					}
					putativeStartIdx <- putativeStartIdx[ putativeStartIdx <= pw ]
					putativeEndIdx <- putativeEndIdx[ putativeEndIdx >= pw ]
					}

				#3-consider the SEcoverage values at the positions
				sCov <- ws[ putativeStartIdx ]
				eCov <- we[ putativeEndIdx ]

				#4-rank coverage values (directional rank)
				rankS <- rank( -sCov, ties.method = 'first' )
				rankE <- rev( rank( rev( -eCov ), ties.method = 'first' ) ) #reversing = apply "first" 3' --> 5'

				#5-consider the cartesian product of the ranks
				cart <- expand.grid( rankS, rankE )
				cartSE <- expand.grid( putativeStartIdx, putativeEndIdx )

				#6-compute theoretical cluster size for each combination	
				clustSize <- apply( cartSE, 1, function(x) x[2] - x[1] + 1 )
				
				#7-rank cluster sizes and augment cart
				rankCs <- rank( clustSize, ties.method = 'min' )
				cart <- cbind( cart, rankCs, cartSE, clustSize )
	
				#8-consider the subset (ideally singleton) of cart that minimizes l2 norm
				l2 <- apply( cart[, 1 : 3] - 1, 1, function(x) sqrt( sum( x ^ 2 ) ) )  #-1, obj opt: (1,1,1)
				candidates <- which( l2 == min( l2 ) )

				#9-if ties, then select the larger window
				selected <- as.numeric( cart[candidates, ][ which.max( cart[ candidates, 6 ] )[ 1 ], 4 : 5 ] )

				#10-add to the list of cluster boundaries
				startVector[ i ] <- window[  selected[1] ]
				endVector[ i ] <- window[ selected[2] ] 

			}
			else {
				startVector[ i ] <- window
				endVector[ i ] <- window	
			}			
		}
	
	cbind( rep( lev, length( startVector ) ), startVector, endVector, as.character( strands ) )
	}
	
	clusters <- GRanges( seqnames = as.vector( clusters[, 1] ), 
			     ranges   = IRanges( as.numeric( as.vector( clusters[, 2] ) ), 
					         as.numeric( as.vector( clusters[, 3] ) ) ),
			     strand   = clusters[, 4] )
	return( clusters )
}


getClustersCWT <- function(highConfSub, coverage, step, snr, cores) {
#maintained pipeline (reproducibility of wavClusteR results), see Sievers et al. (2012), Nucl Acids Res 40(2):e160
	#backend for doMC (windows)
	if ( suppressWarnings( require( "doMC" ) ) ) {
		registerDoMC( cores = cores )
     	} else {
     	}

	highConfSubSplit <- split( highConfSub, seqnames( highConfSub ) )
	chrom <- seqlevels(highConfSub)
	message( 'List of chromosomes in highConfSub: \n', chrom, '\n')
	chromCov <- names( coverage )
	message('List of chromosomes in coverage: \n', chromCov, '\n')
    	message('Chromosomes being processed:', chrom, '\n')
#    	clusters <- data.frame()
	lev <- NULL; rm( lev ) #pass Rcheck
	clusters <- foreach( lev = chrom, .combine = rbind ) %dopar% {
#      	for(lev in chrom) {      
      	message('chromosome being considered = ', lev, '\n')
      	index <- which(chrom == lev)
      	temp.gr <- highConfSubSplit[[index]]
      	index <- which(chromCov == lev)
      	chr.pileup <- coverage[[index]]
      	max.pos <- length( chr.pileup )  
      	pos.diff <- c(1, diff( start(temp.gr) ))
      	start.vector <- c()
      	end.vector <- c()  
      	if(length(temp.gr) > 0) {
		for(i in seq_len(length(temp.gr))) {  #error was at chr18 (n = 10), position = 8.
		tryCatch({
         		if( pos.diff[i] >= step ) {
        			window <- c( max( 1, (start(temp.gr[i]) - 500) ), min( (start(temp.gr[i]) + 500), max.pos) ) #1kb
         			signal <- getSignal( chr.pileup, window ) #chr.pileup = coverage of the chromosome being considered
				len.support <- length(signal)
				if(len.support >= 1000) {
             				res.cwt <- wavCWT( signal )
             				res.tree <- wavCWTTree( res.cwt, n.octave.min = 1 )
            				res.peaks <- wavCWTPeaksCF(res.tree, snr.min = snr, noise.span = 0)
           				if(!is.null(res.tree)) { #res.tree
             					delta <- window[1] #start(temp.gr[i])-500
            					for(peak in res.peaks$x) {
         						start.vector <- c(start.vector, leftBound(peak, signal) + delta - 1)
             						end.vector <- c(end.vector, rightBound(peak, signal) + delta - 1)
          					}
         				}
         			}
   			}
      		}, error = function(e) print(signal))
		}
	}
	cbind(rep(lev, length(start.vector)), start.vector, end.vector)
	}
	clusters <- GRanges(as.vector(clusters[,1]), IRanges(as.numeric(as.vector(clusters[,2])), as.numeric(as.vector(clusters[,3]))))
return(clusters)
}

