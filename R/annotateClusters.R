#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#2013 - Federico Comoglio & Cem Sievers, D-BSSE, ETH Zurich
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

annotateClusters <- function( clusters, txDB = NULL, genome = 'hg19', tablename = 'ensGene', plot = TRUE, verbose = TRUE ) {
# Annotates clusters w.r.t. genomic annotations from UCSC
#
# Args:
#   clusters: GRanges object containing individual clusters as identified by getClusters
#   txDB: TranscriptDb object obtained through the makeTranscriptDbFromUCSC function. Default is NULL, namely the object will be fetched internally
#   genome: character, genome abbreviation used by UCSC and obtained by ucscGenomes()[ , "db"]. Default is human genome (hg19)
#   tablename: character, name of the UCSC table containing the transcript annotations to retrieve, according to supportedUCSCtables(). Default is ensembl gene annotation (ensGene)
#   plot: logical, if TRUE plots dotchart
#   verbose: logical, if TRUE, prints steps. Default is true
#
# Returns:
#   GRanges object witwh clusters, with an additional elementMetadata column containing a character identifier of the genomic feature each cluster maps to. If plot, then a dotchart is returned in addition.
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
		txDB <- makeTranscriptDbFromUCSC(genome = genome, tablename = tablename)
	}
	
	#2-obtaining CDS, introns, 3' and 5'-UTRs from the TranscriptDb object, make them unique and compute their length
	if(verbose)
		message( 'Extracting genomic features from TranscriptDb object...' )

	cds <- cdsBy( txDB )
	ix <- intronsByTranscript( txDB )
	tpUTR <- threeUTRsByTranscript( txDB )
	fpUTR <- fiveUTRsByTranscript( txDB )

	cds <- unique( unlist( cds ) )
	ix <- unique( unlist( ix ) )
	tpUTR <- unique( unlist( tpUTR ) )
	fpUTR <- unique( unlist( fpUTR ) )
	
	nBasesCompartments <- c( nBasesCds = sum( as.numeric( width( cds ) ) ), 
				 nBasesIx = sum( as.numeric( width( ix ) ) ), 
				 nBasesTp = sum( as.numeric( width( tpUTR ) ) ), 
				 nBasesFp = sum( as.numeric( width( fpUTR ) ) ) )
	totBases <- sum( nBasesCompartments )	
	proportionGenome <- nBasesCompartments / totBases * 100

	#3-find overlaps on the sense strand, return a matrix
	if(verbose)
		message( 'Computing overlaps with genomic features on the sense strand...' )
	
	olaps <- cbind( olapsCds   = countOverlaps( clusters, cds ),
      		        olapsIx    = countOverlaps( clusters, ix ),
      		        olapsTpUTR = countOverlaps( clusters, tpUTR ),
     		        olapsFpUTR = countOverlaps( clusters, fpUTR) )
	
	olaps <- olaps > 0 #transform to logical
	rs <- rowSums( olaps )
	whichNotMapped <- which( rs == 0 ) #vector of cluster indices for which no mapping was found	

	#4-prepare output
	n <- length( clusters )  
	emdAnno <- rep( '', n ) #future elementMetadata containing annotation label for each cluster

	#5-flag ambigous/multiple mappings
	multiple <- which( rs > 1 ) #logical comparison, 1 if ambigous
	emdAnno[ multiple ] <- 'Multiple'
	olaps[ multiple, ] <- FALSE  #multiple mappings are cleared

	categoryIdx <- apply( olaps, 2, which, TRUE ) #list of indices of clusters mapping to same feature
	emdAnno[ categoryIdx$olapsCds ] <- 'CDS ss'
	emdAnno[ categoryIdx$olapsIx ] <- 'Introns ss'
	emdAnno[ categoryIdx$olapsTpUTR ] <- '3\' UTR ss'
	emdAnno[ categoryIdx$olapsFpUTR ] <- '5\' UTR ss'
		
	#5-count overlaps by feature
	summ <- colSums( olaps )
	nMultiple <- length( multiple )
	summ <- c( summ, nMultiple, n - sum( summ ) - nMultiple )	#n-sum(summ): not mapped
	
	#6-consider not-mapped ones and map them w.r.t antisense strand
	if( length( whichNotMapped ) > 0 ) {
		if(verbose)
			message( 'Considering non-mapped clusters and computing overlaps with genomic features on the antisense strand...' )

		notMapped <- clusters[ whichNotMapped ]
		plus <- strand( notMapped ) == '+'
		minus <- strand( notMapped ) == '-'
		strand( notMapped[ plus ] ) <- '-'	
		strand( notMapped[ minus ] ) <- '+'
	
		olapsAntisense <- cbind(  olapsCds   = countOverlaps( notMapped, cds ),
      		       			  olapsIx    = countOverlaps( notMapped, ix ),
      		        		  olapsTpUTR = countOverlaps( notMapped, tpUTR ),
     		        		  olapsFpUTR = countOverlaps( notMapped, fpUTR) )
	
		olapsAntisense <- olapsAntisense > 0 #transform to logical

		multiple <- which( rowSums( olapsAntisense ) > 1 ) #logical comparison, > 1 if ambigous
		emdAnno[ whichNotMapped[ multiple ] ] <- 'Multiple'
		olapsAntisense[ multiple, ] <- FALSE

		categoryIdxAntisense <- apply( olapsAntisense, 2, which, TRUE ) #list of indices of clusters mapping to same feature
		emdAnno[ whichNotMapped[ categoryIdxAntisense$olapsCds ] ] <- 'CDS as'
		emdAnno[ whichNotMapped[ categoryIdxAntisense$olapsIx ] ] <- 'Introns as'
		emdAnno[ whichNotMapped[ categoryIdxAntisense$olapsTpUTR ] ] <- '3\' UTR as'
		emdAnno[ whichNotMapped[ categoryIdxAntisense$olapsFpUTR ] ] <- '5\' UTR as'
		emdAnno[ which( emdAnno == '' ) ] <- 'Other' #the remaining unmapped clusters are assigned here
		elementMetadata( clusters )[, 'MapsTo'] <- emdAnno 
	
		#7-count overlaps on the antisense strand by feature
		nNotMapped <- length( notMapped )
		summAntisense <- colSums( olapsAntisense )
		nMultiple <- length( multiple )
		summAntisense <- c( summAntisense, nMultiple, nNotMapped - sum( summAntisense ) - nMultiple )
	}
	else {
		summAntisense <- c( rep(0, 5 ), summ[ 6 ] )	#if all clusters map, just return zeros along with 'other' from sum
	}

	#8-plot the results as dotchart (ggplot)
	if( plot ) {
		names( summ ) <- c('CDS', 'Introns', '3\'-UTR', '5\'-UTR', 'Multiple', 'Other')
		names( summAntisense ) <- c('CDS', 'Introns', '3\'-UTR', '5\'-UTR', 'Multiple', 'Other')
		nSumm <- sum( summ )
		nSummAntisense <- sum( summAntisense )
		summ <- summ / nSumm * 100 
		summAntisense <- summAntisense / nSummAntisense * 100	
		
		summNorm <- summ[ 1 : 4 ] / proportionGenome
		summNorm <- summNorm / sum( summNorm ) * 100	

		summaryDf <- data.frame( group = factor(rep( 0 : 3, times = c( 6, 6, 4, 4) ), 
				   	  labels = c( 'Sense', 'Antisense', 'Transcriptome', 'Normalized') ),
				  	  Compartment = c( rep( names( summ ), 2), rep( names( summ )[ 1 : 4], 2) ),
  				  	  Percentage = c( summ, summAntisense, proportionGenome, summNorm ) )

		#to change facet grid labels appearance
		changeLab <- function( var, value ) {
    			value <- as.character( value )
    			if ( var == 'group' ) { 
        			value[ value == 'Sense'] <- paste0( 'Sense (n=', nSumm, ')' )
        			value[ value == 'Antisense'] <- paste0( 'Antisense (n=', nSummAntisense, ')' )
			}
    			return( value )
		}

		p <- ggplot( data = summaryDf,
			aes( x = Percentage, y = Compartment ) ) + 
 			geom_point( colour = 'royalblue', size = 2 ) + 
 			facet_grid( group ~ ., scales = 'free', labeller = changeLab ) +
			theme_bw()
			labs( title = 'Cluster annotation' ) +
			theme( plot.title = element_text( size = rel( 1 ) ) )
		print( p ) #need explict print
	}

	#10-return clusters with annotation as additional column in elementMetadata
	return( clusters )
}

