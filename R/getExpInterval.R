#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#2013 - Federico Comoglio & Cem Sievers, D-BSSE, ETH Zurich
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

getExpInterval <- function( model, bayes = TRUE, leftProb, rightProb, plot = TRUE ) {
# identify the RSF interval dominated by experimental induction and, if plot, returns model fit diagnostic plots
#
# Args:
#   model: list, mixture model as returned by fitMixtureModel
#   bayes: logical, if TRUE the Bayes criterion (cutoff at 0.5 posterior probabilities) is applied. If FALSE, custom cutoff values should be provided. Default is TRUE.
#   leftProb: numeric, the posterior probability cutoff to be applied to determine the start position of the support
#   rightProb: numeric, the posterior probability cutoff to be applied to determine the end position of the support
#   plot: logical, if TRUE diagnostics plot showing the model components, log odds and the computed posterior with highlighted identified support are returned
#
# Returns:
#   a list with two numeric slot, corresponding to the extreme of the support
#
# Error handling
#   if bayes FALSE and no custom cutoffs provided, raise an error

	if( !bayes & ( missing( leftProb ) | missing( rightProb ) ) ) {
		stop( 'bayes = FALSE and no left and right probability cutoffs provided. Please consider either using the Bayesian criterion or provide custom cutoffs.' )
	}
    	
	#1-assign left and right cutoffs according to user choice
	if( bayes ) {
		leftProb = .5
		rightProb = .5	
	}	
	else {
		leftProb = leftProb
		rightProb = rightProb
	}

	xval <- seq(.001, .999, .001)

	#2-compute log odds and responsibilities (for second component)
	odds <- computelogOdds( model )
	responsibilities <- (model$l2 * model$p2) / (model$l1 * model$p1 + model$l2 * model$p2)

	#3-compute support
	left <- which( responsibilities >= leftProb )[ 1 ]
	right <- which( responsibilities >= rightProb )
	right <- tail( right, 1 )
	supportStart <- xval[ left ]
	supportEnd <- xval[ right ]

	#4-if plot, plot model components, posterior and highlight support
	if( plot ) {
		par( mfrow = c( 1, 2 ) )
		plot(x     = xval, 
		     y     = model$p, 
		     type  = 'l', 
		     lwd   = 2,
		     ylim  = c( 0, max( model$p, odds ) ),
		     col   = 'black', 
		     ylab  = 'density', 
		     xlab  = 'Relative Substitution Frequency',
	 	     main  = 'Model densities' )

		lines( xval, model$p1, col = 'red', lwd = 2 )
		lines( xval, model$p2, col = 'blue3', lwd = 2 ) 
		lines( xval, odds, col = 'green3', lwd = 2 )

		legend( x = 'topright', 
			legend = c( 'p (full density)', 'p1 (non-experimental)', 'p2 (experimental)', 'log odds ratio' ), 
			col = c( 'black', 'red', 'blue3', 'green3'),
			lwd = 2,
			bty = 'n' )

		plot( x    = xval, 
		      y    = responsibilities, 
		      type = 'n',
		      ylim = c( -.1, 1 ), 
		      ylab = 'Posterior probability', 
                      xlab = 'Relative Substitution Frequency', 
                      main = 'Posterior class probability' )

		xCoord <- seq( supportStart, supportEnd, .001 )
		yCoord <- responsibilities[ left : right ]
		polyCoord <- cbind( xCoord, yCoord )
		polyCoord <- rbind( polyCoord, 
				    c( polyCoord[ nrow( polyCoord ), 1 ], 0 ),
				    c( polyCoord[1, 1], 0 ),
				    polyCoord[1, ] )

		polygon( polyCoord, col = 'royalblue1', border = NULL )
		lines( xval, responsibilities, lwd = 2 )
		rect( supportStart, -.1 , supportEnd, -.05, col = 'wheat2', border = NA )
		text( x      = ( supportEnd - supportStart ) / 2, 
		      y      = -.077, 
		      labels = paste( 'Support = [', supportStart, ', ', supportEnd, ']', sep = '' ) )
	}

	return( list( supportStart = supportStart, supportEnd = supportEnd ) )
}
