# Function Name: 	computeRIPScore
# Description: 		x
# Input: 			y
# Output:			z
#
# Author: Yue Li
###############################################################################

computeLogOdd <- function(nbhGR)
{	
	stopifnot(!missing(nbhGR))
	
	stopifnot(!is.null(values(nbhGR)$state_1_postprob))
	
	stopifnot(!is.null(values(nbhGR)$state_2_postprob))
				
	logOdd <- log(values(nbhGR)$state_2_postprob) - log(values(nbhGR)$state_1_postprob)
	
	# in case one of the state has 0 probability, assign it log of teh largest/smallest value
	logOdd[which(logOdd == -Inf)] <- log(.Machine$double.xmin)
	
	logOdd[which(logOdd == Inf)] <- log(.Machine$double.xmax)	
	
	return(logOdd)
}