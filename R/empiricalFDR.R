# Function Name: 	f
# Description: 		x
# Input: 			y
# Output:			z
#
# Author: Yue Li
###############################################################################

empiricalFDR <- function(pval, pvalRIP, pvalCTL)
{
	stopifnot(!missing(pval))
	stopifnot(!missing(pvalRIP))
	stopifnot(!missing(pvalCTL))
	
	falseCount <- length(which(pvalCTL <= pval))
	
	trueCount <- length(which(pvalRIP <= pval))
	
	# control eFDR above pval esp. when it becomes zero (i.e. no control-RIP
	# comparison has pval as significant as the one found in RIP-control 
	eFDR <- max(falseCount / trueCount, pval)
	
	
	# if (false count)/(true count) > 1, then define the FDR as 1
	return(min(1, eFDR))
}
