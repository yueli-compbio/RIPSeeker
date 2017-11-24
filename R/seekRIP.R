# Function Name: 	seekRIP
# Description: 		x
# Input: 			y
# Output:			z
#
# Author: Yue Li
###############################################################################

seekRIP <- function(nbhGRRIP, nbhGRCTL=NULL, padjMethod="BH", logOddCutoff=-Inf,
		pvalCutoff=1, pvalAdjCutoff=1, eFDRCutoff=1)
{
	
	stopifnot(!missing(nbhGRRIP))
	
	if(is.null(nbhGRCTL)) {
		
		# compute logOddScore as difference b/w logOdd in RIP and CTL
		mergedRIP <- logScoreWithoutControl(nbhGRRIP, padjMethod)		
		
	} else {
		
		if(is.list(nbhGRCTL) || class(nbhGRCTL) == "GRangesList") {
			
			chrname <- as.character(runValue(seqnames(nbhGRRIP)))
			
			nbhGRCTL <- nbhGRCTL[[chrname]]									
		}
		
		
		# compute logOddScore as difference b/w logOdd in RIP and CTL
		mergedRIP <- logScoreWithControl(nbhGRRIP, nbhGRCTL, padjMethod)
						
		
		# swap RIP and CTL to compute empirical FDR
		mergedCTL <- logScoreWithControl(nbhGRCTL, nbhGRRIP, 
				padjMethod="none", getControlStats = FALSE)
		
		
		eFDR <- sapply(values(mergedRIP)$pval, empiricalFDR, 
				values(mergedRIP)$pval, values(mergedCTL)$pval)
		
		
		if(existsFunction("mcols_inefficient")) values(mergedRIP) <- cbind(mcols(mergedRIP), eFDR) else values(mergedRIP) <- cbind(as.data.frame(values(mergedRIP)), eFDR)
	}
		
	# construct cutoff rule
	cutoff <- 
			values(mergedRIP)$logOddScore >= logOddCutoff &
			values(mergedRIP)$pval <= pvalCutoff &
			values(mergedRIP)$pvalAdj <= pvalAdjCutoff
	
	
	if(!is.null(nbhGRCTL)) cutoff <- cutoff & values(mergedRIP)$eFDR < eFDRCutoff
	
	names(mergedRIP) <- NULL
	
	return(mergedRIP[cutoff])				
}