# Function Name: 	logScoreWithControl
# Description: 		x
# Input: 			y
# Output:			z
#
# Author: Yue Li
###############################################################################

logScoreWithControl <- function(nbhGRRIP, nbhGRCTL, padjMethod="BH", getControlStats=TRUE)
{	
	
	stopifnot(!missing(nbhGRRIP))
	
	stopifnot(!missing(nbhGRCTL))
			
	
	# significance score on logOdd scores
	logOddRIP <- computeLogOdd(nbhGRRIP)
	
	logOddCTL <- computeLogOdd(nbhGRCTL)
	
		
	logOddScore <- logOddRIP - logOddCTL
	
				
	if(existsFunction("mcols_inefficient")) values(nbhGRRIP) <- cbind(mcols(nbhGRRIP), logOddScore) else values(nbhGRRIP) <- cbind(as.data.frame(values(nbhGRRIP)), logOddScore)
	
	
	enrichIdx <- which(values(nbhGRRIP)$viterbi_state == 2)
	
	# assign 1 as dummy index if no enriched state present
	if(length(enrichIdx) == 0) enrichIdx <- 1	
	
	unmergedRIP <- nbhGRRIP[enrichIdx]
	
	
	mergedRIP <- reduce(unmergedRIP, 
			min.gapwidth = median(width(unmergedRIP) ))
	
	
	overlapIdx <- findOverlaps(mergedRIP, unmergedRIP)
	
	
	mergedRIPList <- lapply(split(overlapIdx, queryHits(overlapIdx)),
			scoreMergedBins, unmergedRIP, mergedRIP)
			
		
	names(mergedRIPList) <- NULL
	
	
	mergedRIP <- do.call(c, mergedRIPList)

		
	# compute p-value of the averaged logOddScore (estimate mean with all bins)
	pval <- pnorm(values(mergedRIP)$logOddScore, mean(logOddScore), 
			sd(logOddScore), lower.tail=FALSE)
  
	# compute p-value of the averaged logOddScore (estimate mean with only merged bins)
# 	pval <- pnorm(values(mergedRIP)$logOddScore,
#                 mean(values(mergedRIP)$logOddScore),
#                 sd(values(mergedRIP)$logOddScore), lower.tail=FALSE)
  
				
	pvalAdj <- p.adjust(pval, padjMethod)
	
		
	if(existsFunction("mcols_inefficient")) values(mergedRIP) <- cbind(mcols(mergedRIP), pval, pvalAdj) else values(mergedRIP) <- cbind(as.data.frame(values(mergedRIP)), pval, pvalAdj)
			
	
	if(getControlStats) {
		
		unmergedCTL <- nbhGRCTL[enrichIdx]
		
		logOddScore <- logOddCTL[enrichIdx]
		
		
		if(existsFunction("mcols_inefficient")) values(unmergedCTL) <- cbind(mcols(unmergedCTL), logOddScore) else values(unmergedCTL) <- cbind(as.data.frame(values(unmergedCTL)), logOddScore) 
		
								
		mergedCTLList <- lapply(split(overlapIdx, queryHits(overlapIdx)),
				scoreMergedBins, unmergedCTL, mergedRIP)
		
		names(mergedCTLList) <- NULL
		
		mergedCTL <- do.call(c, mergedCTLList)
		
		
		RIPcolnames <- names(values(mergedRIP))
		
		
		CTLcolnames <- paste("CTL", names(values(mergedCTL)), sep="_")
						
	
		if(existsFunction("mcols_inefficient")) values(mergedRIP) <- cbind(mcols(mergedRIP), mcols(mergedCTL)) else values(mergedRIP) <- cbind(as.data.frame(values(mergedRIP)), 
					as.data.frame(values(mergedCTL)))
				
		names(values(mergedRIP)) <- c(RIPcolnames, CTLcolnames)
	}
	
	
	return(mergedRIP)	
}
