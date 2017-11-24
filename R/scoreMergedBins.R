# Function Name: 	scoreMergedBins
# Description: 		Given findOverlaps output (2-column indices) of overlaps hit 
#					b/w merged enriched Region (query) and unmerged enriched bins (subjects),
#					average over the logOdd score across bins within the merged region.
#					This function is expected to be called ONLY within seekRIP.
# Input: 			y
# Output:			z
#
# Author: Yue Li
###############################################################################

scoreMergedBins <- function(findOverlapsHits, unmergedGRAll, mergedGRAll) {
				
	# get unmerged regions (bins) based on the subjectHits index
	unmergedGR <- unmergedGRAll[subjectHits(findOverlapsHits)]
		
	
	# the same merge region appear N times, as the number of merged regions
	# thus only one index is needed
	mergedGR <- mergedGRAll[queryHits(findOverlapsHits)[1]]
			
	# summarized values
	count <- sum(values(unmergedGR)$count)
		
	fpk <- 1000 * count/(sum(width(unmergedGR)))
		
	logOddScore <- mean(values(unmergedGR)$logOddScore)
	
								
	# save summarized values
	values(mergedGR) <-  
			data.frame(count=count, fpk=fpk, logOddScore=logOddScore)	
		
	return(mergedGR)
}
