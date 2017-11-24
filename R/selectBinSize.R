# Function Name: 	selectBinSize
# Description: 		Select optimal bin size based on Shimazaki formula
# Input: 			GenomicRanges object, <minBinSize>, <maxBinSize>, <increment>
# Output:			Optimal bin size estimated over chromosomes 
#
# Author: Yue Li
###############################################################################

selectBinSize <- function(alignGR, minBinSize, maxBinSize = 1000, 
		increment=5, getFullResults=FALSE)
{			
	stopifnot(!missing(alignGR))
		
	readlen <- median(width(ranges(alignGR)))
	
	# set min bin size
	if(missing(minBinSize)) minBinSize <- 5 * readlen
	
	
	# set max bin size
	uniqChr <- unique(seqnames(alignGR))
	
	tmp <- maxBinSize
	
	maxBinSize <- min(maxBinSize, seqlengths(alignGR)[uniqChr])
		
	stopifnot(minBinSize <= maxBinSize)
	
				
	if(tmp != maxBinSize) {
						
		warnings(sprintf("Specified/default maxBinSize (%d) is larger than the total sequence length (%d); using sequence length as max bin size instead", 
						maxBinSize, seqlengths(alignGR)[uniqChr]))							
	}
		
	
	binSizes <- seq(minBinSize, maxBinSize, by=increment)
		
	costs <- numeric(length(binSizes))
									
	costs <- lapply(binSizes, evalBinSize, alignGR)
		
	costs <- as.numeric(unlist(costs))
	
	idx <- which.min(costs)
		
	
	if(getFullResults) { # return all results
		list(bestBinSize = binSizes[idx], minCosts = costs[idx],
				costs=costs, binSizes=binSizes)
	} else { 			# default: only return the best bin size
		binSizes[idx]
	}
}
