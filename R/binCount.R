# Function Name: 	binCount
# Description: 		Count reads in nonoverlapping bins across the chromosome
# Input: 			GenomicRanges object, binSize, <seqIRangesList>
# Output:			Vector of bin counts corresponding to nonoverlapping bins
#
# Author: Yue Li
###############################################################################

binCount <- function(alignGR, binSize, returnBinCountOnly=FALSE)
{	
	
	stopifnot(!missing(alignGR))
	stopifnot(!missing(binSize))
	
	
	chrname <- as.character(runValue(seqnames(alignGR)))
	
	seqLen <- seqlengths(alignGR)[chrname]

		
	# Start coordinates of running nonoverlapping bins
	binStarts <- seq(from=1, to=seqLen, by=binSize)
	
	
	# End coordinates of running nonoverlapping bins;
	binEnds <- seq(from=1, to=seqLen, by=binSize) + binSize - 1
	
	# Set the end coordinate of the last bin to sequence length
	binEnds <- c(binEnds[-length(binEnds)], seqLen)
		
	
	# Create GRanges using nonoverlapping bin coordinates above
	binGR <- GRanges(	seqnames=chrname, 
						ranges = IRanges(start=binStarts, end=binEnds),
						seqlengths=seqLen)
				
	binCnt <- countOverlaps(binGR, resize(alignGR, 1))
			
	if(returnBinCountOnly) {
		
		return(binCnt)
		
	} else {
		
		values(binGR) <- binCnt
		
		names(values(binGR)) <- "count"
		
		return(binGR)
	}
}
