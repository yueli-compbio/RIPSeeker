# Function Name: 	galp2gal
# Description: 		Convert GAlignmentPairs to GAlignments using CIGAR to mark flanked portion of the pairs as 'N'					
# Input: 			GAlignmentPairs
# Output:			GAlignments
#
# Author: Yue Li
###############################################################################

galp2gal <- function(galp)
{	
	betweenPairCigar <- paste(abs(start(right(galp)) - end(left(galp)) + 1), "N", sep="")
	
	galcigar <- paste(cigar(left(galp)), betweenPairCigar, cigar(right(galp)), sep="")
  		
	gal <- GAlignments(
				seqnames = seqnames(galp),
			
				pos = start(left(galp)),
			
				cigar = galcigar,

				strand = strand(left(galp)),
			
				names = names(left(galp)),
				
				seqlengths = seqlengths(galp))
  
	# in case where end of a read exceeds chrom length
	# (which occurs for impropoer paired reads with large gap)	
	idx <- which(end(gal) < seqlengths(gal)[as.character(seqnames(gal))])	
  
  if(length(idx) < length(gal))    
    warning(sprintf("%d read-pairs with end larger than chromosome length are discarded", length(gal) - length(idx) + 1))
  
	return(gal[idx,])  	
}
