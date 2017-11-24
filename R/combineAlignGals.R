# Function Name: 	combineAlignGals
# Description: 		Read each BAM/SAM/BED using getAlignGal.R and combine them 
#					(e.g. technical/biological replicates).					
# Input: 			a list or vector of bam file names in full path
# Output:			combined GAlignments objects
#
# Author: Yue Li
###############################################################################

combineAlignGals <- function(bamFiles, ...)
{	
	
	numBAM <- length(bamFiles)
		
	combinedGal <- getAlignGal(bamFiles[1], ...)
			
	if(numBAM == 1) {
	
		message(sprintf("%d BAM files are combined", numBAM))
		return(combinedGal)
	}
					
	for(i in 2:numBAM) {
		
		combinedGal <- 
				append(combinedGal, getAlignGal(bamFiles[i], ...))
	}
	
	message(sprintf("%d BAM files are combined", numBAM))
	return(combinedGal)
}
