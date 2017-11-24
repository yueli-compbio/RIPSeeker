# Function Name: 	f
# Description: 		x
# Input: 			y
# Output:			z
#
# Author: Yue Li
###############################################################################
addPseudoAlignment <- function(alignGR)
{
	emptyChrnames <- names(which(
					unlist(lapply(split(alignGR, seqnames(alignGR)), length)) == 0))
	
	if(length(emptyChrnames) > 0) {
		
		message(sprintf("*** %s do not have any alignment.\n", emptyChrnames))
		
	} else {
		
		message(sprintf("*** All chromosomes have at least one alignment"))
		
		return(alignGR)
	}
	
	pseudoreads <- rep(alignGR[1], length(emptyChrnames))
	
	start(pseudoreads) <- 1
	
	end(pseudoreads) <- 20
	
	emptyChrnames <- factor(emptyChrnames, levels=levels(seqnames(pseudoreads)))
	
	
	seqnames(pseudoreads) <- Rle(emptyChrnames)
	
	
	names(pseudoreads) <- paste("pseudoreads", 1:length(pseudoreads), sep="")
	
	
	alignGR <- append(alignGR, pseudoreads)
	
	message(sprintf("*** %d pseudoreads are appended to the end.", length(emptyChrnames)))
	
	return(alignGR)
}