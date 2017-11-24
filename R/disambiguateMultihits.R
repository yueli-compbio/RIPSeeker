# Function Name: 	disambiguateMultihits
# Description: 		among multiple alignments of the same read (i.e. multihit),
#					select the alignment corresponding to the bin with maximum
#					posterior for the enriched hidden state.
# Input: 			GAlignments with 1-read-(>=1)-alignments
#					posterior decoding from HMM (GRanges) 
# Output:			the same GAlignments with 1-read-1-alignment 
#
# Author: Yue Li
###############################################################################

disambiguateMultihits <- function(alignGal, nbhGRList, postprobCutoff=0)
{
	# unlist GRangesList to search multihits over all chromosomes
	nbhGR <- unlist(nbhGRList)
	
	names(nbhGR) <- NULL
	
	
	# find indices of multihits labeled as FALSE for column "uniqueHit"
	# labelling is expected to have been performed by getAlignGal
	idx <- which(values(alignGal)$uniqueHit == FALSE)
	
	
	# subset alignGal multihit reads
	mhitReads <- alignGal[idx]
	
		
	# for each multihit, get indices of the corresponding genomic regions from nbhGR
	mhitAllOverlaps <- findOverlaps(mhitReads, nbhGR)
	
	
	# index multireads in the order as nbhGR interesected regions occur
	mhitQuery <- mhitReads[queryHits(mhitAllOverlaps)]
	
	# then assign enriched state posterior to each alignment				
	values(mhitQuery) <- values(nbhGR[subjectHits(mhitAllOverlaps)])$state_2_postprob
	
	
	names(values(mhitQuery)) <- "state_2_postprob"

			
	mhitMatrix <- cbind(names(mhitQuery), 1:length(mhitQuery), 							
			values(nbhGR[subjectHits(mhitAllOverlaps)])$state_2_postprob)
	
	
	mhitList <- split(mhitMatrix, mhitMatrix[,1])
	
	
	# central step: split alignment by read names into a list
	# for each list element (containing > 1 alignments for the same read)
	# select the alignment corresponding to the region with max postprob
	# return only the unique index of that alignment in order to unlist the results
	
	desiredIdx <- lapply(mhitList,
			function(x) {
				
				y <- matrix(x, ncol=3)
				
				maxIdx <- which.max(y[,3])
				
				if(max(y[maxIdx,3]) > postprobCutoff) as.numeric(y[maxIdx, 2])
			}
	)
	
	desiredIdx <- unlist(desiredIdx)
	
	
	selectGal <- mhitQuery[desiredIdx, ]
	
	
	values(selectGal) <- FALSE
	
	
	names(values(selectGal)) <- "uniqueHit"
		
		
	alignGalfiltered <- c(alignGal[values(alignGal)$uniqueHit == TRUE], selectGal)
	
	
	message(sprintf("\n%d/%d multihit reads corresponding to %d ambiguous alignments\nhave been assigned to %d unique regions with maximum posterior for the enriched state\n", 
					length(selectGal), length(alignGalfiltered), 
					length(idx), length(selectGal)))
	
		
			
	return(alignGalfiltered)			
}
