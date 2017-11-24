# Function Name: 	mainSeekSingleChrom
# Description: 		x
# Input: 			y
# Output:			z
#
# Author: Yue Li
###############################################################################

mainSeekSingleChrom <- function(alignGR, K=2, binSize = NULL, minReadCount = 10, backupNumBins=10,
		minBinSize = 200, maxBinSize = 1200, increment=5, pathToSavePlotsOfBinSizesVersusCosts, 
		verbose = TRUE, allowSecondAttempt = TRUE, ...)
{	
	stopifnot(!missing(alignGR))
	
	tmp <- gc(); rm(tmp);
	
	useDefinedBinSize <- !is.null(binSize)
							
	chrname <- as.character(runValue(seqnames(alignGR)))
	
	message(sprintf("\n**********\n%s:\n**********\n",chrname))
	
	if(!is.null(binSize) && 
			(!is.numeric(unlist(binSize)) || any(unlist(binSize) < 0) || 
				any(unlist(binSize) != round(unlist(binSize))))
	) {		
			message(print(binSize))
			stop("binSize is/are not positive integers.")		
	}
		
	
	# check whether binSize is a list/vector that contain specific
	# bin sizes for individual chromosomes.
	# Such binSize list/vector *must* has chromosome names the same 
	# as the seqnames in alignGR
	if(!is.null(binSize) && length(binSize) > 1) {
						
		if(is.list(binSize)) {
			tmp <- binSize[[chrname]]
		} else {
			tmp <- binSize[chrname]
		}
		
		if(is.null(tmp) || !is.numeric(tmp) || tmp >= seqlengths(alignGR)[chrname]) {
			
			warning("Invalid bin size stored in the metadata of GRanges object")				
		}				
		
		binSize <- tmp
		
		message(sprintf("***0. Using predefined binSize of %d bp.", binSize))
		
	} else {
		message(sprintf("***0. Using predefined binSize of %d bp.", binSize))
	}
					
	# pseudoreads is a place holder created earlier by addPsuedoAlignment.R
	pseudoreadscnt <- length(grep("pseudoreads", names(alignGR)))
	
	if(pseudoreadscnt > 0 || length(alignGR) < minReadCount ) {
		
		warning(sprintf("Not enough (< %d) alignment found in this chromosome:",
						minReadCount))
				
		if(is.null(binSize) || floor(seqlengths(alignGR)[chrname]/binSize) < backupNumBins) {	
			
			binSize <- floor(seqlengths(alignGR)[chrname]/backupNumBins)
			
			warning(sprintf("binSize has been arbitrarily set to %d (i.e. 1/%d of the total %s length",
							binSize, backupNumBins, seqlengths(alignGR)[chrname], chrname))			
		}								
		
		nbhGR <- addDummyProb(alignGR, K=K, binSize, randomProb=FALSE, runViterbi=list(...)$runViterbi)
		
		return(nbhGR)
	}
	
		
	# finally if binSize is null, then use binSize selection
	if(is.null(binSize)) {
	
		message(sprintf("***0. Computing optimal bin size."))
										
		b <- selectBinSize(alignGR, minBinSize=minBinSize, maxBinSize = maxBinSize, increment=increment,
				getFullResults=!missing(pathToSavePlotsOfBinSizesVersusCosts))
		
		if(!missing(pathToSavePlotsOfBinSizesVersusCosts) && is.list(b)) {
			binSize <- b$bestBinSize
			
			# plot costs as a function of bin sizes 
			pdf(file=sprintf("%s%s.pdf", pathToSavePlotsOfBinSizesVersusCosts, chrname))		
				chrlen <- seqlengths(alignGR)[chrname]
				plot(b$binSizes, b$costs)
				legend("topright", box.lty=0,
						sprintf("%s: 1-%d;\nTotal mapped reads: %d;\nOptimal bin size = %d bp",
								chrname, chrlen, length(alignGR), binSize))
			dev.off()
			
		} else {
			binSize <- b
		}				
		
		message(sprintf("\nOptimal bin size: %d bp\n", binSize))
	}
	
	
	# if binsize is initiallly undefined and is determined to be too large 
	# s.t. there are fewer than backupNumBins (default: 10)
	# then set the bin size to a smaller size to ensure backupNumBins of bins
	if(floor(seqlengths(alignGR)[chrname]/binSize) < backupNumBins) 	{
		
		oldBinSize <- binSize
		
		binSize <- floor(seqlengths(alignGR)[chrname]/backupNumBins)		
		
		warning(sprintf("binSize %d is too large for %s (%d).\nAnd change binSize to %d (i.e. 1/%d of the total chromosome length",
						oldBinSize, chrname, seqlengths(alignGR)[chrname], binSize, backupNumBins)
		)
		
		message(sprintf("\nBin size has been reset to: %d bp\n", binSize))				
	}
	
	
	# Compute Coverage within bins and run HMM (nbh) on bin count vector
	message(sprintf("***1. Traning NB HMM to derive posterior (and Viterbi state sequence:)"))
		

	if(verbose) nbhGR <- nbh(binCount(alignGR, binSize=binSize), K=K, ...)
	if(!verbose) suppressMessages(nbhGR <- nbh(binCount(alignGR, binSize=binSize), K=K, ...))
	
		
	if(class(nbhGR) != "GRanges" && allowSecondAttempt && !useDefinedBinSize) {
		
		i <- 1
		
		warning(sprintf("Main function resulted in None GRanges object on %s.\nThis might be due to bin size. Attempting with different best bin sizes."), chrname)
		
		b <- selectBinSize(alignGR, minBinSize=minBinSize, maxBinSize = maxBinSize, 
				increment=increment, getFullResults = TRUE)
		
		while(class(nbhGR) != "GRanges") {
			binSize <- b$binSizes[order(b$cost)[i]]
			
			# Compute Coverage within bins and run HMM (nbh) on bin count vector
			if(verbose) nbhGR <- nbh(binCount(alignGR, binSize=binSize), K=K, ...)
			if(!verbose) suppressMessages(nbhGR <- nbh(binCount(alignGR, binSize=binSize), K=K, ...))
			
			i <- i + 1
		}
	}
		
	# if all bets are off, then use dummy value with uniform posterior 
	if(class(nbhGR) != "GRanges")
	{
		warning("HMM failed on this chromosome (after trying all bin sizes). Default uniform posterior is used instead")
		
		nbhGR <- addDummyProb(alignGR, K=K, binSize, randomProb=FALSE, runViterbi=list(...)$runViterbi)			
	}
	
							
	return(nbhGR)		
}