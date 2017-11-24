# Function Name: 	mainSeekMC, multicore version of mainSeek using mclapply
# Description: 		runs main for RIP or CTL bam file(s); treat multiple bams as replicate and combine them; actual RIP regions will be predicted in seekRIP outside of this function
# Input: 			bamPath to bam files (either a single dir or the multiple file paths) 
# Output:			GAlignments and GRanges objects representating alignment object after disambiguating multihits (if any) and predicted RIP bins (not regions yet)
#
# Author: Yue Li
###############################################################################

mainSeek <- function(bamFiles, reverseComplement=FALSE, genomeBuild="mm9",
		uniqueHit = TRUE, assignMultihits = TRUE, strandType = NULL,
		paired=FALSE, rerunWithDisambiguatedMultihits = TRUE, silentMain=FALSE,
		multicore = TRUE, returnAllResults = TRUE, ...)	
{	
	if(multicore) {
		
		hasMC <- try(require(parallel))
		
		if(!hasMC) {
			
			warning("multicore package is unavailable; use single-core instead.")
			
			multicore <- FALSE
		}
	}		
	
	stage <- 1
	
	################ Read in BAM ################
	# by default, return all alignments but flag multihits for
	# subsequent subsetByMaxPostprobOverlaps
	
	message(sprintf("\n**%s. Process and combine alignment files\n", LETTERS[stage]))
	
	alignGal <- combineAlignGals(bamFiles, reverseComplement=reverseComplement, genomeBuild=genomeBuild, paired=paired)
	
	if(!is.null(strandType)) {
		
		if(!strandType %in% levels(strand(alignGal))) {
						
			warning(sprintf("*** Strand type (%s) is unknown and ignored.", strandType))
			
		} else {
			
			message(sprintf("*** Only reads from strand %s will be considered.", strandType))
			
			alignGal <- alignGal[strand(alignGal) == strandType]
		}		
	}
	
	
	################ Return unique hits only ################
	if(uniqueHit) {
		
		if(is.logical(values(alignGal)$uniqueHit)) {
			
			message("*** Only unique hits are used to compute read count.")
			
			alignGR <- as(alignGal[values(alignGal)$uniqueHit], "GRanges")
			
		} else {			
			warning("*** uniqueHit is TRUE but uniqueHit flag slot is not set!\nAll alignment will be treated as unique hits;\notherwise try run getAlignGal(..., flagMultiHits=TRUE).")
			
			alignGR <- as(alignGal, "GRanges")
		}		
	} else {
		
		message("*** All reads (including multihits) are used to construct initial HMM model.")
		
		alignGR <- as(alignGal, "GRanges")
	}
	
	################ Add pseudo count ################
	# checking for zero alignment in some chromosome (e.g. chrM)
	alignGR <- addPseudoAlignment(alignGR)
	
	stage <- stage + 1
	
	
	################ start main function to train HMM ################
	
	runViterbi <- all(values(alignGal)$uniqueHit) || !rerunWithDisambiguatedMultihits
	
	message(sprintf("\n**%s. Run NB HMM on each chromosome\n", LETTERS[stage]))
	
	
	if(multicore) {
		
		if(silentMain) {
			suppressMessages(
					nbhGRList <- mclapply(as.list(split(alignGR, seqnames(alignGR))),
							mainSeekSingleChrom, runViterbi = runViterbi, ...))
		} else {			
			nbhGRList <- mclapply(as.list(split(alignGR, seqnames(alignGR))),
					mainSeekSingleChrom, runViterbi = runViterbi, ...)
		}
		
	} else {
		
		if(silentMain) {
			suppressMessages(
					nbhGRList <- lapply(as.list(split(alignGR, seqnames(alignGR))),
							mainSeekSingleChrom, runViterbi = runViterbi, ...))
		} else {			
			nbhGRList <- lapply(as.list(split(alignGR, seqnames(alignGR))),
					mainSeekSingleChrom, runViterbi = runViterbi, ...)
		}
	}
	
	nbhGRList <- GRangesList(nbhGRList)		
	
	tmp <- gc(); rm(tmp)
	
	stage <- stage + 1	
	
	################ Disambiguate multihit if exists ################
	if(!all(values(alignGal)$uniqueHit) && uniqueHit && assignMultihits) {
		
		message(sprintf("\n**%s. Disambiguate multihits based on posterior\n", LETTERS[stage]))
		
		alignGalFiltered <- disambiguateMultihits(alignGal, nbhGRList)
		
		stage <- stage + 1		
		
		if(rerunWithDisambiguatedMultihits) {
			
			message(sprintf("\n**%s. Re-run NB HMM with unique hits + disambiguated multihits.\n",
							LETTERS[stage]))
			
			
			################ Re-run HMM on augmented library ################
			# re-estimating HMM parameters based on augmented read library
			alignGR <- as(alignGalFiltered, "GRanges")
			
			
			# checking for zero alignment in some chromosome (e.g. chrM)
			alignGR <- addPseudoAlignment(alignGR)
			
			
			# re-estimate parameters using augmented reads corrected with
			# disambiguateMultihits function
			# enable viterbi prediction
			
			if(multicore) {
				
				if(silentMain) {
					suppressMessages(
							nbhGRList <- mclapply(as.list(split(alignGR, seqnames(alignGR))), 
									mainSeekSingleChrom, runViterbi = TRUE, ...))
				} else {
					nbhGRList <- mclapply(as.list(split(alignGR, seqnames(alignGR))), 
							mainSeekSingleChrom, runViterbi = TRUE, ...)					
				}
				
			} else {
				
				if(silentMain) {
					suppressMessages(
							nbhGRList <- lapply(as.list(split(alignGR, seqnames(alignGR))), 
									mainSeekSingleChrom, runViterbi = TRUE, ...))
				} else {
					nbhGRList <- lapply(as.list(split(alignGR, seqnames(alignGR))), 
							mainSeekSingleChrom, runViterbi = TRUE, ...)					
				}
			}
						
			
			nbhGRList <- GRangesList(nbhGRList)            
		}		
		
		
		output <- list(nbhGRList=nbhGRList, alignGal=alignGal, alignGalFiltered=alignGalFiltered)
		
	} else {
		
		output <- list(nbhGRList=nbhGRList, alignGal=alignGal)
	}
		
	# assign strandType to the strand RIP regions on each chromosome
	if(!is.null(strandType)) output$nbhGRList <- endoapply(output$nbhGRList,
				function(gr, strandType){strand(gr) <- strandType; gr}, strandType)
									
		
	################ return HMM output (alignGal + alignGalFiltered) ################
	
	if(returnAllResults) return(output)
	if(!returnAllResults) return(output$nbhGRList)		
}
