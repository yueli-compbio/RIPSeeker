# Function Name:	selectBinSize
# Description:		Convert alignment file in BAM or BED format to GAlignments obejct;
# 					If BED file is provided, genomeBuild is needed to get chromosome lengths;
# 					GAlignments object is intended for read coverage
# Input: 			alignFilePath etc
# Output:			GAlignments object 
#
# Author: Yue Li
###############################################################################


getAlignGal <- function(alignFilePath, format, genomeBuild, 
				deleteGeneratedBAM=FALSE, reverseComplement=FALSE, 
				returnDuplicate=FALSE, flagMultiHits=TRUE, 
				returnOnlyUniqueHits=FALSE, paired=FALSE, ...)
{	
	
	# private function used to read single-end or paired-end data 
	getAlignGalPrivate <- function(alignFilePath, paired, param) {
		
		if(!paired) alignGal <- readGAlignments(alignFilePath, use.names=TRUE, param=param)
		
		if(paired) alignGal <- galp2gal(readGAlignmentPairs(alignFilePath, use.names=TRUE, param=param))
		
		return(alignGal)
	}
			
	SUPPORTED_FORMAT <- c("bam", "sam", "bed")
	
	stopifnot(!missing(alignFilePath))
	
	if(missing(format)) {		
		# determine from file extension		
		parts <- unlist(strsplit(alignFilePath, "\\."))
		
		format <- parts[length(parts)]
		
	}
	
	format <- tolower(format)
	
	if(!format %in% SUPPORTED_FORMAT) {
			
		warning("File extension or format does not look like bam, sam, or bed. Try import anyways.")
	}

	
	# set parameters for initial filtering of bam file	
	if(paired) { # paired-end
		
		param <- ScanBamParam(flag=scanBamFlag(isUnmappedQuery=FALSE, isPaired=paired, 
						
						hasUnmappedMate=FALSE, isProperPair=paired, isDuplicate=returnDuplicate))		
		
	} else { # single-end
		
		param <- ScanBamParam(flag=scanBamFlag(isUnmappedQuery=FALSE, isDuplicate=returnDuplicate))
	}

	message(sprintf("Processing %s ... ", alignFilePath), appendLF = FALSE)
	
	if(format == "bam") {
		
		alignGal <- getAlignGalPrivate(alignFilePath, paired, param)
	}
	
	if(format == "sam") {
		
		savedBAM <- sub("\\.sam", "", alignFilePath)
		
		BAMfile <- paste(savedBAM,".bam",sep="")
		
		if(file.exists(BAMfile)) {
			
			stop(sprintf("Looks like %s already exists. Try use it or delete it first.", BAMfile))
		}
		
		message("Convert SAM to BAM ...")
		
		# need Rsamtools library
		asBam(file=alignFilePath, destination=savedBAM, overwrite=FALSE, indexDestination=FALSE)
		
		message(sprintf("BAM file was successfully converted and saved as:\n %s.bam", savedBAM))
		
		alignGal <- getAlignGalPrivate(alignFilePath, paired, param)
				
		if(deleteGeneratedBAM) system(sprintf("rm %s", BAMfile))
		
		message(sprintf("Returned alignGal and deleted %s", BAMfile))
	}
	
	if(format == "bed" || !(format %in% SUPPORTED_FORMAT)) {
		
		if(missing(genomeBuild)) {
			
			stop("genomeBuild is needed to determine chromosome lengths for the GAlignments obejct")
			
		} else { # determien chromosome length from genomeBuild
						
			# need rtracklayer library
			chromInfo <- SeqinfoForUCSCGenome(genomeBuild)
			
			alignGR <- import(file) # convert into GRanges first
			
			if(!returnDuplicate) alignGR <- unique(alignGR)
						
			# construct GAlignments object			
			alignGal <- GAlignments(names = values(alignGR)$name,
						seqnames = seqnames(alignGR),
			 			pos = start(alignGR), strand = strand(alignGR),
					 	cigar=paste(width(alignGR), "M", sep=""))
					 	
			# based on chrom names in alignGR to find their lengths in seqInfo
			# then assign the chrom lengths to alignGR
			seqlengths(alignGal) <- as.data.frame(chromInfo[names(seqlengths(alignGal))])$seqlengths			
			
		}				
	}
	
	# if reads come from the complementary strand, signs needs to be switched
	if(reverseComplement) { 
		
		strandDummy <- runValue(strand(alignGal))
		strandDummy <- sub("\\+", "dummy", strandDummy)
		strandDummy <- sub("\\-", "\\+", strandDummy)
		runValue(strand(alignGal)) <- sub("dummy", "\\-", strandDummy)				
	}
			
	if(!missing(genomeBuild)) genome(alignGal) <- genomeBuild
	
	
	# record processing options in metadata slot
	metadata(alignGal) <- list(reverseComplement=reverseComplement, returnDuplicate=returnDuplicate,
			flagMultiHits=flagMultiHits, returnOnlyUniqueHits=returnOnlyUniqueHits)
	
		
	if((returnOnlyUniqueHits || flagMultiHits) && is.null(names(alignGal))) {
		
		warning(sprintf("No read name provided so cannot detect multihits.\nResults returned without multihit flags"))
		
		return(alignGal)		
	}
	
	if(!flagMultiHits && !returnOnlyUniqueHits) {
		
		message("All reads are returned without flagging multi or unique hits.")
				
		return(alignGal)
	}
		
	
	if(flagMultiHits || returnOnlyUniqueHits) {
				
		readNames <- names(alignGal)
		
		# sort and count read names
		readNamesRle <- Rle(sort(readNames))
		
		# get read names with a single occurence
		nonMultiHitReadNames <- runValue(readNamesRle)[runLength(readNamesRle) == 1]
		
		multiHitReadNames <- runValue(readNamesRle)[runLength(readNamesRle) > 1]
		
		
		# return results
		if(!flagMultiHits && returnOnlyUniqueHits) {
			
			message("Only unique hits are returned without flag.")
			
			return(alignGal[match(nonMultiHitReadNames, readNames)])
		}
				
		
		if(flagMultiHits) {
			# flag multihits (if any)
			values(alignGal) <- FALSE
			
			names(values(alignGal)) <- "uniqueHit"
			
			values(alignGal)[match(nonMultiHitReadNames, readNames), ] <- TRUE
			
			if(length(multiHitReadNames) > 0) {
				values(alignGal)[match(multiHitReadNames, readNames), ] <- FALSE
			}
			
			if(returnOnlyUniqueHits) {
				
				message("Only unique hits are returned with (trivial) flags.")
				
				return(alignGal[match(nonMultiHitReadNames, readNames)])
			} else {
				
				message("All hits are returned with flags.")
				
				return(alignGal)				
			}
		}
		
	}
		
	
}
