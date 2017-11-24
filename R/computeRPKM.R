computeRPKM <- function(bamFiles, RIPSeekerRead=TRUE, paired=FALSE, 
			countMode="IntersectionNotEmpty", featureGRanges, 
			idType="ensembl_transcript_id", featureType="exon", 
			ignore.strand=FALSE, txDbName="biomart", 
			moreGeneInfo=FALSE, saveData, justRPKM=TRUE, ...)
{	
	library(Rsamtools)
	library(GenomicFeatures)
	library(biomaRt)
					
	##################### Make transcript db #####################
	
	if(missing(featureGRanges)) {
		
		stopifnot(txDbName %in% c("biomart", "UCSC"))
		
		stopifnot(featureType %in% c("exon", "intron", "fiveUTR", "threeUTR", "CDS"))
		
		stopifnot(list(...)$by %in% c("tx", "gene"))
		
		# hack functions to ignore unused arguments
		formals(makeTxDbFromBiomart) <- c(formals(makeTxDbFromBiomart), alist(... = ))
		formals(makeTxDbFromUCSC) <- c(formals(makeTxDbFromUCSC), alist(... = ))				
										
		if(txDbName == "biomart") txDb <- makeTxDbFromBiomart(...)
		
		if(txDbName == "UCSC") txDb <- makeTxDbFromUCSC(...)
						
		# group annotation by transcripts ("tx") or gene
		if(list(...)$by == "tx") {
			
			if(featureType == "exon") featureGRanges <- exonsBy(txDb, by="tx", use.names=TRUE)
			
			if(featureType == "CDS") featureGRanges <- cdsBy(txDb, by="tx", use.names=TRUE)
			
		} else {						
			
			if(featureType == "exon") featureGRanges <- exonsBy(txDb, by="gene")
			
			if(featureType == "CDS") featureGRanges <- cdsBy(txDb, by="gene")
			
			if(featureType == "intron") featureGRanges <- intronsByTranscript(txDb, use.names=TRUE)
			
			if(featureType == "fiveUTR") featureGRanges <- fiveUTRsByTranscript(txDb,use.names=TRUE)
			
			if(featureType == "threeUTR") featureGRanges <- threeUTRsByTranscript(txDb,use.names=TRUE)			
		}
		
		
		# ensembl uses numericals for chromosome ID (e.g. 1 rather than chr1)
		# add chr to it to be consistent with the alignment
		if(txDbName == "biomart") seqlevels(featureGRanges) <- paste("chr", seqlevels(featureGRanges), sep="")		
		
	} else {
		
		stopifnot(class(featureGRanges) %in% c("character", "GRanges", "GRangesList"))		
	}

	# assuming featureGRanges is the file name of tab-delim data
	if(!missing(featureGRanges) && is.character(featureGRanges)) {
		
		featureGRanges <- as(read.delim(featureGRanges), "GRanges")
	}
	
	
	##################### Read BAM #####################
	if(RIPSeekerRead) { # read alignment files using RIPSeeker built-in function 
		
		aligns <- combineAlignGals(bamFiles, ...)
		
	} else { # read by directly calling required functions
		
		aligns <- NULL
							
		message(sprintf("%d BAM files are combined", length(bamFiles)))
		
		for(bamFile in bamFiles) {
																
			if(paired) { # paired-end
				
				param <- ScanBamParam(flag=scanBamFlag(isUnmappedQuery=FALSE, isPaired=paired, hasUnmappedMate=FALSE, isProperPair=paired, isDuplicate=FALSE))
				
				if(is.null(aligns)) aligns <- galp2gal(readGAlignmentPairs(bamFile, param=param))
				
				if(!is.null(aligns)) aligns <- append(aligns, galp2gal(readGAlignmentPairs(bamFile, param=param)))
				
			} else { # single-end
				
				param <- ScanBamParam(flag=scanBamFlag(isUnmappedQuery=FALSE, isDuplicate=FALSE))
				
				if(is.null(aligns)) aligns <- readGAlignments(bamFile, param=param)
				
				if(!is.null(aligns)) aligns <- append(aligns, readGAlignments(bamFile, param=param))
			}						
		}
	}
	
	
	##################### Compute RPKM #####################		
	# count reads to return RangedSummarizedExperiment object
	rpkmSEobj <- summarizeOverlaps(features=featureGRanges, reads=aligns, mode=countMode, ignore.strand=ignore.strand)
		
	# get counts	
	counts <- unlist(assays(rpkmSEobj))
	
	# compute the 'K' in RPKM	
	numBases <- sum(width(featureGRanges))	
	geneLengthsInKB <- numBases / 1e3
	
	# comptue the M in RPKM	
	millionsMapped <- length(aligns) / 1e6 	# per million of total mapped reads 

	# counted reads / total reads in millions
	rpm <- counts / millionsMapped

	# reads per million per geneLength in Kb
	rpkm <- rpm / geneLengthsInKB

	# save rpkm in metadata slot of the RangedSummarizedExperiment object
	metadata(rpkmSEobj) <- list(rpkm=rpkm)
	
	
	##################### Save as data.frame #####################		
	# create a data.frame for easy viewing
	
	rpkmDF <- data.frame(count=counts, rpkm=rpkm, totalExonLength=numBases,
		      row.names=names(rowRanges(rpkmSEobj)), check.names=FALSE)
	
	if(!paired) names(rpkmDF) <- c("counts", "rpkm", "totalExonLength")
	
	if(paired) names(rpkmDF) <- c("counts", "fpkm", "totalExonLength")
	
	##################### Annotate Features #####################	
	if(moreGeneInfo) {
		
		featureID <- rownames(rpkmDF)
		
		# hack useMart to ignore unused arguments
		formals(useMart) <- c(formals(useMart), alist(... = ))
		
		mart <- useMart(...)
		
		geneInfo <- getBM(mart=mart, 
				attributes=c("chromosome_name", "start_position", "end_position", "strand",
						"external_gene_id", "ensembl_transcript_id", 
						"ensembl_gene_id", "ucsc", "description"),
				filters=idType, values = featureID)
		
		geneInfo  <- geneInfo[match(featureID, geneInfo[,idType]),]
		
		
		# form data.frame in bed format 
		rpkmDF.annotated <- cbind(geneInfo[,c(1,2,3,5)],
				totalExonLength=rpkmDF$totalExonLength,
				geneInfo[,4, drop=FALSE], uniqueMapCount=rpkmDF$counts,
				ensembl_gene_id=geneInfo$ensembl_gene_id, 
				geneLength=abs(geneInfo$end - geneInfo$start),
				description=geneInfo$description)
    
		if(!paired) rpkmDF.annotated$RPKM <- rpkmDF$rpkm
    
		if(paired) rpkmDF.annotated$FPKM <- rpkmDF$fpkm
		
		rownames(rpkmDF.annotated) <- rownames(rpkmDF)							
		
		if(!("chr" %in% rpkmDF.annotated$chromosome_name)) {
			
			rpkmDF.annotated$chromosome_name <- 
					paste("chr", rpkmDF.annotated$chromosome_name, sep="")
		}
		
		if(is.integer(rpkmDF.annotated$strand)) {
			
			rpkmDF.annotated$strand <- sub("^1$", "+", rpkmDF.annotated$strand)
			
			rpkmDF.annotated$strand <- sub("-1$", "-", rpkmDF.annotated$strand)			
		}
		
		rpkmDF <- rpkmDF.annotated
	}
			
	
	##################### Outputs #####################
	if(!missing(saveData)) save(rpkmDF, file=saveData)
	

	if(justRPKM) {
		
		ruleBasedRIPSeekResults <- rpkmSEobj
		
	} else {
		
		ruleBasedRIPSeekResults <- list(rpkmSEobj=rpkmSEobj, rpkmDF=rpkmDF, featureGRanges=featureGRanges)																							
	}
			
	
	return(ruleBasedRIPSeekResults)		
}

