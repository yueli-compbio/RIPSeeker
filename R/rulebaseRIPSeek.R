rulebaseRIPSeek <- function(bamFiles, cNAME, featureGRanges, 
	rpkmCutoff=0.4, fcCutoff=3, moreRIPGeneInfo=TRUE, 
	idType="ensembl_transcript_id", myMin=.Machine$double.xmin,
	saveRData, ...)
{	
	##################### Collect BAM files #####################
	bamFilesRIP <- grep(cNAME, bamFiles, invert=TRUE, value=TRUE)
	
	bamFilesCTL <- grep(cNAME, bamFiles, invert=FALSE, value=TRUE)
	
		
	##################### Compute RPKM #####################	
	if(!missing(featureGRanges)) {
		
		nRPKM <- computeRPKM(bamFiles=bamFilesRIP, justRPKM=TRUE, featureGRanges=featureGRanges, ...)
		
		dRPKM <- computeRPKM(bamFiles=bamFilesCTL,featureGRanges=featureGRanges, ...)				
		
	} else {
	
		nRPKMList <- computeRPKM(bamFiles=bamFilesRIP, justRPKM=FALSE, ...)
				
		nRPKM <- nRPKMList$rpkmSEobj
		
		featureGRanges <- nRPKMList$featureGRanges
				
		dRPKM <- computeRPKM(bamFiles=bamFilesCTL, featureGRanges=featureGRanges,...)		
	}
		
			
	##################### Compute Foldchange #####################
	RIP_rpkm <- metadata(nRPKM)[["rpkm"]]
	CTL_rpkm <- metadata(dRPKM)[["rpkm"]]
	foldchange <- (RIP_rpkm + myMin) / (CTL_rpkm + myMin)
		
	
	##################### Rule-based Threshold #####################		
	# create a data.frame for easy viewing
	rpkmDF <- data.frame(RIP_count=unlist(assays(nRPKM)), 
			     RIP_rpkm=RIP_rpkm,
			     CTL_count=unlist(assays(dRPKM)), 
			     CTL_rpkm=CTL_rpkm,
			     foldchange=foldchange,
			     row.names=names(rowRanges(nRPKM)), check.names=FALSE)
						
				
	rule <- (rpkmDF$RIP_rpkm >= rpkmCutoff) * (rpkmDF$foldchange >= fcCutoff)
				
	rpkmDF.passed <- rpkmDF[which(rule == 1), ]
	
			
	##################### Annotate Features #####################	
	if(moreRIPGeneInfo) {
						
		featureID <- rownames(rpkmDF.passed)
		
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
		rpkmDF.annotated <- cbind(	geneInfo[,c(1,2,3)], rpkmDF.passed[, c(2,5)],
									geneInfo[,4, drop=FALSE], rpkmDF.passed[, -c(2,5)],
									geneInfo[, -c(1:4)]	)	
							
		rownames(rpkmDF.annotated) <- rownames(rpkmDF.passed)							
		
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
	
	
	
	##################### Output #####################
	ruleBasedRIPSeekResults <- list(nRPKM=nRPKM, dRPKM=dRPKM, rpkmDF=rpkmDF, 
						rpkmCutoff=rpkmCutoff, fcCutoff=fcCutoff,
						featureGRanges=featureGRanges)
						
		
	if(!missing(saveRData)) {
		
		save(ruleBasedRIPSeekResults, file=saveRData)
				
		write.table(rpkmDF.annotated, file=sub("RData$", "txt", saveRData), sep="\t", quote=F, row.names=F)		
	}				
		
	return(ruleBasedRIPSeekResults)			
}
