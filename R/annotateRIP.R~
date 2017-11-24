# Function Name: 	annotateRIP
# Description: 		x
# Input: 			y
# Output:			z
#
# Author: Yue Li
###############################################################################

annotateRIP <- function(sigGRanges, biomaRt_dataset, featureType="TSS", goAnno,
					strandSpecific=FALSE, exportFormat = "txt",
					hasGOdb=!missing(goAnno), goPval=0.1, outDir, ...)
{
	stopifnot(!missing(sigGRanges))
	
	stopifnot(!missing(biomaRt_dataset))
	
	stopifnot(require(biomaRt))
	
	stopifnot(require(ChIPpeakAnno))
	
	# hack useMart to ignore unused arguments
	formals(useMart) <- c(formals(useMart), alist(... = ))	
	
	mart <- useMart(dataset=biomaRt_dataset, ...)			
	
	anno <- getAnnotation(mart, featureType=featureType)
	
	message(sprintf("\n\n*** Annotating %d genomic ranges with %s.\n",
			length(sigGRanges), biomaRt_dataset))

	names(sigGRanges) <- NULL
	
	
	# hack annotatePeakInBatch to ignore unused arguments
	formals(annotatePeakInBatch) <- c(formals(annotatePeakInBatch), alist(... = ))
	
	
	annotatedPeak <- annotatePeakInBatch(RangedData(sigGRanges), 
		AnnotationData = anno, output="both", ...)
			

	if(strandSpecific) subsetByOverlaps(annotatedPeak, sigGRanges, ignore.strand=FALSE)
			
	
	# more useful information based on ensembl gene ID
	geneInfo <- getBM(mart=mart, attributes=c("ensembl_gene_id", "external_gene_name", "description"),filters="ensembl_gene_id", values = annotatedPeak$feature)
			
	geneInfo  <- geneInfo[match(annotatedPeak$feature, geneInfo$ensembl_gene_id),]
		
	
	sigGRangesIdx <- as.numeric(annotatedPeak$peak)
	
		
	sigGRangesAnnotated <- sigGRanges[sigGRangesIdx]
		
	
	combinedInfo <- cbind(geneInfo, as.data.frame(values(sigGRangesAnnotated)),
								as.data.frame(values(annotatedPeak)) )
	
	# remove useless information: "space" (chr)
	combinedInfo <- combinedInfo[, grep("space", colnames(combinedInfo), invert=TRUE)]
						
								
	colnames(combinedInfo)[colnames(combinedInfo) == "strand"] <- "feature_strand"
															
	values(sigGRangesAnnotated) <- combinedInfo
	
	sigGRangesAnnotated <- sigGRangesAnnotated[order(values(sigGRangesAnnotated)$peak)]
	
		
	if(!missing(goAnno)) {
						
		hasGOdb <- require(package=goAnno, character.only=TRUE)
		
		if(hasGOdb) {
			
			message(sprintf("\n\n*** GO analysis for %d associated features with %s.\n",
							length(annotatedPeak$feature), goAnno))
			
			# hack functions to ignore unused arguments
			formals(getEnrichedGO) <- c(formals(getEnrichedGO), alist(... = ))
			
			enrichedGO <- getEnrichedGO(annotatedPeak, orgAnn=goAnno, maxP = goPval,
					multiAdj = TRUE, multiAdjMethod="BH", ...)
			
			# remove redundant sets
			enrichedGO$bp <- unique( enrichedGO$bp[,-11] )
			enrichedGO$mf <- unique( enrichedGO$mf[,-11] )
			enrichedGO$cc <- unique( enrichedGO$cc[,-11] )
			
			enrichedGO <- rbind(enrichedGO$bp[order(enrichedGO$bp$pvalue), ],
					enrichedGO$mf[order(enrichedGO$mf$pvalue), ],
					enrichedGO$cc[order(enrichedGO$cc$pvalue), ])
		} else {
			
			warning(sprintf("%s is not found!", goAnno))
		}
	}
	
	
	
	################ save and export results to outDir ################		
	if(!missing(outDir)) {
		
		# remove backslash to avoid double backslash in the following path names
		outDir <- sub("/$", "", outDir)
		
		
		outfile <- paste(outDir, "/RIPGRanges_annotated.RData", sep="")
		
		message(sprintf("\n\n*** Saving RData to %s\n", outfile))
		
		
		################ save results in RData ################
		
		save(sigGRangesAnnotated, file=outfile)
		
		
		
		################ export RIP regions ################
		outfile <- paste(outDir, "/RIPregions_annotated.", exportFormat, sep="")
		
		message(sprintf("\n\n*** Exporting %s\n", outfile))
		
		exportGRanges(gRanges=sigGRangesAnnotated, outfile=outfile, 
				exportFormat=exportFormat)
		
		
		
		################ export enriched GO ################
		if(hasGOdb) {
			outfile <- paste(outDir, "/RIPregions_enrichedGO.txt", sep="")
			
			message(sprintf("\n\n*** Exporting %s\n", outfile))
			
			write.table(enrichedGO, file=outfile, row.names=F, quote=F, sep="\t")						
		}
	}
	
	if(hasGOdb) return(list(sigGRangesAnnotated=sigGRangesAnnotated, enrichedGO=enrichedGO))
	
	if(!hasGOdb) return(sigGRangesAnnotated)
	
}