# Function Name: 	ripSeek
# Description: 		Grand master function of RIPSeeker: run mainSeek with or without control library
# Input: 			bamPath to bam files (either a single dir or the multiple file paths) 
# Output:			GAlignments and GRanges objects representating alignment object after disambiguating multihits (if any) and predicted RIP regions in a defined format (default: GFF3)
#
# Author: Yue Li
###############################################################################
ripSeek <- function(bamPath, cNAME, binSize=NULL, strandType=NULL, 		
		paired=FALSE, biomaRt_dataset, goAnno, exportFormat = "gff3", 
		annotateFormat = "txt", annotateType = "TSS", outDir, 
		padjMethod="BH", logOddCutoff=0,
		pvalCutoff=1, pvalAdjCutoff=1, eFDRCutoff=1, ...)
{
#	require(GenomicRanges)
#	require(rtracklayer)
#	require(Rsamtools)
							
	stopifnot(!missing(bamPath))
		
	if(!missing(outDir) && !file.exists(outDir)) dir.create(outDir, recursive = TRUE)	
	
	stage <- 1
			
	################ Collect BAM files ################		
	
	stopifnot(all(file.exists(bamPath)))
	
	message(sprintf("\n\n*%s. Collect alignment files\n", as.roman(stage)))
	
	# assuming in this case each bamPath specifies a bam/bed/sam file 
	if(length(bamPath) > 1) {
		
		if(any(file.info(bamPath)$isdir)) {
			
			stop("Some files are directories %s\nTry specify bamPath as either a dir or all files in that dir", bamPath[file.info(bamPath)$isdir])
		}
		
		bamFiles <- bamPath
		
	} else {
		
		if(file.info(bamPath)$isdir) {
			
			bamFiles <- list.files(bamPath, pattern = "\\.bam$", full.names=TRUE, recursive=TRUE)
			
		} else {
			
			bamFiles <- bamPath
		}
	}
	
	if(!missing(cNAME)) {
				
		bamFilesRIP <- grep(cNAME, bamFiles, invert=TRUE, value=TRUE)
		
		bamFilesCTL <- grep(cNAME, bamFiles, invert=FALSE, value=TRUE)
		
		if(length(bamFilesCTL) == 0) {
			
			warning(sprintf("The pattern %s is not found in input file names!\nAll files are treated as RIP.", cNAME))
			
			rm(bamFilesCTL)
		}
		
	} else {
				
		bamFilesRIP <- bamFiles
	}
	
	message(sprintf("\nRIP alignment files:"))
	message(sprintf("\t%s\n", bamFilesRIP))
	
	if(exists("bamFilesCTL")) {
		
		message(sprintf("\nControl alignment files:"))
		message(sprintf("\t%s\n", bamFilesCTL))
		
	} else {
		message(sprintf("No control alignment files specified"))
	}
	
	stage <- stage + 1
	
	################ mainSeek on RIP ################	
	message(sprintf("\n\n*%s. Analyzing RIP library:\n", as.roman(stage)))
		
	mainSeekOutputRIP <- mainSeek(bamFilesRIP, binSize=binSize, strandType=strandType, paired=paired, ...)
		
	################ mainSeek on CTL ################
	if(exists("bamFilesCTL")) {
				
		message(sprintf("\n\n*%s.2. Analyzing control library:\n", as.roman(stage)))
		
		# use defined binSize from RIP
		RIPBinSize <- lapply(mainSeekOutputRIP$nbhGRList, function(x) median(width(x)))
						
		mainSeekOutputCTL <- mainSeek(bamFilesCTL, binSize=RIPBinSize, strandType=strandType, paired=paired, ...)		
	}
	
	stage <- stage + 1
	
	
	################ seek RIP regions ################	
	if(exists("bamFilesCTL")) {
		# with control
		message(sprintf("\n\n*%s. Seek RIP regions with control library:\n", as.roman(stage)))
		
		RIPGRList <- endoapply(mainSeekOutputRIP$nbhGRList, 
				seekRIP, nbhGRCTL=mainSeekOutputCTL$nbhGRList,
				padjMethod=padjMethod, logOddCutoff=logOddCutoff,
				pvalCutoff=pvalCutoff, pvalAdjCutoff=pvalAdjCutoff,
				eFDRCutoff=eFDRCutoff)
		
		genome(RIPGRList) <- genome(mainSeekOutputRIP$alignGal)
				
		
	} else {
		# without control
		message(sprintf("\n\n*%s. Seek RIP regions without control library:\n", as.roman(stage)))		
		
		RIPGRList <- endoapply(mainSeekOutputRIP$nbhGRList, seekRIP,
				padjMethod=padjMethod, logOddCutoff=logOddCutoff,
				pvalCutoff=pvalCutoff, pvalAdjCutoff=pvalAdjCutoff)
		
		genome(RIPGRList) <- genome(mainSeekOutputRIP$alignGal)
	}
	
	stage <- stage + 1
	
	################ Annotate peaks ################
	
	if(!missing(biomaRt_dataset)) {
		message(sprintf("\n\n*%s. Annotate RIP regions via online ensembl database (%s):\n", 
						as.roman(stage), biomaRt_dataset))
		
		sigGRanges <- unlist(RIPGRList)
		
		names(sigGRanges) <- NULL
		
		suppressMessages(
								
				if(!missing(goAnno)) { # also run GO analysis							
					annotatedRIPGR <- annotateRIP(sigGRanges = sigGRanges,
							biomaRt_dataset = biomaRt_dataset, goAnno = goAnno, 
							featureType = annotateType, strandSpecific = !is.null(strandType),...)
				} else {
					
					annotatedRIPGR <- annotateRIP(sigGRanges = sigGRanges,
							biomaRt_dataset = biomaRt_dataset, 
							featureType = annotateType, strandSpecific = !is.null(strandType),...)					
				}
		)
		
		stage <- stage + 1	
	}
	
	if(!exists("bamFilesCTL"))  mainSeekOutputCTL <- NULL
	if(!exists("annotatedRIPGR"))  annotatedRIPGR <- NULL
	
	################ save and export results to outDir ################
	if(!missing(outDir)) {
						
		# remove backslash to avoid double backslash in the following path names
		outDir <- sub("/$", "", outDir)
		
		message(sprintf("\n\n*%s. Save and export all results to %s\n",
						as.roman(stage), outDir))
		
		outfile <- paste(outDir, "seekOutput.RData", sep="/")
		
		message(sprintf("\n**A. Saving RData to %s\n", outfile))
		
						
		################ save results in RData ################			
		save(mainSeekOutputRIP, RIPGRList, mainSeekOutputCTL, annotatedRIPGR, file=outfile)	
			
											
		################ export RIP regions ################		
		outfile <- paste(outDir, "/RIPregions.", exportFormat, sep="")
		
		message(sprintf("\n**B. Exporting %s\n", outfile))
				
		exportRIP <- unlist(RIPGRList)
				
		names(exportRIP) <- NULL
		
		exportGRanges(exportRIP, outfile, exportFormat)
		
		
		################ export annotated RIP regions ################
		if(!is.null(annotatedRIPGR)) {
						
			outfile <- paste(outDir, "/RIPregions_annotated.", exportFormat, sep="")
			
			outfile2 <- paste(outDir, "/RIPregions_annotated.txt", sep="")
			
			message(sprintf("\n\n**C. Exporting %s\n", outfile))
			
			if(is.list(annotatedRIPGR)) {
				sigGRangesAnnotated <- annotatedRIPGR$sigGRangesAnnotated
				enrichedGO <- annotatedRIPGR$enrichedGO
			} else {
				sigGRangesAnnotated <- annotatedRIPGR
				enrichedGO <- NULL
			}
			
			exportGRanges(gRanges=sigGRangesAnnotated, outfile=outfile, 
					exportFormat=exportFormat)
			
			
			exportGRanges(gRanges=sigGRangesAnnotated, outfile=outfile2, 
					exportFormat="txt")
			
			
			
			################ export enriched GO ################
			if(!is.null(enrichedGO)) {
				outfile <- paste(outDir, "/RIPregions_enrichedGO.txt", sep="")
				
				message(sprintf("\n\n**D. Exporting %s\n", outfile))
				
				write.table(enrichedGO, file=outfile, row.names=F, quote=F, sep="\t")
			}
			
		}
		
		stage <- stage + 1
	}
	
		
	################ return results in list ################
	message(sprintf("\n\n*%s. Return results list to R console\n", as.roman(stage)))
		
	seekOut <- list(mainSeekOutputRIP=mainSeekOutputRIP,
			mainSeekOutputCTL=mainSeekOutputCTL, RIPGRList=RIPGRList,
			annotatedRIPGR=annotatedRIPGR)				
	
	
	return(seekOut)
					
}
