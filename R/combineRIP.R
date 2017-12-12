combineRIP <- function(ripPath, pattern="gff3$", combineOption="intersect", 
			pvalCutoff=1, pvalAdjCutoff=1, eFDRCutoff=1, logOddCutoff=-Inf, 
			maxgap=1e3, minIntersect, genomeBuild)
{
	
	ripFiles <- list.files(path=ripPath, pattern=pattern, recursive=TRUE, full.names=TRUE)
		
	stopifnot(length(ripFiles) > 1)
	
	stopifnot(combineOption %in% c("intersect", "merge", "union"))
				
	grlist <- lapply(ripFiles, 
	
		function(x) {
		
			gr <- as(import(x), "GRanges")
				
			         # construct cutoff rule
         			cutoff <- 
            				values(gr)$logOddScore >= logOddCutoff &
            				values(gr)$pval <= pvalCutoff &
            				values(gr)$pvalAdj <= pvalAdjCutoff &
            				values(gr)$eFDR <= eFDRCutoff
         
         			gr[cutoff]
			
		}
	)
	
	grlist <- GRangesList(grlist)
	
	if(!missing(genomeBuild)) {
		
		chromInfo <- as.data.frame(SeqinfoForUCSCGenome(genomeBuild))
		
		seqlengths(grlist) <- chromInfo[as.character(names(seqlengths(grlist))), 2]
	}
	
	names(grlist) <- basename(ripFiles)
			
	if(combineOption == "merge") return(reduce(unlist(grlist)))
	
	if(combineOption == "union") return(unlist(grlist))
		
	if(combineOption == "intersect") {
		
		if(missing(minIntersect)) minIntersect <- length(grlist) - 1
		
		stopifnot(minIntersect <= length(grlist) - 1 & minIntersect >= 1)
																
		grlist2 <- sapply(names(grlist), 
		
			function(grname) {
										
				x <- countOverlaps(grlist[[grname]], maxgap=maxgap, 
						grlist[-which(names(grlist) == grname)])

				return(grlist[[grname]][which(x >= minIntersect)])
			}
		)
    
		return(reduce(unlist(GRangesList(grlist2))))
	}								
}
