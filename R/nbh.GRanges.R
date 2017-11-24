# nbh.GRanges: inheritance function that receives an object of GRanges class 
# with additional column of read counts (for each strand) and call nbh.integer
# to derive the most probable sequence of hidden states

nbh.GRanges <- function(x, K, ...)
{							
#	nbhList <- nbh(values(binGR)$count, K, ...)
	
	binGR <- x
	
	x <- values(x)$count
	
	nbhList <- nbh(x, K, ...)	
		
	if(!is.null(nbhList$viterbi_state)) {
	
		binGRValuesNames <- c("count", paste("state", 1:K, "postprob",sep="_"), "viterbi_state")
				
		
		values(binGR) <- cbind(values(binGR)$count, nbhList$postprob, nbhList$viterbi_state)
	
	} else {
		
		binGRValuesNames <- c("count", paste("state", 1:K, "postprob",sep="_"))
		
		values(binGR) <- cbind(values(binGR)$count, nbhList$postprob)		
	}
	
	names(values(binGR)) <- binGRValuesNames		
	
	metadata(binGR) <- list(alpha = nbhList$alpha, 
			beta=nbhList$beta, TRANS=nbhList$TRANS)
	
	return(binGR)
}
