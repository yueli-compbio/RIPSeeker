# Function Name: 	addDummyProb
# Description: 		Create a dummy GRanges object as a placeholder in case nbh_em fails
# Input: 			y
# Output:			z
#
# Author: Yue Li
###############################################################################

addDummyProb <- function(alignGR, K = 2, randomProb = FALSE, runViterbi=FALSE, ...)
{
	x <- binCount(alignGR, ...)
		
	xcnt <- values(x)$count
	
	if(randomProb) prob <- runif(length(xcnt), 0, 1) # sample from uniform
	
	if(!randomProb) prob <- rep(1/K, length(xcnt))
	
	
	if(runViterbi) {
		
		values(x) <- cbind(xcnt, prob, 1-prob, ((1-prob) > 0.5) + 1)
		
		names(values(x)) <- c("count", paste("state", 1:K, 
					"postprob",sep="_"), "viterbi_state")
	} else {
		
		values(x) <- cbind(xcnt, prob, 1-prob)
		
		names(values(x)) <- c("count", paste("state", 1:K, 
						"postprob",sep="_"))		
	}
	
	# use the default initial paramter setting
	alpha <- tail(quantile(xcnt[xcnt>0], probs=seq(0, 1, 1/K)), K)
	
	beta <- c(1, 1)	
	
	alpha <- matrix(alpha, 1, K)
	
	beta <- matrix(beta, 1, K)
	
	TRANS <- matrix(rep(1, K*K), ncol=K)/K
	
	metadata(x) <- list(alpha = alpha, beta=beta, TRANS=TRANS)
		
	return(x)	
}
