# randindx  Generates random indexes with a specified probability distribution.
# 			I=randindx(p,Total) returns an array of T indexes distributed as specified
# 			by p (which should be a normalized probability vector). By default,
# 			Total=1.
# 			Note: Specifying a third argument (different from zero) turns all
# 			checks (dimension and normalization of p) off.
randindx <- function(p, Total, NO_CHK)
{	
	# Max deviation from 1
	MAX_DEV <- 1e-7
	# Minimum probability
	MIN_PROB <- 1e-10;

	# Default arguments
	if(missing(NO_CHK)) {NO_CHK <- 1}
	
	if(missing(Total)) {Total <- 1}
	
	# Check that p is indeed a probability vector (can be skipped)
	if(!NO_CHK) {
		# Constrain p to be a row vector
		p <- matrix(p, 1, nrow=length(p))
		# Check that p is a probability vector
		if(any(p<0) || any(p>1)) {			
			stop("inconsistent probability")
		}
		if(abs(sum(p)-1) > MAX_DEV) {
			stop("probabilities are not normalized")
		}
	}
	
	# Construct a vector which contains the inverse CDF limits
	# Taking care of the case where there are null probabilities (the inv.
	# cdf table should not contain identical values)
	if (any(p <= MIN_PROB)) {		
		ind <- which(p > MIN_PROB)
		p <- p[ind]
		NULL_FLAG <- 1		
	} else {		
		NULL_FLAG <- 0
	}
	p_p <- cumsum(p)
	p_m <- c(0, p_p[1:(length(p_p)-1)])
	
	# Generates random numbers
	R <- matrix(runif(Total), nrow=Total)
	I <- matrix(0, nrow=Total)
	
	if(Total > 1){
		for(i in 1:Total){			
			I[i] <- which( (R[i] >= p_m) & (R[i] < p_p) )			
		}				
	} else {
		# Try to do something slightly more efficient here
		I <- 1
		while(R >= p_p[I]){
			I <- I + 1			
		}				
	}
	if(NULL_FLAG) {		
		# Revert to the true indexes in case of null probabilities
		I <- ind(I)
	}
	
	I
}
