# nbm_chk   Checks the parameters of a negative binomial mixture model
#           and returns its dimension.
#           Use: N <- nbm_chk(alpha,beta,mixprop)
nbm_chk <- function(alpha, beta, wght)
{
	# Constants
	# Minimum deviation from 1 (beware parameters may have been computed in float)
	MAX_DEV <- 1e-6
	
	# Check input arguments
	stopifnot(!missing(alpha))
	stopifnot(!missing(beta))
	stopifnot(!missing(wght))
	
	# Number of mixture components
	N <- length(wght)
	
	if(any(wght < 0) | any(wght > 1)) {
		stop("Mixture weights must be between 0 and 1.")
	}
	
	if(abs(sum(wght) - 1) > MAX_DEV) {
		stop("Mixture weigths are not normalized.")		
	}
	
	if(length(alpha) != N) {
		stop("Vectors of mixture weigths and scale must have the same size.")
	}
	
	if(any(alpha <= 0)) {
		stop("Shape parameters must be positive.")
	}		
	
	if(length(beta) != N) {
		stop("Vectors of mixture weigths and intensity must have the same size.")
	}
	
	if(any(beta <= 0)) {
		stop("Inverse scales must be positive.")
	}
	
	return(N)
}