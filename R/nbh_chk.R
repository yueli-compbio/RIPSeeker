# nbh_chk   Checks the parameters of a negative binomial HMM
#           and returns its dimension.
#           Use: N <- nbh_chk(TRANS,alpha,beta).
nbh_chk <- function(TRANS, alpha, beta)
{
	# Constants
	# Max deviation from 1 (beware parameters may have been computed in float)
	MAX_DEV <- 1e-6
	
	# Check input arguments
	stopifnot(!missing(TRANS))
	stopifnot(!missing(alpha))
	stopifnot(!missing(beta))


	# Check transition matrix
	N <- dim(TRANS)
	
	if(N[1] != N[2]) {
		stop("Inconsistent number in transition matrix.")
	}
	if(any(abs(apply(TRANS, 1, sum) - 1) > MAX_DEV)) {
		stop("Transition matrix is not normalized.")
	}
	if(length(alpha) != N[1]) {
		stop("Transition matrix and shape parameters must have the compatible sizes.")
	}
	if(any(alpha <= 0)) {
		stop("Shape parameters must be positive.")
	}
	if(length(beta) != N[1]) {
		stop("Transition matrix and inverse scales must have the compatible sizes.")
	}
	if(any(beta <= 0)) {
		stop("Inverse scales must be positive.")
	}
	
	return(N[1])
}