# statdis   Returns the stationary distribution of a Markov chain.
# 			Use : w = statdis(A), where A is the transition matrix.
statdis <- function(A)
{
	# Max deviation from 1 (beware parameters may have
	# been computed in float)
	MAX_DEV <- 1e-6
	
	if(missing(A)) {stop("Transition matrix A is missing.")}
		
	# Check that A is a transition matrix
	N <- dim(A)[1]
	Nc <- dim(A)[2]
	
	if(Nc != N) {
		stop("Transition matrix must be square.")
	}
	
	if(any(A<0) || any(A>1)){
		stop("Inconsistent number in transition matrix.")
	}
	
	if(any(abs(apply(A, 1, sum) - 1) > MAX_DEV)) {
		stop("Transition matrix is not normalized.")
	}
	
	# Compute left eigenvalues
	e <- eigen(t(A))
	V <- e$vector
	d <- matrix(e$values)
	n1 <- (d > 1-MAX_DEV)
	
	if(sum(n1) == 1){
		# Stationnary distribution is given by the left eigenvector corresponding
		# to the eigenvalue 1
		w <- t(V[,n1])
		w <- w/sum(w)		
	} else {
		
		if(sum(n1) > 1) {
			warning(sprintf("Warning: transition matrix is not irreducible\n"))
			w <- t(V[,n1])
			
			for(i in 1:length(w[,1])) {				
				w[i,] <- w[i,]/sum(w[i,])
			}
		} else {
			# This should never happen if A is a transition matrix
			stop("matrix does not seem to be a transition matrix")
		}		
	}
	
	return(w)	
}
