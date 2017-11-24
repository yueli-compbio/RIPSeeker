# nbm_em	Estimates the parameters of a negative binomial mixture model using EM.
# 			Use: [alpha,beta,wght, logl] =
# 			nbm_em(count,alpha_0,beta_0,wght_0,NBM_NIT_MAX,TOL) where wght, alpha
# 			and beta are the estimated model parameters, logl contains
# 			the log-likehood values for the successive iterations,


# The function is adapted from the Matlab code nbh_em.m from
# H2M/cnt Toolbox, Version 2.0
# Olivier Cappé, 31/12/97 - 22/08/2001
# ENST Dpt. TSI / LTCI (CNRS URA 820), Paris
# Further Refs.:
# - X. L. Meng, D. B. Rubin, Maximum likelihood estimation via the
#   ECM algorithm: A general framework, Biometrika, 80(2):267-278 (1993).
# - J. A. Fessler, A. O. Hero, Space-alternating generalized
#   expectation-maximization algorithm, IEEE Tr. on Signal
#   Processing, 42(10):2664 -2677 (1994).

nbm_em <- function(count, alpha, beta, wght, NBM_NIT_MAX=250, NBM_TOL=1e-2)
{
	# Inputs arguments
	stopifnot(!missing(count))
	stopifnot(!missing(alpha))
	stopifnot(!missing(beta))
	stopifnot(!missing(wght))

	
	# Data length
	Total <- length(count)	
	if(any(count < 0) || any(count != round(count))){
		stop("Data does not contain positive integers.")
	}
	
	count <- matrix(count, nrow=Total, 1)
	# Number of mixture components
	N <- nbm_chk(alpha, beta, wght)
	wght <- matrix(wght, 1, N)
	alpha <- matrix(alpha, 1, N)
	beta <- matrix(beta, 1, N)
	
	
	# Save initial alpha and beta in case error occurs in the first EM
	wght0 <- wght
	alpha0 <- alpha
	beta0 <- beta


	# Compute log(count!), the second solution is usually much faster
	# except if max(count) is very large
	cm <- max(count)
	if(cm > 50000){
		dnorm <- as.matrix(lgamma(count + 1))
	} else {
		tmp <- cumsum(rbind(0, log(as.matrix(1:max(count)))))
		dnorm <- as.matrix(tmp[count+1])
	}
	
	# Variables
	logl <- matrix(0, ncol=NBM_NIT_MAX)
	postprob <- matrix(0, Total, N)
	
		
	# Main loop of the EM algorithm
	for(nit in 1:NBM_NIT_MAX){
		
		# 1: E-Step, compute density values
		postprob <- exp( matrix(1, nrow=Total) %*% (alpha * log(beta/(1+beta)) - lgamma(alpha))
				- count %*% log(1+beta) + lgamma(count %*% matrix(1, ncol=N) + matrix(1, nrow=Total) %*% alpha)
				- dnorm %*% matrix(1, ncol=N) )
    
    # set zero value to the minimum double to avoid -inf when applying log
    # due to large dnorm (or essential large count)
		postprob <- apply(postprob, 2, function(x) {x[x==0] <- .Machine$double.xmin; x})
        
		postprob <- postprob * (matrix(1, Total, 1) %*% wght)
    																
		# Compute log-likelihood
		logl[nit] <- sum(log(apply(postprob, 1, sum)))
		
		message(sprintf('Iteration %d:\t%.3f', (nit-1), logl[nit]))
				
		postprob <- postprob / (apply(postprob, 1, sum) %*% matrix(1, ncol=N))

		
		# 4: M-Step, reestimation of the mixture weights
		wght <- apply(postprob, 2, sum)		
		wght <- wght / sum(wght)
				
								
		# 5: CM-Step 1, reestimation of the inverse scales beta with alpha fixed
		
		eq_count <- apply(postprob, 2, sum)		
		
		mu <- (t(count) %*% postprob) / eq_count
				
		beta <- alpha / mu
		
		
		# 5: CM-Step 2, reestimation of the shape parameters with beta fixed
		# Use digamma and trigamma function to perfom a Newton step on
		# the part of the intermediate quantity of EM that depends on alpha
		# Compute first derivative for all components
		# Use a Newton step for updating alpha
		# Compute first derivative		
		grad <- eq_count * (log(beta / (1+beta)) - digamma(alpha)) +
					apply(postprob * digamma(count %*% matrix(1,ncol=N) +
					matrix(1,nrow=Total) %*% alpha), 2, sum)
					
		# and second derivative
		hess <- -eq_count * trigamma(alpha) +
    			apply(postprob * trigamma(count %*% matrix(1,ncol=N) +
    			matrix(1, nrow=Total) %*% alpha), 2, sum)
    
    
		# Newton step
		tmp_step <- - grad / hess
		tmp <- alpha + tmp_step
		
		# erroneous update occurs, give up and return the previous trained parameters
		if(any(is.na(tmp))) {
			
			warning("Updated alpha becomes NA probably due to bad initial alpha or insuff. data")
			return(list(wght=wght0, alpha=alpha0, beta=beta0, 
							logl=logl, postprob=postprob))
			
		}
		
		# When performing the Newton step, one should check that the intermediate
		# quantity of EM indeed increases and that alpha does not become negative. In
		# practise this is almost never needed but the code below may help in some
		# cases (when using real bad initialization values for the parameters for
		# instance)
		while (any(tmp <= 0)){
			warning(sprintf("Alpha (%.4f) became negative! Try smaller (10%s) Newton step ...\n", tmp_step,"%"))
			tmp_step <- tmp_step/10
			tmp <- alpha + tmp_step
		}
		
		alpha <- tmp
		
		# stop iteration if improvement in logl is less than TOL (default 10^-5)
		if(nit > 1 && abs((logl[nit] - logl[nit-1])/logl[nit-1]) < NBM_TOL){
			logl <- logl[1:nit]			
			break
		}
	}

	# return the trained parameters as an data frame object
	list(wght=wght, alpha=alpha, beta=beta, 
		logl=logl, postprob=postprob)
				
}
