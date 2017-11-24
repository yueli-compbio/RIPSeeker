# nbh_em	Estimates the parameters of a negative binomial HMM using EM.
# 			Use: [TRANS,alpha,beta,logl,postprob,dens] =
# 			nbh_em(count,TRANS_0,alpha_0,beta_0,NIT_MAX,TOL) where TRANS, alpha
# 			and beta are the estimated model parameters, logl contains
# 			the log-likehood values for the successive iterations,
# 			postprob the marginal posterior probabilities at the last
# 			iteration and dens the negative binomial probabilities
# 			computed also at the last iteration.

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

nbh_em <- function(count, TRANS, alpha, beta, NBH_NIT_MAX=250, NBH_TOL=1e-5,
		MAXALPHA=1e7, MAXBETA=1e7)
{
	# Inputs arguments
	stopifnot(!missing(count))
	stopifnot(!missing(TRANS))
	stopifnot(!missing(alpha))
	stopifnot(!missing(beta))

	
	# Data length
	Total <- length(count)	
	if(any(count < 0) || any(count != round(count))){
		stop("Data does not contain positive integers.")
	}
	
	count <- matrix(count, nrow=Total, 1)
	# Number of mixture components
	N <- nbh_chk(TRANS, alpha, beta)
	alpha <- matrix(alpha, 1, N)
	beta <- matrix(beta, 1, N)
	
	
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
	logl <- matrix(0, ncol=NBH_NIT_MAX)
	dens <- matrix(0, Total, N)
	scale <- rbind(1, matrix(0, nrow=Total-1))
	# The forward and backward variables
	forwrd <- matrix(0, Total, N)
	bckwrd <- matrix(0, Total, N)
	
	
	# Save initial alpha and beta in case error occurs in the first EM
	alpha0 <- alpha
	beta0 <- beta
	TRANS0 <- TRANS
	postprob0 <- matrix(1, Total, N)/N
	dens0 <- dens
	logl0 <- logl
	
		
	# Main loop of the EM algorithm
	for(nit in 1:NBH_NIT_MAX) {
		
		# 1: E-Step, compute density values
		dens <- exp( matrix(1, nrow=Total) %*% (alpha * log(beta/(1+beta)) - lgamma(alpha))
				- count %*% log(1+beta) + lgamma(count %*% matrix(1, ncol=N) + matrix(1, nrow=Total) %*% alpha)
				- dnorm %*% matrix(1, ncol=N) )
    
		# set zero value to the minimum double to avoid -inf when applying log
		# due to large dnorm (or essential large count)
		dens <- apply(dens, 2, function(x) {x[x==0] <- .Machine$double.xmin; x})				
        	
				
		# 2: E-Step, forward recursion and likelihood computation
		# Use a uniform a priori probability for the initial state
		forwrd[1,] <- dens[1,]/N;
		
		for(t in 2:Total){
			forwrd[t,] <- (forwrd[t-1,] %*% TRANS) * dens[t,]
			# Systematic scaling
			scale[t] <- sum(forwrd[t,])
			forwrd[t,] <- forwrd[t,] / scale[t]
		}
				
		# Compute log-likelihood
		logl[nit] <- log(sum(forwrd[Total,])) + sum(log(scale))
		message(sprintf('Iteration %d:\t%.3f', (nit-1), logl[nit]))
		
		if(is.nan(logl[nit])) {
			
			warning("NaN logl data detected. Returning the previous training results")
			
			logl[nit] <- 0
			
			TRANS=TRANS0; alpha=alpha0; beta=beta0
			
			logl=logl0; bckwrd=postprob0; dens=dens0
			
			break
		}
		
		
		# 3: E-Step, backward recursion
		# Scale the backward variable with the forward scale factors (this ensures
		# that the reestimation of the transition matrix is correct)
		bckwrd[Total,] <- matrix(1, ncol=N)
				
		for(t in (Total-1):1) {			
			bckwrd[t,] <- (bckwrd[t+1,] * dens[t+1,]) %*% t(TRANS)
			# Apply scaling
			bckwrd[t,] <- bckwrd[t,] / scale[t]
		}
		
		
		# 4: M-Step, reestimation of the transition matrix
		# Compute unnormalized transition probabilities (this is indeed still the
		# end of the E-step, which explains that TRANS appears on the right-hand
		# side below)
		TRANS <- TRANS * (t(forwrd[1:(Total-1),]) %*% (dens[2:Total,] * bckwrd[2:Total,]))
		# Normalization of the transition matrix
		TRANS <- TRANS / (apply(TRANS, 1, sum) %*% matrix(1,ncol=N))
								
		# 5: CM-Step 1, reestimation of the inverse scales beta with alpha fixed
		# Compute a posteriori probabilities (and store them in matrix
		# bckwrd to save some space)
		bckwrd <- forwrd * bckwrd
		bckwrd <- bckwrd / ( apply(bckwrd, 1, sum) %*% matrix(1,ncol=N) )
		# Reestimate shape parameters beta conditioning on alpha
		eq_count <- apply(bckwrd, 2, sum)
		beta <- alpha / ( (t(count) %*% bckwrd) / eq_count )
		
		
		# 5: CM-Step 2, reestimation of the shape parameters with beta fixed
		# Use digamma and trigamma function to perfom a Newton step on
		# the part of the intermediate quantity of EM that depends on alpha
		# Compute first derivative for all components
		grad <- eq_count * (log(beta / (1+beta)) - digamma(alpha)) + 
					apply(bckwrd * digamma(count %*% matrix(1,ncol=N) +
					matrix(1,nrow=Total) %*% alpha), 2, sum)
				
					
		# And second derivative
		hess <- -eq_count * trigamma(alpha) +
    			apply(bckwrd * trigamma(count %*% matrix(1,ncol=N) +
    			matrix(1, nrow=Total) %*% alpha), 2, sum)
		# Newton step
		tmp_step <- - grad / hess
		tmp <- alpha + tmp_step
		
		
		############ Check any abnormal trained values ############
		# in case tmp becomes NA for bad initialization or insufficient data
		if(any(is.na(tmp))) {
			
			warning(sprintf("Updated alpha becomes NA probably %s",
							"due to bad initial alpha or insuff. data"))
			
									
			TRANS=TRANS0; alpha=alpha0; beta=beta0
			
			logl=logl0; bckwrd=postprob0; dens=dens0
			
			break
		}
		
		# in case tmp becomes NA for bad initialization or insufficient data
		if(any(alpha > MAXALPHA) || any(beta > MAXBETA)) {
			
			warning(sprintf("Updated alpha (%f) or beta (%f) becomes too large probably %s",
							alpha, beta, "due to bad initial alpha or large count"))
			
			
			TRANS=TRANS0; alpha=alpha0; beta=beta0
			
			logl=logl0; bckwrd=postprob0; dens=dens0
		
			break
		}
		############ Check END ############
		
			
				
		
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
		if(nit > 1 && abs((logl[nit] - logl[nit-1])/logl[nit-1]) < NBH_TOL){
			logl <- logl[1:nit]
			break
		}
		
		# backup trained values before M-step
		TRANS0 <- TRANS
		alpha0 <- alpha
		beta0 <- beta
		postprob0 <- bckwrd
		dens0 <- dens
		logl0 <- logl		
	}

	# return the trained parameters as an data frame object
	list(TRANS=TRANS, alpha=alpha, beta=beta, 
		logl=logl, postprob=bckwrd, dens=dens)				
}
